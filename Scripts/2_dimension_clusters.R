# dimension_clusters.R
# Analysis pipeline for phenotype clustering.

# ------------------------------------------------------------------------------
# Data Loading and Preprocessing
# ------------------------------------------------------------------------------

load_dx_wide <- function(path, include_nodiag = FALSE, min_prev = 0, max_prev = 1) {
  # Load raw binary matrix and sanitize column names
  DX0 <- suppressWarnings(readr::read_csv2(path, show_col_types = FALSE))
  stopifnot(is.data.frame(DX0), ncol(DX0) >= 2)
  names(DX0) <- trimws(gsub("\ufeff", "", names(DX0)))
  
  if (!identical(names(DX0)[1], "participant_id")) {
    stop('First column must be "participant_id" (exact).')
  }
  stopifnot(
    !anyDuplicated(DX0$participant_id),
    all(nzchar(as.character(DX0$participant_id)))
  )
  
  # Enforce binary encoding (0/1) for diagnosis columns
  dx_cols <- setdiff(names(DX0), "participant_id")
  for (nm in dx_cols) {
    v <- suppressWarnings(as.numeric(DX0[[nm]]))
    if (!all(v %in% c(0, 1, NA))) stop(sprintf("Non-binary column: %s", nm))
    DX0[[nm]] <- as.integer(v)
  }
  
  # Compute summary flags (ANY_DX, NODIAG)
  any_dx <- as.integer(rowSums(DX0[, dx_cols] == 1L, na.rm = TRUE) > 0L)
  if (!"ANY_DX" %in% names(DX0)) DX0$ANY_DX <- any_dx
  if (include_nodiag && !"NODIAG" %in% names(DX0)) DX0$NODIAG <- 1L - any_dx
  
  # Filter diagnoses based on prevalence thresholds
  n <- nrow(DX0)
  pos <- vapply(dx_cols, function(nm) sum(DX0[[nm]] == 1L, na.rm = TRUE), integer(1))
  prev <- pos / n
  keep <- names(prev)[prev >= min_prev & prev <= max_prev]
  keep <- union(keep, intersect(c("ANY_DX", "NODIAG"), names(DX0)))
  
  DX1 <- DX0[, c("participant_id", keep), drop = FALSE]
  
  # Optional exclusion of Not Otherwise Specified (NOS) codes
  if (DX_DENY_NOS) {
    kk <- setdiff(names(DX1), "participant_id")
    kk <- kk[!grepl("\\bNOS\\b$", kk, ignore.case = TRUE)]
    DX1 <- DX1[, c("participant_id", kk), drop = FALSE]
  }
  DX1
}

# ------------------------------------------------------------------------------
# Deduplication
# ------------------------------------------------------------------------------

dedup_dx <- function(DX) {
  # Collapse participants with identical diagnosis signatures into weighted unique rows
  # to reduce computational complexity for distance matrix calculation.
  cols <- setdiff(names(DX), "participant_id")
  X <- as.data.frame(DX[, cols, drop = FALSE])
  for (nm in cols) {
    v <- X[[nm]]
    v[is.na(v)] <- 0L
    X[[nm]] <- as.integer(v)
  }
  
  sig <- apply(X, 1, paste0, collapse = "")
  tab <- as.data.frame(table(sig), stringsAsFactors = FALSE)
  first <- tapply(seq_len(nrow(DX)), sig, `[`, 1)
  
  DXu <- DX[unname(first), , drop = FALSE]
  DXu$mult <- as.integer(tab$Freq[match(names(first), tab$sig)])
  
  list(DXu = DXu, sig_all = sig)
}

build_counts <- function(DX) {
  # Compute global prevalence statistics
  cols <- setdiff(names(DX), "participant_id")
  n1 <- vapply(cols, function(v) sum(DX[[v]] == 1L, na.rm = TRUE), integer(1))
  n0 <- nrow(DX) - n1
  data.frame(
    dx = cols, n1 = as.integer(n1), n0 = as.integer(n0),
    prev = ifelse(n1 + n0 > 0, n1 / (n1 + n0), NA_real_), stringsAsFactors = FALSE
  )
}

# ------------------------------------------------------------------------------
# Distance Metrics and Graph Construction
# ------------------------------------------------------------------------------

gower_dist_dx <- function(DXdf) {
  # Compute Gower distance. For asymmetric binary data, this is equivalent to Jaccard distance.
  cols <- setdiff(names(DXdf), c("participant_id", "mult", "ANY_DX"))
  if (length(cols) < 2) return(NULL)
  
  # Remove invariant columns
  vary <- vapply(cols, function(nm) {
    v <- DXdf[[nm]]
    u <- unique(v[is.finite(v)])
    length(u) >= 2 && all(u %in% c(0L, 1L))
  }, logical(1))
  cols <- cols[vary]
  
  if (length(cols) < 2) return(NULL)
  
  suppressWarnings(cluster::daisy(DXdf[, cols, drop = FALSE],
                                  metric = "gower",
                                  type = list(asymm = cols), stand = FALSE
  ))
}

knn_graph_from_dist <- function(D, k = 9, variant = c("mutual", "union"),
                                local_scale = TRUE, mult = NULL,
                                kernel = c("gaussian2", "laplace")) {
  # Construct a k-Nearest Neighbor graph using adaptive bandwidth scaling.
  # Handles multiplicity (deduplicated weights) during edge weighting.
  variant <- match.arg(variant)
  kernel <- match.arg(kernel)
  M <- as.matrix(D)
  n <- nrow(M)
  diag(M) <- Inf
  
  # Identify k-th nearest neighbor distance for adaptive scaling
  kth <- function(r, k) {
    rf <- r[is.finite(r)]
    if (!length(rf)) return(Inf)
    sort(rf, partial = k)[k]
  }
  kd <- apply(M, 1, kth, k = k)
  kd[!is.finite(kd) | kd <= 0] <- median(M[is.finite(M)])
  
  nbr_idx <- lapply(seq_len(n), function(i) {
    di <- M[i, ]
    ok <- which(is.finite(di))
    ok[order(di[ok])[seq_len(min(k, length(ok)))]]
  })
  
  edges <- list()
  is_union <- variant == "union"
  
  for (i in seq_len(n)) {
    js <- nbr_idx[[i]]
    partners <- if (is_union) {
      unique(c(js, which(vapply(nbr_idx, function(x) i %in% x, FALSE))))
    } else {
      intersect(js, which(vapply(nbr_idx, function(x) i %in% x, FALSE)))
    }
    partners <- partners[partners > i]
    
    if (!length(partners)) next
    
    # Compute similarity kernel
    if (kernel == "gaussian2") {
      sprod <- pmax(1e-12, kd[i] * kd[partners])
      w <- exp(-(M[i, partners]^2) / (2 * sprod))
    } else {
      sprod <- sqrt(pmax(1e-12, kd[i] * kd[partners]))
      w <- exp(-M[i, partners] / sprod)
    }
    
    # Adjust weights for multiplicity of unique signatures
    if (!is.null(mult)) {
      mi <- ifelse(is.finite(mult[i]) & mult[i] > 0, mult[i], 1)
      mj <- ifelse(is.finite(mult[partners]) & mult[partners] > 0, mult[partners], 1)
      w <- w * sqrt(pmax(1, mi) * pmax(1, mj))
    }
    
    edges[[length(edges) + 1L]] <- data.frame(
      from = rep(i, length(partners)),
      to = partners,
      weight = w
    )
  }
    
  if (!is.null(mult)) {
    w_self <- ifelse(mult > 1, mult, 0)
  } else {
    w_self <- rep(0, n)
  }
  
  # Only add self-loops if they have weight > 0
  if (any(w_self > 0)) {
    edges[[length(edges) + 1L]] <- data.frame(
      from = seq_len(n), 
      to = seq_len(n), 
      weight = w_self
    )
  }
  
  # Combine all edges
  if (length(edges) == 0) {
    # Fallback for completely empty graph to avoid crash
    df_edges <- data.frame(from=integer(), to=integer(), weight=numeric())
    warning("No edges were produced, empty graph.")
  } else {
    df_edges <- do.call(rbind, edges)
  }
  
  g <- igraph::graph_from_data_frame(do.call(rbind, edges), directed = FALSE, vertices = data.frame(name = seq_len(n)))
  igraph::simplify(g, remove.multiple = TRUE, remove.loops = FALSE, edge.attr.comb = list(weight = "sum"))
}

# ------------------------------------------------------------------------------
# Community Detection
# ------------------------------------------------------------------------------

community_labels <- function(g,
                             algo = COMMUNITY_ALGO,
                             objective = LEIDEN_OBJECTIVE,
                             gamma = LEIDEN_GAMMA,
                             iters = LEIDEN_ITERS) {
  # Wrapper for Leiden or Louvain community detection algorithms
  a <- tolower(algo)
  if (a == "leiden") {
    if (!"cluster_leiden" %in% getNamespaceExports("igraph")) {
      stop("Leiden not available. Install latest igraph or use algo='louvain'.")
    }
    obj <- if (tolower(objective) == "modularity") "modularity" else "CPM"
    fn <- igraph::cluster_leiden
    
    # Handle version discrepancies in igraph argument names
    if ("resolution" %in% names(formals(fn))) {
      memb <- fn(g,
                 weights = igraph::E(g)$weight, objective_function = obj,
                 resolution = gamma, n_iterations = as.integer(iters)
      )
    } else {
      memb <- fn(g,
                 weights = igraph::E(g)$weight, objective_function = obj,
                 resolution_parameter = gamma, n_iterations = as.integer(iters)
      )
    }
    return(igraph::membership(memb))
  }
  igraph::membership(igraph::cluster_louvain(g, weights = igraph::E(g)$weight))
}

# ------------------------------------------------------------------------------
# Validation Metrics (Null Models & ID)
# ------------------------------------------------------------------------------

degree_null_Q <- function(g, memb, reps = 200L) {
  # Generate null distribution of modularity via degree-preserving rewiring
  Q <- numeric(reps)
  w <- igraph::E(g)$weight
  m <- igraph::gsize(g)
  for (b in seq_len(reps)) {
    g2 <- igraph::rewire(g, with = igraph::keeping_degseq(niter = 100 * max(1L, m)))
    igraph::E(g2)$weight <- sample(w, length(w), replace = FALSE)
    Q[b] <- igraph::modularity(g2, memb, weights = igraph::E(g2)$weight)
  }
  Q
}

silhouette_graph <- function(g, memb) {
  # Compute silhouette width using graph geodesic distance (inverse weight)
  g2 <- igraph::set_edge_attr(g, "length", value = 1 / pmax(igraph::E(g)$weight, 1e-9))
  Dg <- igraph::distances(g2, weights = igraph::E(g2)$length)
  cluster::silhouette(memb, dist(Dg))
}

bootstrap_ari <- function(Dm, k, variant, local_scale, algo, mult, B = 200) {
  # Assess clustering stability via bootstrap resampling and Adjusted Rand Index
  Dm <- as.matrix(Dm)
  n <- nrow(Dm)
  diag(Dm) <- Inf
  p <- as.numeric(mult)
  p[!is.finite(p) | p <= 0] <- 1
  p <- p / sum(p)
  
  g0 <- knn_graph_from_dist(Dm, k = k, variant = variant, local_scale = local_scale, mult = mult)
  labs0 <- as.integer(community_labels(g0, algo = algo))
  
  ari <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    idx <- sort(unique(sample.int(n, n, replace = TRUE, prob = p)))
    Db <- Dm[idx, idx, drop = FALSE]
    
    g <- try(knn_graph_from_dist(Db, k = k, variant = variant, local_scale = local_scale, mult = mult[idx]),
             silent = TRUE
    )
    if (inherits(g, "try-error") || igraph::ecount(g) == 0L) next
    labs <- as.integer(community_labels(g, algo = algo))
    
    # Map bootstrapped samples back to original space via 1-NN for ARI calculation
    nn_idx <- apply(Dm[, idx, drop = FALSE], 1L, function(r) {
      j <- which.min(r)
      if (length(j) && is.finite(r[j])) idx[j] else idx[1]
    })
    labs_back <- labs[match(nn_idx, idx)]
    
    a <- as.integer(labs0)
    bvec <- as.integer(labs_back)
    ok <- is.finite(a) & is.finite(bvec)
    if (length(unique(a[ok])) >= 2L && length(unique(bvec[ok])) >= 2L) {
      ari[sum(!is.na(ari)) + 1L] <- aricode::ARI(a[ok], bvec[ok])
    }
  }
  ari
}

# ------------------------------------------------------------------------------
# Intrinsic Dimension Estimation
# ------------------------------------------------------------------------------

twonn_id_from_dist <- function(D, eps = 1e-8, trim = 0.02) {
  # TWO-NN estimator (Facco et al.)
  M <- as.matrix(D)
  n <- nrow(M)
  if (n < 3) return(NA_real_)
  diag(M) <- Inf
  
  r <- vapply(seq_len(n), function(i) {
    di <- M[i, ]
    di <- di[is.finite(di)]
    if (length(di) < 2) return(NA_real_)
    ds <- sort(di, partial = 2)[1:2]
    d1 <- max(ds[1], eps)
    d2 <- max(ds[2], d1 + eps)
    d2 / d1
  }, numeric(1))
  
  r <- r[is.finite(r) & r > 1]
  if (!length(r)) return(NA_real_)
  
  logr <- sort(log(r))
  k <- floor(trim * length(logr))
  if (k > 0 && 2 * k < length(logr)) logr <- logr[(k + 1):(length(logr) - k)]
  1 / mean(logr)
}

lb_mle_id <- function(Dm, k_lo = 5, k_hi = 15) {
  # Levina-Bickel Maximum Likelihood Estimator
  Dm <- as.matrix(Dm)
  n <- nrow(Dm)
  diag(Dm) <- Inf
  if (n <= k_lo) return(NA_real_)
  
  k_hi <- max(k_lo, min(k_hi, n - 1))
  ids <- sapply(k_lo:k_hi, function(k) {
    nn <- t(apply(Dm, 1L, function(r) {
      rf <- r[is.finite(r)]
      if (length(rf) < k) return(rep(NA_real_, k))
      sort(rf, partial = k)[1:k]
    }))
    if (!nrow(nn)) return(NA_real_)
    l <- log(nn[, k, drop = TRUE] / nn[, 1:(k - 1), drop = FALSE])
    d <- 1 / rowMeans(l, na.rm = TRUE)
    mean(d[is.finite(d)], na.rm = TRUE)
  })
  mean(ids, na.rm = TRUE)
}

core_band_idx <- function(D, k = 10, band = c(0.15, 0.85)) {
  # Identify core points based on centrality (k-th neighbor distance)
  M <- as.matrix(D)
  diag(M) <- Inf
  kth <- function(r) {
    rf <- r[is.finite(r)]
    if (!length(rf)) return(NA_real_)
    k_eff <- min(k, length(rf))
    sort(rf, partial = k_eff)[k_eff]
  }
  rk <- apply(M, 1, kth)
  ok <- is.finite(rk)
  if (!any(ok)) return(integer(0))
  q <- stats::quantile(rk[ok], band, na.rm = TRUE)
  which(ok & rk >= q[1] & rk <= q[2])
}

# ------------------------------------------------------------------------------
# Characterization Pillars (A/B/C)
# ------------------------------------------------------------------------------

# Pillar A: Enrichment (Guards + BH FDR)
pillar_A <- function(DX, clusters, alpha_fdr = 0.05, min_prev_in = 0, min_or = 2.0,
                     min_in = MIN_CASES_IN, min_out = MIN_CASES_OUT, min_total = MIN_CASES_TOTAL) {
  # Identify enriched diagnoses using Fisher's Exact Test with global BH correction
  stopifnot(all(c("participant_id", "cluster") %in% names(clusters)))
  cols <- setdiff(names(DX), c("participant_id", "ANY_DX"))
  df <- merge(DX, clusters, by = "participant_id", all.x = TRUE)
  df$cluster[is.na(df$cluster)] <- 0L
  cl_list <- setdiff(sort(unique(df$cluster)), 0L)
  
  out <- list()
  for (cl in cl_list) {
    in_c <- df$cluster == cl
    for (v in cols) {
      a <- sum(df[[v]][in_c] == 1L, na.rm = TRUE)
      b <- sum(df[[v]][!in_c] == 1L, na.rm = TRUE)
      c0 <- sum(df[[v]][in_c] == 0L, na.rm = TRUE)
      d <- sum(df[[v]][!in_c] == 0L, na.rm = TRUE)
      
      if ((a + b) < min_total || a < min_in || b < min_out) next
      
      prev_in <- a / max(1, a + c0)
      prev_all <- (a + b) / max(1, a + b + c0 + d)
      OR <- if (b == 0 || c0 == 0) Inf else (a * d) / (b * c0)
      p <- tryCatch(fisher.test(matrix(c(a, b, c0, d), 2, byrow = TRUE), alternative = "greater")$p.value,
                    error = function(e) NA_real_
      )
      out[[length(out) + 1L]] <- data.frame(
        cluster = cl, dx = v, prev_in = prev_in, prev_all = prev_all,
        OR = OR, p = p, stringsAsFactors = FALSE
      )
    }
  }
  
  if (!length(out)) return(list(enrichment = NULL, majors = character(0)))
  
  res <- dplyr::bind_rows(out)
  res$FDR <- p.adjust(res$p, "BH") # Global correction
  majors <- unique(res$dx[res$FDR <= alpha_fdr & res$prev_in >= min_prev_in & res$OR >= min_or])
  list(enrichment = res[order(res$cluster, res$FDR, -res$OR), ], majors = majors)
}

# Pillar B: Localization (Assortativity + kNN Purity)
pillar_B <- function(g, DXu_id, B = 1000, n_pos_min = 10, n_neg_min = 10) {
  # Calculate topological localization of diagnoses on the graph vs label permutation nulls
  el <- if (igraph::ecount(g) > 0L) igraph::as_edgelist(g, names = FALSE) else matrix(integer(), ncol = 2)
  nbrs <- igraph::adjacent_vertices(g, igraph::V(g), mode = "all")
  
  frac_pos_neigh <- function(i, lab) {
    nb <- as.integer(nbrs[[i]])
    nb <- setdiff(nb, i)
    if (!length(nb)) return(NA_real_)
    mean(lab[nb], na.rm = TRUE)
  }
  
  cols <- setdiff(names(DXu_id), c("participant_id", "mult", "ANY_DX"))
  out <- lapply(cols, function(v) {
    pos_tot <- sum(DXu_id[[v]] == 1L, na.rm = TRUE)
    neg_tot <- sum(DXu_id[[v]] == 0L, na.rm = TRUE)
    if (pos_tot < n_pos_min || neg_tot < n_neg_min) return(NULL)
    
    z <- as.integer(DXu_id[[v]] == 1L)
    z[is.na(z)] <- 0L
    
    # Assortativity coefficient
    if (nrow(el)) {
      zi <- z[el[, 1]]; zj <- z[el[, 2]]
      tab <- matrix(0, 2, 2)
      for (i in seq_len(nrow(el))) {
        a <- zi[i] + 1L; b <- zj[i] + 1L
        tab[a, b] <- tab[a, b] + 1L; tab[b, a] <- tab[b, a] + 1L
      }
      e <- tab / sum(tab)
      a <- rowSums(e); b <- colSums(e)
      denom <- 1 - sum(a * b)
      r_obs <- if (abs(denom) < 1e-12) NA_real_ else (sum(diag(e)) - sum(a * b)) / denom
      
      # Permutation null for assortativity
      r_null <- replicate(B, {
        zsh <- sample(z)
        zsi <- zsh[el[, 1]]; zsj <- zsh[el[, 2]]
        tb0 <- matrix(0, 2, 2)
        for (i in seq_len(nrow(el))) {
          aa <- zsi[i] + 1L; bb <- zsj[i] + 1L
          tb0[aa, bb] <- tb0[aa, bb] + 1L; tb0[bb, aa] <- tb0[bb, aa] + 1L
        }
        e0 <- tb0 / sum(tb0)
        a0 <- rowSums(e0); b0 <- colSums(e0)
        d0 <- 1 - sum(a0 * b0)
        if (abs(d0) < 1e-12) NA_real_ else (sum(diag(e0)) - sum(a0 * b0)) / d0
      })
      r_null <- r_null[is.finite(r_null)]
      p_r <- if (!is.finite(r_obs) || !length(r_null)) NA_real_ else (sum(r_null >= r_obs) + 1) / (length(r_null) + 1)
    } else {
      r_obs <- NA_real_; p_r <- NA_real_
    }
    
    # Local Purity
    pos_idx <- which(z == 1L)
    pur_obs <- if (length(pos_idx)) mean(vapply(pos_idx, frac_pos_neigh, numeric(1), lab = z), na.rm = TRUE) else NA_real_
    
    pur_null <- replicate(B, {
      zsh <- sample(z)
      pos_sh <- which(zsh == 1L)
      if (!length(pos_sh)) return(NA_real_)
      mean(vapply(pos_sh, frac_pos_neigh, numeric(1), lab = zsh), na.rm = TRUE)
    })
    pur_null <- pur_null[is.finite(pur_null)]
    p_pur <- if (!is.finite(pur_obs) || !length(pur_null)) NA_real_ else (sum(pur_null >= pur_obs) + 1) / (length(pur_null) + 1)
    
    data.frame(dx = v, assort = r_obs, assort_p = p_r, purity = pur_obs, purity_p = p_pur, stringsAsFactors = FALSE)
  })
  dplyr::bind_rows(Filter(Negate(is.null), out))
}

# Pillar C: Predictability (AUC on kNN score)
auc_one_vs_rest_knn_weighted <- function(DXu_id, target, k = 10, pos_min = 10, neg_min = 10) {
  # Compute One-vs-Rest AUC using graph transition probabilities
  cols <- setdiff(names(DXu_id), c("participant_id", "mult", "ANY_DX", target))
  if (length(cols) < 2) return(NA_real_)
  
  D <- cluster::daisy(DXu_id[, cols, drop = FALSE], metric = "gower", type = list(asymm = cols), stand = FALSE)
  g <- knn_graph_from_dist(D, k = k, variant = "mutual", local_scale = TRUE, mult = DXu_id$mult)
  A <- igraph::as_adj(g, type = "both", attr = "weight", sparse = FALSE)
  P <- A / pmax(rowSums(A), 1e-9) # Transition matrix
  
  y <- as.integer(DXu_id[[target]] == 1L)
  y[is.na(y)] <- 0L
  w <- as.numeric(DXu_id$mult)
  w[!is.finite(w) | w <= 0] <- 1
  
  if (sum(w[y == 1]) < pos_min || sum(w[y == 0]) < neg_min) return(NA_real_)
  
  s <- as.numeric(P %*% y) # Score = smooth of labels
  auc_weighted(y, s, w)
}

pillar_C <- function(DXu_id, k = 10, pos_min = 10, neg_min = 10) {
  cols <- setdiff(names(DXu_id), c("participant_id", "mult", "ANY_DX"))
  out <- lapply(cols, function(v) {
    data.frame(dx = v, auc = auc_one_vs_rest_knn_weighted(DXu_id, v, k, pos_min, neg_min))
  })
  dplyr::bind_rows(Filter(Negate(is.null), out))
}

# ------------------------------------------------------------------------------
# Heuristic for Unclustered Nodes
# ------------------------------------------------------------------------------

absorb_leftovers <- function(g, D, memb, DXu_id, mult,
                                        tau = 0.70, k_min = 3, shared_min = 1,
                                        B = 200, alpha = 0.05, seed = 42) {
  # Conservatively assign unclustered nodes to communities if:
  # 1. Neighbor purity > tau
  # 2. Distance is within core distribution of target community
  # 3. Permutation test of neighborhood weight is significant
  
  set.seed(seed)
  kept <- sort(setdiff(unique(memb), 0L))
  if (!length(kept)) return(list(membership = memb, absorbed = integer(0), q = numeric(0)))
  
  M <- as.matrix(D)
  diag(M) <- NA_real_
  
  # Prepare neighborhood weights
  get_w <- function(i) {
    inc <- igraph::incident(g, i, mode = "all")
    if (!length(inc)) return(list(j = integer(0), w = numeric(0)))
    ed <- igraph::ends(g, inc, names = FALSE)
    js <- ifelse(ed[, 1] == i, ed[, 2], ed[, 1])
    list(j = js, w = igraph::E(g)$weight[inc])
  }
  
  cols <- setdiff(names(DXu_id), c("participant_id", "mult"))
  X <- as.matrix(data.frame(lapply(DXu_id[, cols, drop = FALSE], function(v) {
    v <- as.integer(v == 1L); v[is.na(v)] <- 0L; v
  })))
  
  w_mult <- as.numeric(mult)
  w_mult[!is.finite(w_mult) | w_mult <= 0] <- 1
  cand <- which(memb == 0L)
  if (!length(cand)) return(list(membership = memb, absorbed = integer(0), q = numeric(0)))
  
  assign_to <- rep(0L, length(cand)); names(assign_to) <- as.character(cand)
  pvals <- rep(NA_real_, length(cand))
  kept_idx <- which(memb > 0L)
  kept_labels <- memb[kept_idx]
  
  for (ii in seq_along(cand)) {
    i <- cand[ii]
    nb <- get_w(i)
    if (!length(nb$j)) next
    
    lab_nb <- memb[nb$j]
    ok <- which(lab_nb %in% kept)
    if (length(ok) < k_min) next
    
    j_ok <- nb$j[ok]
    w_ok <- nb$w[ok] * w_mult[j_ok]
    
    supp <- tapply(w_ok, lab_nb[ok], sum)
    supp[setdiff(as.character(kept), names(supp))] <- 0
    supp <- supp[as.character(kept)]
    tot <- sum(supp)
    if (tot <= 0) next
    
    c_star <- kept[which.max(supp)]
    purity <- as.numeric(max(supp) / tot)
    if (purity < tau) next
    
    # Distance guard: Check if within 75th percentile of target intra-cluster distances
    js_c <- j_ok[lab_nb[ok] == c_star]
    med_d <- stats::median(M[i, js_c], na.rm = TRUE)
    q75 <- {
      idx <- which(memb == c_star)
      if (length(idx) < 2) Inf else stats::quantile(as.dist(M[idx, idx]), 0.75, na.rm = TRUE)
    }
    if (!is.finite(med_d) || med_d > q75) next
    
    xi <- X[i, ]
    xN <- colSums(X[js_c, , drop = FALSE]) > 0
    if (sum(xi & xN) < shared_min) next
    
    # Permutation test
    s_obs <- max(supp)
    s_null <- numeric(B)
    for (b in seq_len(B)) {
      sh <- sample(kept_labels, length(kept_labels), replace = FALSE)
      memb_sh <- memb; memb_sh[kept_idx] <- sh
      lab_nb_sh <- memb_sh[nb$j]
      okb <- which(lab_nb_sh %in% kept)
      
      if (!length(okb)) { s_null[b] <- 0; next }
      
      j_okb <- nb$j[okb]
      w_okb <- nb$w[okb] * w_mult[j_okb]
      s_b <- tapply(w_okb, lab_nb_sh[okb], sum)
      s_null[b] <- if (length(s_b)) max(s_b) else 0
    }
    pvals[ii] <- (1 + sum(s_null >= s_obs, na.rm = TRUE)) / (1 + sum(is.finite(s_null)))
    assign_to[ii] <- c_star
  }
  
  tested <- which(is.finite(pvals) & assign_to > 0)
  if (!length(tested)) return(list(membership = memb, absorbed = integer(0), q = numeric(0)))
  
  q <- rep(NA_real_, length(cand))
  q[tested] <- p.adjust(pvals[tested], "BH")
  accept <- tested[which(q[tested] <= alpha)]
  
  m2 <- memb
  if (length(accept)) m2[as.integer(names(assign_to)[accept])] <- assign_to[accept]
  
  list(membership = m2, absorbed = as.integer(names(assign_to)[accept]), q = stats::setNames(q[accept], names(assign_to)[accept]))
}

# ------------------------------------------------------------------------------
# Null Models (Column Shuffle)
# ------------------------------------------------------------------------------

null_column_shuffle <- function(DX, KNN_K, KNN_VARIANT, MIN_CLUSTER_SIZE, reps = 300L, seed = 42) {
  # Preserve marginal prevalence of diagnoses but destroy correlation structure.
  cols <- setdiff(names(DX), "participant_id")
  do_one <- function(r) {
    set.seed(seed + r)
    DXs <- DX
    for (v in cols) DXs[[v]] <- sample(DXs[[v]])
    
    dd <- dedup_dx(DXs)
    DXu <- dd$DXu
    keep <- rowSums(DXu[, setdiff(names(DXu), c("participant_id", "mult")), drop = FALSE], na.rm = TRUE) > 0L
    DXu_id <- DXu[keep, , drop = FALSE]
    
    if (nrow(DXu_id) < 5) return(list(Q = NA, S = NA, ID = NA))
    
    D <- gower_dist_dx(DXu_id)
    if (is.null(D)) return(list(Q = NA, S = NA, ID = NA))
    
    g <- try(knn_graph_from_dist(D, k = KNN_K, variant = KNN_VARIANT, local_scale = TRUE, mult = DXu_id$mult), silent = TRUE)
    if (inherits(g, "try-error") || igraph::ecount(g) == 0L) return(list(Q = NA, S = NA, ID = NA))
    
    memb <- community_labels(g)
    Q <- igraph::modularity(g, memb, weights = igraph::E(g)$weight)
    
    S <- try({
      ss <- cluster::silhouette(memb, D)
      mean(ss[, "sil_width"])
    }, silent = TRUE)
    S <- if (inherits(S, "try-error")) NA_real_ else S
    
    ID <- suppressWarnings(twonn_id_from_dist(D))
    list(Q = Q, S = S, ID = ID)
  }
  
  progressr::handlers(global = TRUE)
  progressr::with_progress({
    p <- progressr::progressor(steps = reps)
    out <- future.apply::future_lapply(seq_len(reps), function(r) {
      z <- do_one(r)
      p()
      z
    }, future.seed = TRUE)
    as.data.frame(do.call(rbind, lapply(out, as.data.frame)))
  })
}

# ------------------------------------------------------------------------------
# Utilities
# ------------------------------------------------------------------------------

auc_weighted <- function(y, s, w = NULL) {
  # Calculate AUC with support for instance weights and tie handling
  y <- as.integer(y); s <- as.numeric(s)
  if (is.null(w)) w <- rep(1, length(y)) else w <- as.numeric(w)
  
  ok <- is.finite(y) & is.finite(s) & is.finite(w) & w > 0
  y <- y[ok]; s <- s[ok]; w <- w[ok]
  
  if (!length(y)) return(NA_real_)
  wp_total <- sum(w[y == 1L])
  wn_total <- sum(w[y == 0L])
  if (wp_total <= 0 || wn_total <= 0) return(NA_real_)
  
  o <- order(s, method = "radix")
  y <- y[o]; s <- s[o]; w <- w[o]
  
  r <- rle(s)
  ends <- cumsum(r$lengths)
  starts <- c(1L, head(ends, -1L) + 1L)
  
  numer <- 0
  cum_wn <- 0
  for (i in seq_along(starts)) {
    idx <- starts[i]:ends[i]
    wn <- sum(w[idx][y[idx] == 0L])
    wp <- sum(w[idx][y[idx] == 1L])
    numer <- numer + wp * cum_wn + 0.5 * wp * wn # Handle ties
    cum_wn <- cum_wn + wn
  }
  numer / (wp_total * wn_total)
}

expand_membership <- function(memb_u, sig_u, sig_all, ids_all) {
  # Map deduplicated cluster labels back to all participants
  lu <- setNames(memb_u, sig_u)
  ma <- lu[as.character(sig_all)]
  data.frame(participant_id = ids_all, cluster = as.integer(ma))
}

# ==============================================================================
# Pipeline Execution
# ==============================================================================

# 1. Data Ingestion
DX <- load_dx_wide(DIAG_CSV, include_nodiag = INCLUDE_NODIAG, min_prev = DX_MIN_PREV, max_prev = DX_MAX_PREV)
stopifnot(nrow(DX) >= MIN_N_DX)

dd <- dedup_dx(DX)
DXu <- dd$DXu
colsU <- setdiff(names(DXu), c("participant_id", "mult"))
keepU <- rowSums(DXu[, colsU, drop = FALSE], na.rm = TRUE) > 0L
DXu_id <- DXu[keepU, , drop = FALSE]

# Prepare metadata for reconstruction
sigU <- apply(DXu[, colsU, drop = FALSE], 1, paste0, collapse = "")
sigU_id <- sigU[keepU]
multU <- as.numeric(DXu_id$mult)
multU[!is.finite(multU) | multU <= 0] <- 1

counts_all <- build_counts(DX)

# 2. Distance Matrix & Intrinsic Dimension
D_dx <- gower_dist_dx(DXu_id)
if (is.null(D_dx)) stop("DX after dedup has <2 varying diagnosis columns.")

Dm <- as.matrix(D_dx)
diag(Dm) <- Inf

# Estimate ID (global and core)
ID_twonn_all <- suppressWarnings(twonn_id_from_dist(D_dx))
idx_core <- core_band_idx(D_dx, k = max(5, min(20, KNN_K)), band = c(0.15, 0.85))
ID_twonn_core <- if (length(idx_core) >= 5) suppressWarnings(twonn_id_from_dist(dist(Dm[idx_core, idx_core, drop = FALSE]))) else NA_real_
ID_lbmle_core <- if (length(idx_core) >= 5) suppressWarnings(lb_mle_id(Dm[idx_core, idx_core, drop = FALSE], 5, 15)) else NA_real_

# 3. kNN Graph & Initial Clustering
g <- knn_graph_from_dist(D_dx, k = KNN_K, variant = KNN_VARIANT, local_scale = LOCAL_SCALE, mult = multU)
memb0 <- community_labels(g)

# Filter communities by size and weight guards
tab_size <- table(memb0)
keep_c <- as.integer(names(tab_size)[tab_size >= MIN_CLUSTER_SIZE])
w_by_c <- tapply(multU, memb0, sum)
keep_c <- intersect(keep_c, as.integer(names(w_by_c)[w_by_c >= MIN_CLUSTER_WEIGHT]))

kept_idx <- which(memb0 %in% keep_c)

if (length(kept_idx) >= 2L) {
  g_kept <- igraph::induced_subgraph(g, vids = kept_idx)
  memb_kept <- as.integer(factor(memb0[kept_idx], levels = keep_c))
  Q <- igraph::modularity(g_kept, memb_kept, weights = igraph::E(g_kept)$weight)
  
  # Validation: Null distribution of Modularity (degree-preserving)
  Q_null_deg <- degree_null_Q(g_kept, memb_kept, reps = N_PERM)
  write_csv(data.frame(Q_null = Q_null_deg, Q_obs = Q), "modularity_degree_null.csv")
  
  # Validation: Silhouette width (on kept subgraph)
  D_kept <- stats::as.dist(Dm[kept_idx, kept_idx, drop = FALSE])
  S_obs <- try({
    ss <- cluster::silhouette(memb_kept, D_kept)
    mean(ss[, "sil_width"])
  }, silent = TRUE)
  S_obs <- if (inherits(S_obs, "try-error")) NA_real_ else S_obs
} else {
  Q <- NA_real_; Q_null_deg <- rep(NA_real_, N_PERM); S_obs <- NA_real_
  write_csv(data.frame(Q_null = Q_null_deg, Q_obs = Q), "modularity_degree_null.csv")
}

# 4. Absorb leftovers
memb_keep_full <- integer(igraph::vcount(g))
if (length(kept_idx)) {
  memb_keep_full[kept_idx] <- as.integer(factor(memb0[kept_idx], levels = keep_c))
}
write_csv(data.frame(), "cluster_silhouette.csv") 

abs <- absorb_leftovers(
  g, D_dx, memb_keep_full, DXu_id, multU,
  tau = 0.70, k_min = 3, shared_min = 1, B = 200, alpha = 0.05, seed = SEED
)
memb2 <- abs$membership
if (length(abs$absorbed)) {
  write_csv(data.frame(v_idx = abs$absorbed, q = as.numeric(abs$q)), "absorbed_vertices_q.csv")
}
KEPT_CLUSTERS <- sort(setdiff(unique(memb2), 0L))

# 5. Validation Metrics
# Full silhouette table
sil_all <- try(silhouette_graph(g, memb2), silent = TRUE)
if (!inherits(sil_all, "try-error")) {
  sil_df <- data.frame(i = seq_len(nrow(sil_all)), cluster = sil_all[, 1], sil_width = sil_all[, 3])
  write_csv(sil_df, "cluster_silhouette.csv")
}

# Bootstrap ARI for stability check
ari <- bootstrap_ari(Dm, KNN_K, KNN_VARIANT, LOCAL_SCALE, COMMUNITY_ALGO, mult = multU, B = N_BOOT)
write_csv(data.frame(boot = seq_along(ari), ari = ari), "cluster_bootstrap_ari.csv")

# Null Model: Column-shuffle (prevalence preserving)
old_plan <- future::plan()
on.exit(future::plan(old_plan), add = TRUE)
future::plan(future::multisession, workers = max(1, parallel::detectCores() - 1))

nulls_col <- null_column_shuffle(DX, KNN_K, KNN_VARIANT, MIN_CLUSTER_SIZE, reps = COL_SHUFFLE_REPS, seed = SEED)
Q_p_col <- mean(nulls_col$Q >= Q, na.rm = TRUE)
S_p_col <- mean(nulls_col$S >= S_obs, na.rm = TRUE)
ID_p_high <- mean(nulls_col$ID >= ID_twonn_all, na.rm = TRUE)
ID_p_low <- mean(nulls_col$ID <= ID_twonn_all, na.rm = TRUE)
ID_p_two <- 2 * min(ID_p_high, ID_p_low)
write_csv(nulls_col, "column_shuffle_nulls.csv")

# 6. Characterization (Pillars A/B/C)
clusters_all <- expand_membership(memb2, sigU_id, dd$sig_all, DX$participant_id)
clusters_all <- subset(clusters_all, cluster %in% KEPT_CLUSTERS)
write_csv(clusters_all, "cluster_membership_all_participants.csv")

# Pillar A: Enrichment
enr <- pillar_A(DX, clusters_all,
                alpha_fdr = ALPHA_FDR, min_prev_in = MIN_PREV_IN_CL, min_or = MIN_OR,
                min_in = MIN_CASES_IN, min_out = MIN_CASES_OUT, min_total = MIN_CASES_TOTAL
)
if (!is.null(enr$enrichment)) {
  write_csv(enr$enrichment, "cluster_diagnosis_enrichment.csv")
  write_csv(data.frame(major_diagnosis = enr$majors), "selected_major_diagnoses.csv")
} else {
  write_csv(data.frame(), "cluster_diagnosis_enrichment.csv")
  write_csv(data.frame(major_diagnosis = character(0)), "selected_major_diagnoses.csv")
}

# Pillar B: Localization
loc_tab <- pillar_B(g, DXu_id, B = 1000, n_pos_min = MIN_CASES_TOTAL, n_neg_min = MIN_CASES_TOTAL)
write_csv(loc_tab, "dx_label_localization.csv")
maj_B <- with(loc_tab, dx[is.finite(assort_p) & is.finite(purity_p) & assort_p <= ALPHA_LOCALIZE & purity_p <= ALPHA_LOCALIZE])

# Pillar C: Predictability
auc_tab <- pillar_C(DXu_id, k = 10, pos_min = MIN_CASES_TOTAL, neg_min = MIN_CASES_TOTAL)
write_csv(auc_tab, "dx_predictability_auc_knn.csv")
maj_C <- subset(auc_tab, is.finite(auc) & auc >= AUC_MIN)$dx

# Determine major diagnoses via union of criteria
eligible <- with(counts_all, dx[prev >= PREV_MIN | n1 >= NCASE_MIN])
inc_nodiag_majors <- if (INCLUDE_NODIAG) c("ANY_DX") else c("ANY_DX", "NODIAG")
maj_union <- sort(setdiff(intersect(unique(c(enr$majors, maj_B, maj_C)), eligible), inc_nodiag_majors))
if (INCLUDE_NODIAG && "NODIAG" %in% names(DX)) {
  maj_union <- union(maj_union, "NODIAG")
}
if (DX_DENY_NOS) maj_union <- maj_union[!grepl("NOS$", maj_union, ignore.case = TRUE)]

# ------------------------------------------------------------------------------
# Sensitivity Analysis (Hyperparameter Sweep)
# ------------------------------------------------------------------------------

build_dist <- function(M) {
  if (TRUE) { # Fallback to local Gower if loaded
    tmp <- data.frame(participant_id = rownames(M), M, check.names = FALSE)
    gower_dist_dx(tmp)
  } else if (requireNamespace("vegan", quietly = TRUE)) {
    vegan::vegdist(M, method = "jaccard", binary = TRUE)
  } else {
    stats::dist(M, method = "binary")
  }
}

get_membership <- function(D, k, variant, mult, B = SWEEP_BOOT_B, deg_reps = DEG_REPS,
                           null_scope = NULL_SCOPE, seed = SEED_BOOT, local_scale = LOCAL_SCALE) {
  # Generate clustering solution for a specific parameter set
  null_scope <- match.arg(null_scope, c("kept", "full"))
  stopifnot(inherits(D, "dist") || is.matrix(D))
  M <- if (inherits(D, "dist")) as.matrix(D) else D
  vnames <- rownames(M) %||% as.character(seq_len(nrow(M)))
  
  g <- knn_graph_from_dist(D, k = k, variant = variant, local_scale = local_scale, mult = mult)
  lv <- withr::with_seed(as.integer(seed), community_labels(g, algo = COMMUNITY_ALGO, obj = LEIDEN_OBJECTIVE, gamma = LEIDEN_GAMMA, iters = LEIDEN_ITERS))
  
  tab <- table(lv)
  keep <- as.integer(names(tab)[tab >= MIN_CLUSTER_SIZE])
  w_by <- tapply(mult, lv, sum)
  keep <- intersect(keep, as.integer(names(w_by)[w_by >= MIN_CLUSTER_WEIGHT]))
  
  kept_idx <- which(lv %in% keep)
  if (length(kept_idx) < 2L) {
    return(list(m = setNames(integer(igraph::vcount(g)), vnames), kept_idx = integer(0), kept_names = character(0), Q = NA_real_, Q_null_w = rep(NA_real_, B), Q_null_deg = rep(NA_real_, deg_reps), z_w = NA_real_, z_deg = NA_real_, n_kept = 0L))
  }
  
  gk <- igraph::induced_subgraph(g, vids = kept_idx)
  lv_k <- as.integer(factor(lv[kept_idx], levels = keep))
  Q <- igraph::modularity(gk, lv_k, weights = igraph::E(gk)$weight)
  
  m <- setNames(integer(igraph::vcount(g)), vnames)
  m[kept_idx] <- lv_k
  
  # Null models for specific sweep run
  Q_null_w <- rep(NA_real_, B)
  Q_null_deg <- rep(NA_real_, deg_reps)
  
  if (null_scope == "kept") {
    Wk <- igraph::E(gk)$weight
    for (b in seq_len(B)) {
      withr::with_seed(as.integer(seed) + b, {
        idx <- sample.int(length(Wk), replace = FALSE)
        igraph::E(gk)$weight <- Wk[idx]
        mb_k <- community_labels(gk, algo = COMMUNITY_ALGO, obj = LEIDEN_OBJECTIVE, gamma = LEIDEN_GAMMA, iters = LEIDEN_ITERS)
        mb_k <- as.integer(factor(mb_k))
        Q_null_w[b] <- igraph::modularity(gk, mb_k, weights = igraph::E(gk)$weight)
      })
    }
    igraph::E(gk)$weight <- Wk
    Q_null_deg <- suppressWarnings(withr::with_seed(as.integer(seed) + 10^6L, degree_null_Q(gk, lv_k, reps = deg_reps)))
  } else { 
    # Scope: Full graph
    W <- igraph::E(g)$weight
    for (b in seq_len(B)) {
      withr::with_seed(as.integer(seed) + b, {
        idx <- sample.int(length(W), replace = FALSE)
        igraph::E(g)$weight <- W[idx]
        mb <- community_labels(g, algo = COMMUNITY_ALGO, obj = LEIDEN_OBJECTIVE, gamma = LEIDEN_GAMMA, iters = LEIDEN_ITERS)
        tabb <- table(mb)
        keepb <- as.integer(names(tabb)[tabb >= MIN_CLUSTER_SIZE])
        wb <- tapply(mult, mb, sum)
        keepb <- intersect(keepb, as.integer(names(wb)[wb >= MIN_CLUSTER_WEIGHT]))
        kept_b <- which(mb %in% keepb)
        if (length(kept_b) >= 2L) {
          gb <- igraph::induced_subgraph(g, vids = kept_b)
          mb_k <- as.integer(factor(mb[kept_b], levels = keepb))
          Q_null_w[b] <- igraph::modularity(gb, mb_k, weights = igraph::E(gb)$weight)
        } else {
          Q_null_w[b] <- NA_real_
        }
      })
    }
    igraph::E(g)$weight <- W
    if (NULL_SCOPE == "kept") {
      Q_null_deg <- suppressWarnings(withr::with_seed(as.integer(seed) + 10^6L, degree_null_Q(gk, lv_k, reps = deg_reps)))
    }
  }
  
  mu_w <- mean(Q_null_w, na.rm = TRUE)
  sd_w <- stats::sd(Q_null_w, na.rm = TRUE)
  z_w <- if (is.finite(mu_w) && is.finite(sd_w) && sd_w > 0) (Q - mu_w) / sd_w else NA_real_
  
  mu_d <- mean(Q_null_deg, na.rm = TRUE)
  sd_d <- stats::sd(Q_null_deg, na.rm = TRUE)
  z_d <- if (is.finite(mu_d) && is.finite(sd_d) && sd_d > 0) (Q - mu_d) / sd_d else NA_real_
  
  list(m = m, kept_idx = kept_idx, kept_names = vnames[kept_idx], Q = Q, Q_null_w = Q_null_w, Q_null_deg = Q_null_deg, z_w = z_w, z_deg = z_d, n_kept = length(unique(lv_k)))
}

align_ari <- function(a, b) {
  # Compute ARI between two clusterings (filtering unassigned nodes)
  if (is.null(names(a))) names(a) <- as.character(seq_along(a))
  if (is.null(names(b))) names(b) <- as.character(seq_along(b))
  ii <- intersect(names(a), names(b))
  aa <- as.integer(a[ii])
  bb <- as.integer(b[ii])
  keep <- (aa != 0L) & (bb != 0L)
  aa <- aa[keep]; bb <- bb[keep]
  if (length(aa) < 2 || length(unique(aa)) < 2 || length(unique(bb)) < 2) return(NA_real_)
  aricode::ARI(aa, bb)
}

compute_silhouette_with_colshuffle <- function(DX, ids_all, kept_idx, labs_kept, reps = SIL_REPS, seed = SEED + 1e5L) {
  # Evaluate silhouette significance via column-shuffled nulls
  if (length(kept_idx) < 3L || length(unique(labs_kept)) < 2L) return(list(S_obs = NA_real_, S_p = NA_real_))
  
  dx_cols <- setdiff(names(DX), c("participant_id", "ANY_DX", "mult"))
  Xall <- as.matrix(DX[match(ids_all, DX$participant_id), dx_cols, drop = FALSE])
  if (anyNA(Xall)) Xall[is.na(Xall)] <- 0L 
  storage.mode(Xall) <- "numeric"
  rownames(Xall) <- ids_all
  
  Dk <- as.dist(as.matrix(D_dx)[kept_idx, kept_idx, drop = FALSE])
  Sil <- try(cluster::silhouette(labs_kept, Dk), silent = TRUE)
  S_obs <- if (!inherits(Sil, "try-error")) mean(Sil[, "sil_width"]) else NA_real_
  if (!is.finite(S_obs)) return(list(S_obs = NA_real_, S_p = NA_real_))
  
  S_null <- rep(NA_real_, reps)
  set.seed(as.integer(seed))
  for (b in seq_len(reps)) {
    Xb <- Xall
    for (j in seq_len(ncol(Xb))) Xb[, j] <- Xb[sample.int(nrow(Xb)), j]
    Db <- build_dist(Xb)
    
    nz_kept <- rowSums(Xb[kept_idx, , drop = FALSE]) > 0
    idx_k <- kept_idx[nz_kept]
    labs_k <- labs_kept[nz_kept]
    
    if (length(idx_k) >= 2L && length(unique(labs_k)) >= 2L) {
      Db_k <- as.dist(as.matrix(Db)[idx_k, idx_k, drop = FALSE])
      Sb <- try(cluster::silhouette(labs_k, Db_k), silent = TRUE)
      S_null[b] <- if (!inherits(Sb, "try-error")) mean(Sb[, "sil_width"]) else NA_real_
    } else {
      S_null[b] <- NA_real_
    }
  }
  
  B <- sum(is.finite(S_null))
  S_p <- if (B > 0) (sum(S_null >= S_obs, na.rm = TRUE) + 1) / (B + 1) else NA_real_
  list(S_obs = S_obs, S_p = S_p)
}

if (isTRUE(DO_SWEEP)) {
  # Baseline for ARI comparison
  base <- get_membership(
    D = D_dx, k = KNN_K, variant = KNN_VARIANT,
    mult = multU, B = N_PERM, deg_reps = DEG_REPS,
    null_scope = "kept", seed = SEED, local_scale = LOCAL_SCALE
  )
  m_base <- base$m
  
  grid_legacy <- expand.grid(K = K_GRID, variant = VARIANT_GRID, stringsAsFactors = FALSE)
  seed_vec <- SEED + seq_len(nrow(grid_legacy))
  
  old_plan <- future::plan()
  on.exit(future::plan(old_plan), add = TRUE)
  future::plan(future::multisession, workers = NCORES_PAR)
  
  sens <- future.apply::future_lapply(seq_len(nrow(grid_legacy)), function(i) {
    k <- grid_legacy$K[i]
    v <- grid_legacy$variant[i]
    gi <- get_membership(D = D_dx, k = k, variant = v, mult = multU, B = N_PERM, deg_reps = DEG_REPS, null_scope = "kept", seed = seed_vec[i], local_scale = LOCAL_SCALE)
    
    p_pref <- {
      Bq <- sum(is.finite(gi$Q_null_deg))
      if (Bq > 0) (sum(gi$Q_null_deg >= gi$Q, na.rm = TRUE) + 1) / (Bq + 1) else NA_real_
    }
    
    Q_p_two <- {
      Bw <- sum(is.finite(gi$Q_null_w))
      if (Bw > 0) {
        mu <- mean(gi$Q_null_w, na.rm = TRUE)
        sd <- stats::sd(gi$Q_null_w, na.rm = TRUE)
        if (is.finite(mu) && is.finite(sd) && sd > 0) 2 * stats::pnorm(-abs((gi$Q - mu) / sd)) else NA_real_
      } else { NA_real_ }
    }
    
    S <- NA_real_
    lab_all <- as.integer(gi$m)
    kk <- gi$kept_idx
    if (length(kk) > 1L) {
      lab_k <- as.integer(factor(lab_all[kk]))
      tabk <- table(lab_k)
      keep_lab <- as.integer(names(tabk)[tabk >= 2L]) 
      sel_k <- lab_k %in% keep_lab
      if (sum(sel_k) >= 2L && length(keep_lab) >= 2L) {
        Dk <- as.matrix(D_dx)[kk, kk, drop = FALSE]
        Sil <- try(cluster::silhouette(lab_k[sel_k], as.dist(Dk[sel_k, sel_k, drop = FALSE])), silent = TRUE)
        if (!inherits(Sil, "try-error")) S <- mean(Sil[, "sil_width"], na.rm = TRUE)
      }
    }
    
    S_p <- NA_real_
    if (is.finite(S)) {
      ids_all <- as.character(DXu_id$participant_id)
      cs <- compute_silhouette_with_colshuffle(DX, ids_all, kept_idx = kk[sel_k], labs_kept = lab_k[sel_k], reps = SIL_REPS, seed = seed_vec[i] + 1e5L)
      S <- cs$S_obs
      S_p <- cs$S_p
    }
    
    data.frame(
      K = k, variant = v, Q = gi$Q, p_pref = p_pref, Q_p_two = Q_p_two,
      S = S, S_p = S_p, z_w = gi$z_w, z_deg = gi$z_deg, n_kept = gi$n_kept,
      ARI_vs_base = align_ari(m_base, gi$m), stringsAsFactors = FALSE
    )
  }, future.seed = TRUE) |> dplyr::bind_rows()
  
  write_csv(sens, "clustering_sensitivity_grid.csv")
  
  # Summary of stability across sweep
  grid_new <- expand.grid(k = K_GRID, variant = VARIANT_GRID, stringsAsFactors = FALSE)
  rows <- vector("list", nrow(grid_new))
  for (i in seq_len(nrow(grid_new))) {
    k <- grid_new$k[i]
    v <- grid_new$variant[i]
    Q_rep <- numeric(SWEEP_BOOT_B)
    K_rep <- integer(SWEEP_BOOT_B)
    for (b in seq_len(SWEEP_BOOT_B)) {
      withr::with_seed(SEED + 1000L * i + b, {
        g_i <- knn_graph_from_dist(D_dx, k = k, variant = v, local_scale = LOCAL_SCALE, mult = multU)
        lab_i <- community_labels(g_i)
      })
      tab_i <- table(lab_i)
      keep_i <- as.integer(names(tab_i)[tab_i >= MIN_CLUSTER_SIZE])
      w_i <- tapply(multU, lab_i, sum)
      keep_i <- intersect(keep_i, as.integer(names(w_i)[w_i >= MIN_CLUSTER_WEIGHT]))
      if (length(keep_i) >= 1L) {
        kept_idx <- which(lab_i %in% keep_i)
        gk <- igraph::induced_subgraph(g_i, vids = kept_idx)
        labs_k <- as.integer(factor(lab_i[kept_idx], levels = keep_i))
        Q_rep[b] <- igraph::modularity(gk, labs_k, weights = igraph::E(gk)$weight)
        K_rep[b] <- length(unique(labs_k))
      } else {
        Q_rep[b] <- NA_real_; K_rep[b] <- 0L
      }
    }
    rows[[i]] <- data.frame(
      k = k, variant = v, Q = mean(Q_rep, na.rm = TRUE), Q_sd = stats::sd(Q_rep, na.rm = TRUE),
      n_clusters = as.integer(stats::median(K_rep, na.rm = TRUE)), n_clusters_iqr = stats::IQR(K_rep, na.rm = TRUE),
      reps = SWEEP_BOOT_B
    )
  }
  sweep_df <- dplyr::bind_rows(rows)
  write_csv(sweep_df, "clustering_sensitivity_grid_summary.csv")
}

# ------------------------------------------------------------------------------
# Delta Modularity Analysis
# ------------------------------------------------------------------------------

if (DO_DELTAQ) {
  Q_base <- Q
  dqs <- list()
  j_idx <- 0L
  for (j in setdiff(names(DXu_id), c("participant_id", "mult", "ANY_DX"))) {
    j_idx <- j_idx + 1L
    DXm <- DX; DXm[[j]] <- 0L
    
    dd0 <- dedup_dx(DXm); DXu0 <- dd0$DXu
    keepr <- rowSums(DXu0[, setdiff(names(DXu0), c("participant_id", "mult")), drop = FALSE], na.rm = TRUE) > 0L
    DXu0 <- DXu0[keepr, , drop = FALSE]
    if (nrow(DXu0) < 5) next
    
    D0 <- gower_dist_dx(DXu0)
    if (is.null(D0)) next
    
    g0 <- knn_graph_from_dist(D0, k = KNN_K, variant = KNN_VARIANT, local_scale = LOCAL_SCALE, mult = DXu0$mult)
    lab0 <- withr::with_seed(SEED + j_idx, community_labels(g0))
    
    tab0 <- table(lab0)
    keep0 <- as.integer(names(tab0)[tab0 >= MIN_CLUSTER_SIZE])
    w0 <- tapply(DXu0$mult, lab0, sum)
    keep0 <- intersect(keep0, as.integer(names(w0)[w0 >= MIN_CLUSTER_WEIGHT]))
    
    if (length(keep0) >= 1L) {
      kept_idx0 <- which(lab0 %in% keep0)
      g0k <- igraph::induced_subgraph(g0, vids = kept_idx0)
      lab0k <- as.integer(factor(lab0[kept_idx0], levels = keep0))
      Q0 <- igraph::modularity(g0k, lab0k, weights = igraph::E(g0k)$weight)
    } else {
      Q0 <- NA_real_
    }
    dqs[[length(dqs) + 1L]] <- data.frame(dx = j, dQ = Q0 - Q_base)
  }
  dq_df <- dplyr::bind_rows(dqs)
  write_csv(dq_df, "dx_deltaQ.csv")
}

# ------------------------------------------------------------------------------
# Over-representation Analysis (Figures)
# ------------------------------------------------------------------------------

compute_overrep <- function(DX, clusters_all,
                            min_in = NIN_MIN %||% 5L,
                            min_out = NOUT_MIN %||% 5L,
                            dx_filter = NULL,
                            bh_scope = c("global", "per_cluster"),
                            min_pos_in = 0L, min_pos_out = 0L, min_total_pos = 0L) {
  # Compute enrichment lift and Fisher p-values with BH correction
  bh_scope <- match.arg(bh_scope)
  stopifnot(all(c("participant_id", "cluster") %in% names(clusters_all)))
  
  DX <- as.data.frame(DX)
  ids_ref <- intersect(as.character(DX$participant_id), as.character(clusters_all$participant_id))
  if (!length(ids_ref)) return(data.frame())
  
  clmap <- setNames(clusters_all$cluster, clusters_all$participant_id)
  clv <- as.integer(clmap[ids_ref])
  clv[is.na(clv)] <- 0L
  
  cols <- setdiff(names(DX), c("participant_id", "ANY_DX"))
  if (!is.null(dx_filter)) cols <- intersect(cols, dx_filter)
  if (!length(cols)) return(data.frame())
  
  DXb <- DX[match(ids_ref, DX$participant_id), cols, drop = FALSE]
  
  out <- list()
  for (cid in sort(setdiff(unique(clv), 0L))) {
    inidx <- which(clv == cid); outidx <- which(clv != cid)
    if (length(inidx) < min_in || length(outidx) < min_out) next
    n_in <- length(inidx); n_out <- length(outidx)
    
    for (dx in cols) {
      v <- as.integer(DXb[[dx]] == 1L); v[is.na(v)] <- 0L
      a <- sum(v[inidx] == 1L); b <- n_in - a
      c <- sum(v[outidx] == 1L); d <- n_out - c
      if ((a + c) == 0 || (b + d) == 0) next
      
      in_prev <- a / n_in
      out_prev <- if ((c + d) > 0) c / (c + d) else NA_real_
      lift <- if (is.finite(out_prev) && out_prev > 0) in_prev / out_prev else Inf
      
      eligible <- ((a + c) >= min_total_pos) && (a >= min_pos_in) && (c >= min_pos_out)
      dir_sign <- sign(a * n_out - c * n_in)
      alt <- if (dir_sign < 0) "less" else if (dir_sign > 0) "greater" else "two.sided"
      
      p <- if (eligible) {
        tryCatch(suppressWarnings(fisher.test(matrix(c(a, b, c, d), 2, byrow = TRUE), alternative = alt)$p.value), error = function(e) NA_real_)
      } else { NA_real_ }
      
      out[[length(out) + 1L]] <- data.frame(cluster = cid, dx = dx, n_in = n_in, n_out = n_out, a = a, b = b, c = c, d = d, in_prev = in_prev, out_prev = out_prev, lift = lift, p = p, stringsAsFactors = FALSE)
    }
  }
  T <- dplyr::bind_rows(out)
  if (!nrow(T)) return(T)
  
  if (bh_scope == "per_cluster") {
    T$q <- NA_real_
    for (cid in sort(unique(T$cluster))) {
      i <- which(T$cluster == cid & is.finite(T$p))
      if (length(i)) T$q[i] <- p.adjust(T$p[i], "BH")
    }
  } else {
    T$q <- ifelse(is.finite(T$p), p.adjust(T$p, "BH"), NA_real_)
  }
  
  T$star <- dplyr::case_when(
    is.finite(T$q) & T$q < 0.01 ~ "**",
    is.finite(T$q) & T$q < 0.05 ~ "*",
    !is.finite(T$q) ~ "",
    TRUE ~ {if (isTRUE(SHOW_NS_LABELS)) "ns" else ""}
  )
  
  T$log2_lift <- log2(pmax(T$lift, 1e-12))
  T[order(T$cluster, -T$log2_lift, T$q), ]
}

tab_sig <- compute_overrep(DX, clusters_all, min_in = NIN_MIN, min_out = NOUT_MIN)
write_csv(tab_sig, "dx_overrep_full.csv")

# Export for heatmap
if (nrow(tab_sig) > 0) {
  arch_wide <- tab_sig %>%
    dplyr::select(dx, cluster, lift) %>%
    dplyr::mutate(cluster = paste0("C", cluster)) %>%
    tidyr::pivot_wider(names_from = "cluster", values_from = "lift", values_fill = 1.0)
  
  expected_cols <- paste0("C", sort(unique(clusters_all$cluster)))
  missing_cols <- setdiff(expected_cols, names(arch_wide))
  if (length(missing_cols)) {
    for (mc in missing_cols) arch_wide[[mc]] <- 1.0
  }
  write_csv(arch_wide, "cluster_archetypes_dx_Kbest.csv")
} else {
  write_csv(data.frame(dx = character(0)), "cluster_archetypes_dx_Kbest.csv")
}

# NODIAG Sensitivity check
c_nodiag <- NA_integer_
if ("NODIAG" %in% names(DX) && nrow(tab_sig)) {
  nd <- subset(tab_sig, dx == "NODIAG" & (lift > 1 | !is.finite(lift)))
  if (nrow(nd)) c_nodiag <- nd$cluster[which.max(nd$lift)]
}
if (is.finite(c_nodiag)) {
  clusters_no_nodiag <- subset(clusters_all, cluster != c_nodiag)
  tab_sig2 <- compute_overrep(DX, clusters_no_nodiag, min_in = NIN_MIN, min_out = NOUT_MIN)
  write_csv(tab_sig2, "dx_overrep_noNODIAG.csv")
}

# ------------------------------------------------------------------------------
# Seed-Bootstrap Sweep (Full Distribution)
# ------------------------------------------------------------------------------

baseline <- get_membership(
  D_dx, k = KNN_K, variant = KNN_VARIANT, mult = multU,
  B = SWEEP_BOOT_B, deg_reps = DEG_REPS, null_scope = NULL_SCOPE, seed = SEED, local_scale = LOCAL_SCALE
)
base_m <- baseline$m

grid <- expand.grid(K = K_GRID, variant = VARIANT_GRID, stringsAsFactors = FALSE)
future::plan(future::multisession, workers = NCORES_PAR)

jobs <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  data.frame(i = i, K = grid$K[i], variant = grid$variant[i], seed = SEED_BOOT + seq_len(SWEEP_BOOT_B), stringsAsFactors = FALSE)
}))

progressr::handlers(global = TRUE)
boot_raw <- progressr::with_progress({
  p <- progressr::progressor(steps = nrow(jobs))
  res <- future_lapply(seq_len(nrow(jobs)), function(j) {
    k <- jobs$K[j]; v <- jobs$variant[j]; s <- jobs$seed[j]
    gi <- get_membership(D_dx, k = k, variant = v, mult = multU, B = N_PERM, deg_reps = DEG_REPS, null_scope = NULL_SCOPE, seed = s, local_scale = LOCAL_SCALE)
    p()
    data.frame(K = k, variant = v, Q = gi$Q, z_w = gi$z_w, n_kept = gi$n_kept, ARI_vs_base = align_ari(base_m, gi$m), stringsAsFactors = FALSE)
  }, future.seed = TRUE)
  dplyr::bind_rows(res)
})

write_csv(boot_raw, "clustering_sensitivity_seed_bootstrap_long.csv")

summ_boot <- boot_raw %>%
  dplyr::group_by(K, variant) %>%
  dplyr::summarise(n_runs = dplyr::n(), Q_med = stats::median(Q, na.rm = TRUE), Q_sd = stats::sd(Q, na.rm = TRUE), z_w_med = stats::median(z_w, na.rm = TRUE), ARI_med = stats::median(ARI_vs_base, na.rm = TRUE), ARI_IQR = IQR(ARI_vs_base, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(Q_med), dplyr::desc(ARI_med))

write_csv(summ_boot, "clustering_sensitivity_seed_bootstrap_summary.csv")

# ------------------------------------------------------------------------------
# Major Diagnosis Stability
# ------------------------------------------------------------------------------

dx_pool <- setdiff(names(DX), c("participant_id", "ANY_DX"))
zero_vec <- function() setNames(integer(length(dx_pool)), dx_pool)

boot_one <- function(b) {
  success <- zero_vec(); trials <- zero_vec(); valid <- 0L
  
  withr::with_seed(SEED_BOOT + b, { idx <- sample.int(nrow(DX), replace = TRUE) })
  DX_b <- DX[idx, , drop = FALSE]
  dd_b <- dedup_dx(DX_b)
  DXu_b <- dd_b$DXu
  colsB <- setdiff(names(DXu_b), c("participant_id", "mult", "ANY_DX"))
  keepB <- rowSums(DXu_b[, colsB, drop = FALSE], na.rm = TRUE) > 0L
  DXu_id_b <- DXu_b[keepB, , drop = FALSE]
  
  if (nrow(DXu_id_b) < 5) return(list(success = success, trials = trials, valid = valid))
  
  D_b <- gower_dist_dx(DXu_id_b)
  if (is.null(D_b)) return(list(success = success, trials = trials, valid = valid))
  
  g_b <- try(knn_graph_from_dist(D_b, k = KNN_K, variant = KNN_VARIANT, local_scale = TRUE, mult = DXu_id_b$mult), silent = TRUE)
  if (inherits(g_b, "try-error") || igraph::ecount(g_b) == 0L) return(list(success = success, trials = trials, valid = valid))
  
  mb <- withr::with_seed(SEED_BOOT + 10^5L + b, community_labels(g_b, algo = COMMUNITY_ALGO, obj = LEIDEN_OBJECTIVE, gamma = LEIDEN_GAMMA, iters = LEIDEN_ITERS))
  
  tab <- table(mb)
  keep <- as.integer(names(tab)[tab >= MIN_CLUSTER_SIZE])
  w <- tapply(DXu_id_b$mult, mb, sum)
  keep <- intersect(keep, as.integer(names(w)[w >= MIN_CLUSTER_WEIGHT]))
  mb[!mb %in% keep] <- 0L
  valid <- 1L
  
  sigU_b <- apply(DXu_id_b[, colsB, drop = FALSE], 1, paste0, collapse = "")
  cl_all_b <- expand_membership(mb, sigU_b, dd_b$sig_all, rownames(DX_b))
  
  enr_b <- pillar_A(DX_b, cl_all_b, alpha_fdr = ALPHA_FDR, min_prev_in = MIN_PREV_IN_CL, min_or = MIN_OR, min_in = NIN_MIN, min_out = NOUT_MIN, min_total = NCASE_MIN)
  mA <- if (!is.null(enr_b$enrichment)) enr_b$majors else character(0)
  
  loc_b <- pillar_B(g_b, DXu_id_b, B = PILLAR_B_REPS, n_pos_min = NCASE_MIN, n_neg_min = NCASE_MIN)
  mB <- subset(loc_b, pmin(assort_p, purity_p) <= ALPHA_LOCALIZE)$dx
  
  auc_b <- pillar_C(DXu_id_b, k = 10, pos_min = NCASE_MIN, neg_min = NCASE_MIN)
  mC <- subset(auc_b, is.finite(auc) & auc >= AUC_MIN)$dx
  
  counts_b <- build_counts(DX_b)
  eligible_b <- with(counts_b, dx[prev >= PREV_MIN | n1 >= NCASE_MIN])
  
  elig_here <- intersect(eligible_b, dx_pool)
  if (length(elig_here)) trials[elig_here] <- trials[elig_here] + 1L
  
  majors_b <- sort(intersect(unique(c(mA, mB, mC)), eligible_b))
  sel_here <- intersect(majors_b, dx_pool)
  if (length(sel_here)) success[sel_here] <- success[sel_here] + 1L
  
  list(success = success, trials = trials, valid = valid)
}

future::plan(future::multisession, workers = NCORES_PAR)
progressr::handlers(global = TRUE)
res <- progressr::with_progress({
  p <- progressr::progressor(steps = MAJORS_BOOT_B)
  parts <- future.apply::future_lapply(seq_len(MAJORS_BOOT_B), function(b) {
    x <- boot_one(b); p(); x
  }, future.seed = TRUE)
  parts
})
future::plan(future::sequential)

success_by_dx <- Reduce(`+`, lapply(res, `[[`, "success"), init = zero_vec())
trials_by_dx <- Reduce(`+`, lapply(res, `[[`, "trials"), init = zero_vec())
n_valid <- sum(vapply(res, `[[`, integer(1), "valid"))

stab <- data.frame(dx = names(success_by_dx), count = as.integer(success_by_dx), trials = as.integer(trials_by_dx), stringsAsFactors = FALSE)
stab <- subset(stab, trials > 0)
stab$freq <- with(stab, count / trials)
stab <- stab[order(stab$freq, decreasing = TRUE), ]
stab$lo <- stats::qbinom(0.025, size = stab$trials, prob = stab$freq) / stab$trials
stab$hi <- stats::qbinom(0.975, size = stab$trials, prob = stab$freq) / stab$trials

write_csv(stab, "major_dx_bootstrap_frequency_by_eligibility.csv")

# ------------------------------------------------------------------------------
# Final Output Summary
# ------------------------------------------------------------------------------

run_summary <- data.frame(
  n_total = nrow(DX), n_unique = nrow(DXu), n_used_nonzero = nrow(DXu_id),
  K_KNN = KNN_K, Q_modularity = Q, silhouette = S_obs,
  ARI_median = median(ari, na.rm = TRUE), ARI_IQR = IQR(ari, na.rm = TRUE),
  ID_TwoNN_all = ID_twonn_all, ID_TwoNN_core = ID_twonn_core, ID_LB_core = ID_lbmle_core,
  n_core = if (length(idx_core)) length(idx_core) else NA_integer_,
  n_clusters_kept = length(KEPT_CLUSTERS),
  p_Q_colshuffle = Q_p_col, p_S_colshuffle = S_p_col, p_ID_two = ID_p_two
)
write_csv(run_summary, "dx_space_run_summary.csv")
print(run_summary)

# ------------------------------------------------------------------------------
# Visualization
# ------------------------------------------------------------------------------

if (DO_PLOTS) {
  # kNN Degree Distribution
  deg <- igraph::degree(g)
  p <- ggplot2::ggplot(data.frame(deg = deg), ggplot2::aes(deg)) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = "kNN graph degree", x = "degree", y = "count")
  print(p)
  save_plot_gg("FIG_knn_degree", p, width = 6.0, height = 4.2, dpi = PNG_DPI)
  
  # Cluster Sizes
  cs <- table(memb2[memb2 != 0L])
  p2 <- ggplot2::ggplot(data.frame(cluster = names(cs), size = as.numeric(cs)), ggplot2::aes(cluster, size)) +
    ggplot2::geom_col() + ggplot2::theme_minimal() + ggplot2::labs(title = "Cluster sizes (kept)")
  print(p2)
  save_plot_gg("FIG_cluster_sizes", p2, width = 6.0, height = 4.2, dpi = PNG_DPI)
  
  # Degree-null validation
  BQ <- sum(is.finite(Q_null_deg))
  Q_p_deg <- (sum(Q_null_deg >= Q, na.rm = TRUE) + 1) / (BQ + 1)
  write_csv(data.frame(Q_null_deg = Q_null_deg, Q_obs = Q, Q_p_deg = Q_p_deg), "modularity_degree_null_kept.csv")
  
  # Table Exports
  loc_v1 <- dplyr::rename(loc_tab, assort_r = assort, knn_purity = purity, knn_p = purity_p)
  write_csv(loc_tab, "dx_label_localization.csv")
  write_csv(auc_tab, "dx_predictability_auc_knn.csv")
  write_csv(data.frame(dx = maj_union), "selected_major_diagnoses_union.csv")
  df_ari <- data.frame(boot = seq_along(ari), ARI = ari)
  write_csv(dplyr::rename(df_ari, ari = ARI), "cluster_bootstrap_ari.csv")
  
  # Helper: Wilson Score Interval
  wilson_ci <- function(k, n, z = 1.96) {
    p <- k / n
    denom <- 1 + z^2 / n
    centre <- (p + z^2 / (2 * n)) / denom
    half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
    cbind(lo = pmax(0, centre - half), hi = pmin(1, centre + half))
  }
  
  plot_overrep_heatmap <- function(tab, file = "FIG_dx_overrep_heatmap.png", top_k_per_cluster = 12L, cap = 2, clusters_on_y = TRUE) {
    if (!nrow(tab)) return(invisible(NULL))
    if (!"log2_lift" %in% names(tab)) tab$log2_lift <- log2(pmax(tab$lift, 1e-12))
    if (!"star" %in% names(tab)) {
      tab$star <- dplyr::case_when(
        is.finite(tab$q) & tab$q < 0.01 ~ "**", is.finite(tab$q) & tab$q < 0.05 ~ "*", is.finite(tab$q) & tab$q >= 0.05 & show_ns ~ "ns", TRUE ~ ""
      )
    }
    
    TOP <- tab %>% dplyr::group_by(cluster) %>% dplyr::arrange(dplyr::desc(log2_lift), .by_group = TRUE) %>% dplyr::filter(is.finite(log2_lift)) %>% dplyr::slice_head(n = top_k_per_cluster) %>% dplyr::ungroup()
    keep_dx <- unique(TOP$dx)
    if (!length(keep_dx)) return(invisible(NULL))
    
    H <- tab %>% dplyr::filter(dx %in% keep_dx)
    H$val <- pmin(pmax(H$log2_lift, -cap), cap)
    cl_lev <- sort(unique(H$cluster))
    H$cluster_f <- factor(paste0("C", H$cluster), levels = paste0("C", cl_lev))
    H$dx_f <- factor(H$dx, levels = sort(unique(keep_dx)))
    
    if (clusters_on_y) {
      H$x_i <- as.numeric(H$dx_f); H$y_i <- as.numeric(H$cluster_f)
      x_breaks <- seq_along(levels(H$dx_f)); x_labels <- levels(H$dx_f)
      y_breaks <- seq_along(levels(H$cluster_f)); y_labels <- levels(H$cluster_f)
      ax_x <- NULL; ax_y <- "cluster"
    } else {
      H$x_i <- as.numeric(H$cluster_f); H$y_i <- as.numeric(H$dx_f)
      x_breaks <- seq_along(levels(H$cluster_f)); x_labels <- levels(H$cluster_f)
      y_breaks <- seq_along(levels(H$dx_f)); y_labels <- levels(H$dx_f)
      ax_x <- "cluster"; ax_y <- NULL
    }
    H <- dplyr::distinct(H, x_i, y_i, .keep_all = TRUE)
    eps <- 0.0017
    H$alpha_lab <- ifelse(H$star == "ns", 0.35, 1)
    
    p <- ggplot2::ggplot(H) +
      ggplot2::geom_tile(ggplot2::aes(x = x_i, y = y_i, fill = val), width = 1 + 2 * eps, height = 1 + 2 * eps, linewidth = 0) +
      ggplot2::geom_text(ggplot2::aes(x = x_i, y = y_i, label = star, alpha = alpha_lab), size = 3, fontface = "bold") +
      ggplot2::scale_alpha_identity() +
      ggplot2::scale_x_continuous(breaks = x_breaks, labels = x_labels, expand = c(0, 0)) +
      ggplot2::scale_y_continuous(breaks = y_breaks, labels = y_labels, expand = c(0, 0)) +
      scale_shared_div(name = "log2(lift)", limits = c(-cap, cap)) +
      ggplot2::labs(x = ax_x, y = ax_y, title = "Over-representation by cluster") +
      ggplot2::theme_minimal(11) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1))
    
    nm <- tools::file_path_sans_ext(basename(file))
    save_plot_gg(nm, p, width = 8.5, height = 5)
    invisible(p)
  }
  
  plot_overrep_lift_table <- function(tab_out, file = "FIG_overrep_lift_table.png", digits = 2) {
    stopifnot(all(c("cluster", "dx", "lift") %in% names(tab_out)))
    df <- tab_out
    if (!"label" %in% names(df)) df$label <- if ("star" %in% names(df)) df$star else ""
    df$cluster <- factor(df$cluster, levels = sort(unique(df$cluster)))
    if (!is.factor(df$dx)) df$dx <- factor(df$dx, levels = unique(df$dx))
    
    cap_lift <- 99.9
    lift_c <- pmin(df$lift, cap_lift)
    lift_txt <- ifelse(is.finite(lift_c), sprintf(paste0("%.", digits, "f"), lift_c), "")
    sig_txt <- ifelse(df$label %in% c("*", "**"), paste0(" ", df$label), "")
    df$cell <- paste0(lift_txt, sig_txt)
    
    p <- ggplot2::ggplot(df, ggplot2::aes(x = df$dx, y = df$cluster)) +
      ggplot2::geom_tile(fill = "white", colour = "grey85") +
      ggplot2::geom_text(ggplot2::aes(label = df$cell), size = 3.2) +
      ggplot2::labs(title = "Lift (in / overall) with significance", x = NULL, y = "cluster", caption = "Cells show lift; * FDR≤0.05, ** FDR≤0.01.") +
      ggplot2::theme_minimal(12) +
      ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, vjust = 1))
    
    nm <- tools::file_path_sans_ext(basename(file))
    save_plot_gg(nm, p, width = 8.5, height = 7)
  }
  
  dx_for_heatmap <- intersect(maj_union, setdiff(names(DX), c("participant_id", "ANY_DX")))
  tab_sig_maj <- subset(tab_sig, dx %in% dx_for_heatmap)
  plot_overrep_heatmap(tab_sig_maj, "FIG_dx_overrep_heatmap.png", top_k_per_cluster = 12, cap = 2, clusters_on_y = TRUE)
  
  # Cluster Detail Plots
  ids_cur <- as.character(intersect(DX$participant_id, clusters_all$participant_id))
  cl_map <- setNames(clusters_all$cluster, clusters_all$participant_id)
  cl_vec <- as.integer(cl_map[ids_cur])
  cl_ids <- sort(unique(cl_vec[cl_vec != 0]))
  
  plot_cluster_detail_one <- function(cid, tab_cluster, ids_cur, cl_vec, top_n = 15) {
    T <- tab_cluster %>% dplyr::filter(cluster == cid) %>% dplyr::mutate(log2_lift = if (!"log2_lift" %in% names(.)) log2(pmax(lift, 1e-12)) else log2_lift) %>% dplyr::arrange(dplyr::desc(log2_lift)) %>% dplyr::filter(is.finite(in_prev), is.finite(out_prev))
    if (!nrow(T)) return(invisible(NULL))
    
    n_in <- sum(cl_vec == cid); n_out <- length(ids_cur) - n_in
    k_in <- round(pmin(pmax(T$in_prev, 0), 1) * n_in)
    k_out <- round(pmin(pmax(T$out_prev, 0), 1) * n_out)
    T$dx <- factor(T$dx, levels = rev(head(T$dx, top_n)))
    T <- subset(T, dx %in% levels(T$dx))
    ci_in <- wilson_ci(k_in, n_in); ci_out <- wilson_ci(k_out, n_out)
    
    D <- tibble::tibble(dx = rep(T$dx, each = 2), group = rep(c("in-cluster", "out-of-cluster"), times = nrow(T)), prev = as.numeric(t(cbind(T$in_prev, T$out_prev))), lo = as.numeric(t(cbind(ci_in[, "lo"], ci_out[, "lo"]))), hi = as.numeric(t(cbind(ci_in[, "hi"], ci_out[, "hi"]))))
    S <- tibble::tibble(dx = T$dx, prev_in = T$in_prev, prev_out = T$out_prev)
    
    p_det <- ggplot2::ggplot(D, ggplot2::aes(y = dx)) +
      ggplot2::geom_segment(data = S, ggplot2::aes(y = dx, yend = dx, x = prev_out, xend = prev_in), linewidth = 0.5, alpha = 0.6, colour = "grey40", inherit.aes = FALSE) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, xmax = hi, colour = group), height = 0, alpha = 0.8) +
      ggplot2::geom_point(ggplot2::aes(x = prev, colour = group), size = 2.2) +
      ggplot2::scale_x_continuous(labels = scales::percent, limits = c(0, 1)) +
      ggplot2::scale_colour_manual(values = c("out-of-cluster" = "#606060", "in-cluster" = "#1b7cff"), name = NULL) +
      ggplot2::labs(x = "Prevalence", y = NULL, title = paste0("Cluster C", cid, ": Dx enrichment (top ", min(top_n, nrow(T)), ")")) +
      ggplot2::theme_minimal(11) + ggplot2::theme(legend.position = "top", panel.grid.minor = ggplot2::element_blank())
    
    save_plot_gg(sprintf("FIG_dx_overrep_C%d_detail", cid), p_det, width = 7.5, height = 6.5)
    invisible(p_det)
  }
  
  invisible(lapply(cl_ids, function(cid) { plot_cluster_detail_one(cid, tab_sig_maj, ids_cur, cl_vec, top_n = 15) }))
  
  # MDS Validation & CV R2
  cv_r2_mds <- function(D, k_grid = 1:8, K = 5, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    Dm <- as.matrix(D); n <- nrow(Dm); folds <- sample(rep(1:K, length.out = n))
    res <- vector("list", length(k_grid))
    for (kk in seq_along(k_grid)) {
      k <- k_grid[kk]; r2_fold <- rep(NA_real_, K)
      for (f in 1:K) {
        te <- which(folds == f); tr <- setdiff(seq_len(n), te)
        fit <- stats::cmdscale(as.dist(Dm[tr, tr]), k = k, eig = TRUE, add = TRUE)
        Xtr <- as.matrix(fit$points)
        Dte2 <- Dm[te, tr, drop = FALSE]^2; Dtr2 <- Dm[tr, tr, drop = FALSE]^2
        ntr <- length(tr)
        Htr <- diag(ntr) - matrix(1, ntr, ntr) / ntr
        Btr <- -0.5 * Htr %*% Dtr2 %*% Htr
        E <- eigen(Btr, symmetric = TRUE)
        pos <- which(E$values > 1e-12); m <- min(k, length(pos))
        vals <- E$values[pos][1:m]; V <- E$vectors[, pos, drop = FALSE][, 1:m, drop = FALSE]
        cmtr <- colMeans(Dtr2); rmte <- rowMeans(Dte2); mu <- mean(Dtr2)
        Gte <- -0.5 * (Dte2 - matrix(rmte, nrow(Dte2), ncol(Dte2), byrow = FALSE) - matrix(cmtr, nrow(Dte2), ncol(Dte2), byrow = TRUE) + mu)
        Xte <- sweep(Gte %*% V, 2, 1 / sqrt(vals), `*`)
        X <- rbind(Xtr, Xte); idx <- c(tr, te)
        dhat <- as.numeric(dist(X)); dref <- as.numeric(as.dist(Dm[idx, idx, drop = FALSE]))
        r2_fold[f] <- suppressWarnings(cor(dhat, dref)^2)
      }
      res[[kk]] <- data.frame(k = k, r2 = mean(r2_fold, na.rm = TRUE))
    }
    do.call(rbind, res)
  }
  
  fit2 <- stats::cmdscale(D_dx, k = 2, eig = TRUE, add = TRUE)
  XY <- fit2$points; colnames(XY) <- c("MDS1", "MDS2")
  dfp <- data.frame(MDS1 = XY[, 1], MDS2 = XY[, 2], cluster = factor(memb2))
  dfp <- subset(dfp, cluster != 0)
  dfp$cluster_f <- factor(paste0("C", memb2[memb2 != 0]), levels = paste0("C", sort(unique(memb2[memb2 != 0]))))
  p2 <- ggplot2::ggplot(dfp, ggplot2::aes(MDS1, MDS2, colour = cluster_f)) +
    ggplot2::geom_point(size = 1.9, alpha = 0.95) +
    ggplot2::coord_equal() +
    ggplot2::scale_colour_manual(values = cluster_colours(levels(dfp$cluster_f)))
  save_plot_gg("FIG_dxspace_clusters_mds2", p2, width = 7, height = 6, dpi = 150)
  
  r2k <- cv_r2_mds(D_dx, k_grid = 1:8, K = 5, seed = 42)
  p_r2 <- ggplot2::ggplot(r2k, ggplot2::aes(k, r2)) +
    ggplot2::geom_line() + ggplot2::geom_point() + ggplot2::scale_x_continuous(breaks = r2k$k) +
    ggplot2::labs(x = "MDS dimensionality (k)", y = expression(R^2 ~ "(distances, 5-fold CV)"), title = "Out-of-sample distance preservation vs k") +
    ggplot2::theme_minimal(12)
  save_plot_gg("FIG_mds_cv_r2", p_r2, width = 6.5, height = 4.5, dpi = 150)
  
  n_valid <- max(stab$trials)
  stab$dx_lab <- ifelse(stab$trials < n_valid, sprintf("%s (%d/%d)", stab$dx, stab$trials, n_valid), stab$dx)
  p_ci <- ggplot2::ggplot(stab, ggplot2::aes(x = reorder(dx_lab, freq), y = freq)) +
    ggplot2::geom_col(fill = "grey30") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = lo, ymax = hi), width = .2) +
    ggplot2::coord_flip() +
    ggplot2::scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    ggplot2::labs(x = NULL, y = "bootstrap frequency (among eligible resamples)", title = "Major Dx stability (per-dx denominator)") +
    ggplot2::theme_minimal(12)
  save_plot_gg("FIG_major_dx_bootstrap_frequency_CI_by_eligibility", p_ci, width = 7, height = 5, dpi = 150)
}