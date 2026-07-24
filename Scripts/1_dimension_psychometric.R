# ==============================================================================
#                               dimension_psychometric.R
# ==============================================================================

# ==============================================================================
# 1. Helper Functions
# ==============================================================================

if (!BASE_DIM %in% c(2L, 3L)) {
  stop("BASE_DIM must be 2 or 3, got: ", BASE_DIM)
}

# --- Data Preparation for Gower Distance ---
.robust01 <- function(x, q = c(0.01, 0.99)) {
  x <- as.numeric(x)
  qs <- stats::quantile(x, probs = q, na.rm = TRUE, type = 8)
  lo <- qs[[1]]; hi <- qs[[2]]
  if (!is.finite(lo) || !is.finite(hi) || hi <= lo) return(x)  # fallback
  xw <- pmin(pmax(x, lo), hi)
  (xw - lo) / (hi - lo)
}

# winsorise + quantile bins -> numeric in [0,1] + tiny jitter (NOT factor)
.bin_continuous_numeric <- function(x,
                                    bins = 8L,
                                    q_winsor = c(0.05, 0.95),
                                    min_unique = 4L,
                                    tie_jitter = 1e-9,
                                    out_jitter = 1e-6,
                                    seed = NULL) {
  
  x <- as.numeric(x)
  if (!any(is.finite(x))) return(x)
  
  qs <- stats::quantile(x, probs = q_winsor, na.rm = TRUE, type = 8)
  lo <- qs[[1]]; hi <- qs[[2]]
  if (is.finite(lo) && is.finite(hi) && hi > lo) x <- pmin(pmax(x, lo), hi)
  
  u <- sort(unique(na.omit(x)))
  if (length(u) < min_unique) {
    lev <- u
    idx <- match(x, lev)
    return((idx - 1) / max(1, length(lev) - 1))
  }
  
  sdv <- stats::sd(x, na.rm = TRUE)
  
  xj <- if (is.finite(sdv) && sdv > 0) {
    eps <- tie_jitter * sdv
    if (!is.null(seed)) {
      x + .with_seed(seed, stats::rnorm(length(x), 0, eps))
    } else {
      x + stats::rnorm(length(x), 0, eps)
    }
  } else x
  
  probs <- seq(0, 1, length.out = as.integer(bins) + 1L)
  br <- stats::quantile(xj, probs = probs, na.rm = TRUE, type = 8)
  br <- unique(as.numeric(br))
  
  if (length(br) < 3L) return(.robust01(x, q = q_winsor))
  
  f <- cut(xj, breaks = br, include.lowest = TRUE, right = TRUE)
  k <- as.integer(f)
  B_eff <- max(k, na.rm = TRUE)
  z <- (k - 1) / max(1, B_eff - 1)
  z[is.na(k)] <- NA_real_
  
  step <- 1 / max(1, B_eff - 1)
  j_sd <- out_jitter * step
  if (is.finite(j_sd) && j_sd > 0) {
    if (!is.null(seed)) {
      z <- z + .with_seed(as.integer(seed + 1L), stats::rnorm(length(z), 0, j_sd))
    } else {
      z <- z + stats::rnorm(length(z), 0, j_sd)
    }
  }
  
  z <- pmin(pmax(z, 0), 1)
  z
}

dist_health_console <- function(D, name = "D", eps0 = 1e-12) {
  M <- as.matrix(D)
  diag(M) <- Inf
  n <- nrow(M)
  
  # first and second NN
  negM <- -M
  i1 <- max.col(negM, ties.method = "first")
  d1 <- M[cbind(seq_len(n), i1)]
  negM[cbind(seq_len(n), i1)] <- -Inf
  i2 <- max.col(negM, ties.method = "first")
  d2 <- M[cbind(seq_len(n), i2)]
  
  # basic scale & tie diagnostics
  v <- as.numeric(D)
  v <- v[is.finite(v)]
  uq <- length(unique(round(v, 10)))
  prop0 <- mean(v <= eps0)
  
  # near-collision mass
  qd1 <- quantile(d1[is.finite(d1)], probs = c(0, 0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 1), na.rm = TRUE)
  eps_adapt <- max(1e-12, as.numeric(qd1[[2]]) * 0.1)
  
  d1c <- pmax(d1, eps_adapt)
  d2c <- pmax(d2, d1c + eps_adapt)
  r <- d2c / d1c
  lr <- log(r[is.finite(r) & r > 1])
  
  cat("\n================ DIST HEALTH:", name, "================\n")
  cat(sprintf("n=%d | unique_dist(rounded 1e-10)=%d | prop(dist<=1e-12)=%.4f\n", n, uq, prop0))
  cat("d1 quantiles:\n"); print(qd1)
  cat(sprintf("eps_adapt=%.3g | mean(log r)=%.3f | sd(log r)=%.3f | TwoNN=%.3f\n",
              eps_adapt,
              mean(lr, na.rm=TRUE), sd(lr, na.rm=TRUE),
              1/mean(lr, na.rm=TRUE)))
  cat(sprintf("tie rate (d1 <= eps_adapt): %.4f\n", mean(d1 <= eps_adapt, na.rm=TRUE)))
}

knn_from_dist <- function(D, k = 15) {
  M <- as.matrix(D)
  diag(M) <- Inf
  t(apply(M, 1, function(r) order(r)[1:k]))
}

mean_jaccard_knn <- function(idxA, idxB) {
  stopifnot(nrow(idxA) == nrow(idxB), ncol(idxA) == ncol(idxB))
  n <- nrow(idxA); k <- ncol(idxA)
  s <- 0
  for (i in 1:n) {
    a <- idxA[i, ]; b <- idxB[i, ]
    s <- s + length(intersect(a, b)) / length(union(a, b))
  }
  s / n
}

constant_profile_value <- function(Xdf) {
  n <- nrow(Xdf)
  if (!n || !ncol(Xdf)) return(rep(NA_character_, n))
  
  ref <- as.character(Xdf[[1L]])
  same <- rep(TRUE, n)
  if (ncol(Xdf) > 1L) {
    for (j in 2:ncol(Xdf)) {
      equal_j <- as.character(Xdf[[j]]) == ref
      equal_j[is.na(equal_j)] <- FALSE
      same <- same & equal_j
    }
  }
  
  out <- rep(NA_character_, n)
  out[same] <- ref[same]
  out
}

warn_diag_subsample <- function(tag, n_used, n_total, cap, pool_label = "rows") {
  warning(sprintf(
    "[%s] Using %d/%d %s because full distance diagnostics materialise dist objects and dense matrices (O(n^2) memory/time); cap=%d.",
    tag, n_used, n_total, pool_label, cap
  ))
}

prep_X_for_gower <- function(X,
                             rare_prop = 0.01,
                             do_jitter = TRUE,
                             robust_cont = TRUE,
                             prep_shcnt = SHT_CONT,
                             prep_cntbins = CONT_BINS,
                             prep_cntwinsor_q = CONT_WINSOR_Q,
                             seed = NULL,
                             treat_ordinals_as_nominal = TREAT_ORDINALS_AS_NOMINAL) {
  
  X1 <- as.data.frame(X, check.names = TRUE, stringsAsFactors = FALSE)
  binned_continuous_cols <- character(0)
  for (nm in names(X1)) if (is.character(X1[[nm]])) X1[[nm]] <- factor(X1[[nm]])
  
  drop_rare <- function(f, prop) {
    if (!is.factor(f) || is.ordered(f)) return(f)
    tb <- prop.table(table(f))
    keep <- names(tb)[tb >= prop]
    out <- as.character(f)
    out[!is.na(out) & !out %in% keep] <- NA_character_
    f <- factor(out, levels = keep, exclude = NA)
    droplevels(f)
  }
  X1 <- as.data.frame(lapply(X1, drop_rare, prop = rare_prop), stringsAsFactors = FALSE)
  
  if (isTRUE(robust_cont) && !isTRUE(prep_shcnt)) {
    for (nm in names(X1)) if (is.numeric(X1[[nm]])) X1[[nm]] <- .robust01(X1[[nm]])
  }
  
  if (isTRUE(prep_shcnt)) do_jitter <- FALSE
  
  if (isTRUE(do_jitter)) {
    for (nm in names(X1)) {
      if (is.numeric(X1[[nm]])) {
        sdv <- stats::sd(X1[[nm]], na.rm = TRUE)
        if (is.finite(sdv) && sdv > 0) {
          eps <- 1e-6 * sdv
          if (!is.null(seed)) {
            s_nm <- .seed_from_key(seed, paste0("prep_jitter::", nm))
            X1[[nm]] <- X1[[nm]] + .with_seed(s_nm, stats::rnorm(length(X1[[nm]]), 0, eps))
          } else {
            X1[[nm]] <- X1[[nm]] + stats::rnorm(length(X1[[nm]]), 0, eps)
          }
        }
      }
    }
  }
  
  if (isTRUE(prep_shcnt)) {
    for (nm in names(X1)) {
      if (is.numeric(X1[[nm]])) {
        s_nm <- if (!is.null(seed)) .seed_from_key(seed, paste0("bin::", nm)) else NULL
        
        if (isTRUE(CONT_AS_ORD)) {
          z <- .bin_continuous_numeric(
            X1[[nm]],
            bins       = prep_cntbins,
            q_winsor   = prep_cntwinsor_q,
            out_jitter = 0,
            seed       = s_nm
          )
          step <- 1 / max(1, prep_cntbins - 1)
          k <- round(z / step) + 1L
          X1[[nm]] <- ordered(k)
          binned_continuous_cols <- union(binned_continuous_cols, nm)
        } else {
          z <- .bin_continuous_numeric(
            X1[[nm]],
            bins       = prep_cntbins,
            q_winsor   = prep_cntwinsor_q,
            out_jitter = CONT_JITTER_FRAC,
            seed       = s_nm
          )
          X1[[nm]] <- as.numeric(z)
        }
      }
    }
  }
  
  if (isTRUE(treat_ordinals_as_nominal)) {
    treat_binned_cont_as_nominal <- TREAT_BINNED_CONTINUOUS_AS_NOMINAL
    X1 <- as.data.frame(Map(function(v, nm) {
      is_binned_cont <- nm %in% binned_continuous_cols
      if (is.ordered(v) && (isTRUE(treat_binned_cont_as_nominal) || !is_binned_cont)) {
        factor(v, levels = levels(v), ordered = FALSE, exclude = NULL)
      } else {
        v
      }
    }, X1, names(X1)), stringsAsFactors = FALSE)
  }
  
  ord_cols <- names(X1)[vapply(X1, is.ordered, logical(1))]
  fac_cols <- names(X1)[vapply(X1, function(z) is.factor(z) && !is.ordered(z), logical(1))]
  .is_binary <- function(x) length(unique(na.omit(x))) == 2
  bin_cols <- if (isTRUE(treat_ordinals_as_nominal)) {
    character(0)
  } else {
    fac_cols[vapply(X1[fac_cols], .is_binary, logical(1))]
  }
  
  type_list <- list()
  if (length(bin_cols)) type_list$asymm <- bin_cols
  if (length(ord_cols)) type_list$ordratio <- ord_cols
  
  w <- setNames(rep(1, ncol(X1)), names(X1))
  list(
    X = X1,
    type = type_list,
    weights = w,
    binned_continuous_cols = binned_continuous_cols
  )
}

# --- Gower Distance Wrapper ---
gower_dist <- function(Xdf, type_list = NULL, weights = NULL) {
  if (!is.null(type_list)) {
    type_list <- lapply(type_list, function(cols) intersect(cols, names(Xdf)))
    type_list <- type_list[lengths(type_list) > 0]
    if (!length(type_list)) type_list <- NULL
  }
  
  if (is.null(weights)) {
    weights <- rep(1, ncol(Xdf))
  } else if (length(weights) == 1) {
    weights <- rep(weights, ncol(Xdf))
  } else if (!is.null(names(weights))) {
    weights <- weights[names(Xdf)]
  }
  
  stopifnot(length(weights) == ncol(Xdf))
  cluster::daisy(Xdf, metric = "gower", type = type_list, weights = weights)
}

# --- Nearest Neighbor Helpers (Vectorized/Fast) ---

# Fast single pass: first and second NN using max.col
.two_nn_from_distvec <- function(D, eps = NULL) {
  M <- as.matrix(D)
  n <- nrow(M)
  if (n < 2L) return(list(d1 = numeric(0), d2 = numeric(0)))
  diag(M) <- Inf
  
  negM <- -M
  idx1 <- max.col(negM, ties.method = "first")
  d1 <- M[cbind(seq_len(n), idx1)]
  
  negM[cbind(seq_len(n), idx1)] <- -Inf
  idx2 <- max.col(negM, ties.method = "first")
  d2 <- M[cbind(seq_len(n), idx2)]
  
  if (is.null(eps)) {
    d1f <- d1[is.finite(d1)]
    eps <- if (length(d1f)) max(1e-6, as.numeric(stats::quantile(d1f, 0.05, na.rm = TRUE)) * 0.05) else 1e-6
  }
  
  d1 <- pmax(d1, eps)
  d2 <- pmax(d2, d1 + eps)
  list(d1 = d1, d2 = d2)
}

# TwoNN intrinsic dimension without materialising the full matrix
twonn_id_from_dist <- function(D, eps = 1e-8, trim = 0.02) {
  n <- attr(D, "Size")
  if (n < 3L) return(NA_real_)
  
  nn <- .two_nn_from_distvec(D, eps)
  r  <- nn$d2 / nn$d1
  r  <- r[is.finite(r) & r > 1]
  
  if (!length(r)) return(NA_real_)
  
  lr <- sort(log(r))
  k  <- floor(trim * length(lr))
  if (k > 0 && 2 * k < length(lr)) lr <- lr[(k + 1):(length(lr) - k)]
  
  1 / mean(lr)
}

first_nn_d1 <- function(D) {
  n <- attr(D, "Size")
  if (n < 2L) return(rep(Inf, n))
  .two_nn_from_distvec(D)$d1
}

# --- Dedup & Coreset Helpers ---

collapse_curve <- function(D, eps_grid) {
  n0 <- attr(D, "Size")
  
  get_ngroups <- function(eps) {
    length(unique(hclust(D, method = "complete") |> cutree(h = eps)))
  }
  
  data.frame(
    eps = eps_grid,
    n_groups = sapply(eps_grid, get_ngroups),
    prop_retained = NA_real_
  ) |>
    transform(prop_retained = n_groups / n0)
}

complete_groups <- function(D, eps) {
  hclust(D, method = "complete") |> cutree(h = eps)
}

dedup_key_col <- function(v, mode = c("hash_exact", "hash_round"),
                          digits = 6L, na_token = "<NA>") {
  mode <- match.arg(mode)
  
  if (is.numeric(v)) {
    x <- as.numeric(v)
    
    if (mode == "hash_exact") {
      out <- ifelse(is.na(x), na_token, sprintf("%.17g", x))
      out[out %in% c("-0", "-0.0")] <- "0"
      return(out)
    }
    
    x <- round(x, digits = digits)
    fmt <- paste0("%.", digits, "f")
    out <- ifelse(is.na(x), na_token, sprintf(fmt, x))
    
    zero_tok <- sprintf(fmt, 0)
    out[grepl("^-0(?:\\.0+)?$", out)] <- zero_tok
    return(out)
  }
  
  if (is.factor(v) || is.ordered(v) || is.character(v) || is.logical(v)) {
    out <- as.character(v)
    out[is.na(out)] <- na_token
    return(out)
  }
  
  out <- as.character(v)
  out[is.na(out)] <- na_token
  out
}

dedup_groups_from_hash <- function(Xdf, mode = c("hash_exact", "hash_round"),
                                   digits = 6L, na_token = "<NA>") {
  mode <- match.arg(mode)
  
  X_key <- as.data.frame(
    lapply(Xdf, dedup_key_col, mode = mode, digits = digits, na_token = na_token),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  key <- do.call(paste, c(X_key, sep = "\r"))
  split(seq_len(nrow(Xdf)), key, drop = TRUE)
}

# Core-set finder without matrix expansion
twonn_core_by_slope <- function(D, min_frac = 0.30, w = 20, slope_tol = 0.08, rmse_tol = 0.10) {
  n <- attr(D, "Size")
  if (n < 8L) return(seq_len(n))
  
  nn <- .two_nn_from_distvec(D, eps = .Machine$double.eps)
  mu <- pmax(nn$d2 / pmax(nn$d1, .Machine$double.eps), 1 + 1e-12)
  
  ord <- order(mu)
  x <- log(mu[ord])
  m <- length(x)
  
  u <- (seq_len(m) - 0.5) / (m + 1)
  y <- log(1 - u)
  
  k0 <- max(20L, floor(min_frac * m))
  slope <- rmse <- rep(NA_real_, m)
  
  for (k in k0:(m - 2L)) {
    fit <- stats::lm(y[1:k] ~ x[1:k])
    slope[k] <- coef(fit)[2]
    rmse[k] <- sqrt(mean(residuals(fit)^2))
  }
  
  ok <- which(is.finite(slope) & is.finite(rmse))
  if (!length(ok)) return(ord[seq_len(max(3L, k0))])
  
  pick <- function(k) {
    L <- max(k0, k - w + 1L)
    s <- slope[L:k]
    (max(s, na.rm = TRUE) - min(s, na.rm = TRUE) <= slope_tol) && (rmse[k] <= rmse_tol)
  }
  
  ks <- ok[vapply(ok, pick, logical(1))]
  k_star <- if (length(ks)) max(ks) else floor(0.6 * m)
  
  ord[seq_len(max(3L, min(k_star, m - 2L)))]
}

core_band_idx <- function(D, k = 10, band = c(0.20, 0.70)) {
  M <- as.matrix(D)
  diag(M) <- Inf
  
  kth <- function(r, k) {
    rf <- r[is.finite(r)]
    if (!length(rf)) return(NA_real_)
    k_eff <- min(k, length(rf))
    sort(rf, partial = k_eff)[k_eff]
  }
  
  rk <- apply(M, 1, kth, k = k)
  ok <- is.finite(rk)
  if (!any(ok)) return(integer(0))
  
  q <- stats::quantile(rk[ok], band, na.rm = TRUE)
  which(ok & rk >= q[1] & rk <= q[2])
}

# TwoNN + LB estimators
lb_mle_id <- function(Dm, k_lo = 5, k_hi = 15) {
  Dm <- as.matrix(Dm)
  n <- nrow(Dm)
  diag(Dm) <- Inf
  if (n <= k_lo) return(NA_real_)
  
  k_hi <- max(k_lo, min(k_hi, n - 1))
  
  ids <- sapply(k_lo:k_hi, function(k) {
    nn <- t(apply(Dm, 1L, function(r) {
      rf <- r[is.finite(r)]
      m <- length(rf)
      if (m < k) return(rep(NA_real_, k))
      sort(rf, partial = k)[1:k]
    }))
    
    if (!nrow(nn)) return(NA_real_)
    
    l <- log(nn[, k, drop = TRUE] / nn[, 1:(k - 1), drop = FALSE])
    d <- 1 / rowMeans(l, na.rm = TRUE)
    mean(d[is.finite(d)], na.rm = TRUE)
  })
  
  mean(ids, na.rm = TRUE)
}

# --- Weight Optimization Helpers ---

make_NS_cache <- function(Xdf, type = NULL) {
  p <- ncol(Xdf)
  N_list <- vector("list", p)
  S_list <- vector("list", p)
  
  for (j in seq_len(p)) {
    w1 <- rep(0, p)
    w1[j] <- 1
    Dj <- cluster::daisy(Xdf, metric = "gower", type = type, weights = w1)
    vD <- as.numeric(Dj)
    ok <- as.numeric(is.finite(vD))
    vD[!is.finite(vD)] <- 0
    N_list[[j]] <- vD
    S_list[[j]] <- ok
  }
  list(N = N_list, S = S_list, n = nrow(Xdf))
}

gower_subset_support <- function(Xdf,
                                 min_prop = GOWER_SUBSET_MIN_LEVEL_PROP,
                                 min_n = GOWER_SUBSET_MIN_LEVEL_N,
                                 binned_continuous_cols = character(0),
                                 min_binned_unique = 2L) {
  n <- nrow(Xdf)
  if (!ncol(Xdf)) {
    stop("[optim] subset support check got 0 columns; check input delimiter and predictor selection.")
  }
  min_count <- max(as.integer(min_n), ceiling(min_prop * n))
  
  rows <- lapply(names(Xdf), function(nm) {
    v <- Xdf[[nm]]
    
    if (is.numeric(v)) {
      ok <- is.finite(v)
      n_observed <- sum(ok)
      vv <- v[ok]
      n_unique <- length(unique(vv))
      tab <- table(vv)
      is_low_cardinality_numeric <- n_unique > 0L && n_unique <= LOW_SUPPORT_NUMERIC_MAX_UNIQUE
      min_level_count <- if (is_low_cardinality_numeric) {
        min(as.integer(tab))
      } else if (n_unique >= 2L) {
        NA_integer_
      } else {
        n_observed
      }
      support_ok <- if (is_low_cardinality_numeric) {
        n_unique >= 2L && min_level_count >= min_count
      } else {
        n_unique >= 2L
      }
      support_rule <- if (is_low_cardinality_numeric) "numeric_low_cardinality" else "numeric_continuous"
    } else {
      vv <- v[!is.na(v)]
      n_observed <- length(vv)
      tab <- table(vv)
      n_unique <- length(tab)
      min_level_count <- if (n_unique) min(as.integer(tab)) else 0L
      is_binned_cont <- nm %in% binned_continuous_cols
      if (is_binned_cont) {
        support_ok <- n_unique >= as.integer(min_binned_unique)
        support_rule <- "binned_continuous_ordered"
      } else {
        support_ok <- n_unique >= 2L && min_level_count >= min_count
        support_rule <- if (is.ordered(v)) "ordered" else "categorical"
      }
    }
    
    data.frame(
      var = nm,
      support_rule = support_rule,
      n_observed = as.integer(n_observed),
      n_unique = as.integer(n_unique),
      min_level_count = as.integer(ifelse(is.na(min_level_count), NA_integer_, min_level_count)),
      min_count_required = as.integer(min_count),
      support_ok = isTRUE(support_ok),
      stringsAsFactors = FALSE
    )
  })
  
  out <- do.call(rbind, rows)
  if (is.null(out) || !nrow(out) || !"support_ok" %in% names(out) || !is.logical(out$support_ok)) {
    stop("[optim] subset support check did not produce a logical support_ok column.")
  }
  out
}

optimise_gower_weights_constrained <- function(
    X, init_weights, allow_update,
    objective    = "TwoNN_all",
    w_min        = W_MIN,
    step_grid    = W_STEP_GRID,
    batch_k      = W_BATCH_K,
    batch_factor = W_BATCH_FACTOR,
    max_iter     = W_MAX_ITERS,
    n_rows_sub   = N_ROWS_SUB,
    ncores       = NULL,
    seed_jitter  = SEED_JITTER,
    reps_idx     = NULL,
    core_idx_rep = NULL,
    lambda_l1    = GOWER_L1,
    lambda_l2    = GOWER_L2,
    verbose      = TRUE,
    plot_progress = TRUE,
    base_seed    = SEED_GLOBAL,
    progress_fun = NULL
) {
  seed_sub  <- .seed_from_key(base_seed, "optim::subsample_rows")
  seed_iter <- .seed_from_key(base_seed, "optim::iter_sampling")
  seed_prep <- seed_jitter
  
  # Internal ID calc from numerator/denominator
  calc_id_fast <- function(num, den, n_rows) {
    Dvec <- num / pmax(den, .Machine$double.eps)
    attr(Dvec, "Size")  <- n_rows
    attr(Dvec, "Diag")  <- FALSE
    attr(Dvec, "Upper") <- FALSE
    class(Dvec) <- "dist"
    twonn_id_from_dist(Dvec)
  }
  
  # L1 + L2 penalty centred at w_min, only for excess weight above the floor
  penalty_val <- function(w) {
    v <- pmax(w - w_min, 0)  # excess over floor
    p1 <- if (lambda_l1 > 0) lambda_l1 * mean(v)   else 0
    p2 <- if (lambda_l2 > 0) lambda_l2 * mean(v^2) else 0
    p1 + p2
  }
  
  # Combine ID and penalty into a single scalar
  calc_obj_fast <- function(num, den, n_rows, w) {
    id_val <- calc_id_fast(num, den, n_rows)
    if (!is.finite(id_val)) id_val <- Inf
    obj_val <- id_val + penalty_val(w)
    list(id = id_val, obj = obj_val)
  }
  
  # Prepare subset
  px  <- prep_X_for_gower(X, rare_prop = RARE_LEVEL_MIN_PROP, do_jitter = TRUE, seed = seed_prep)
  X0  <- px$X
  typ <- px$type
  
  row_pool <- if (!is.null(reps_idx)) reps_idx else seq_len(nrow(X0))
  if (is.null(n_rows_sub)) {
    eff_n_sub <- length(row_pool)
  } else if (!is.finite(n_rows_sub) || as.integer(n_rows_sub) <= 0L) {
    # Explicit full-pool request (e.g. n_rows_sub = Inf)
    eff_n_sub <- length(row_pool)
  } else {
    eff_n_sub <- min(as.integer(n_rows_sub), length(row_pool))
  }
  
  ix_sub_from_reps <- if (length(row_pool) > eff_n_sub) {
    if (isTRUE(FIX_REP_SUBSET)) {
      head(row_pool, eff_n_sub)
    } else {
      .with_seed(seed_sub, sample(row_pool, eff_n_sub))
    }
  } else {
    row_pool
  }
  
  Xs    <- X0[ix_sub_from_reps, , drop = FALSE]
  n_sub <- nrow(Xs)
  
  vars <- colnames(Xs)
  w    <- init_weights[vars]
  allow_update <- allow_update[vars]
  
  subset_support <- gower_subset_support(
    Xs,
    binned_continuous_cols = px$binned_continuous_cols
  )
  if (!is.data.frame(subset_support) || !"var" %in% names(subset_support) || !"support_ok" %in% names(subset_support)) {
    stop("[optim] subset support check returned an invalid table.")
  }
  if (!is.logical(subset_support$support_ok)) {
    stop("[optim] subset support_ok is not logical.")
  }
  unsupported <- subset_support$var[!subset_support$support_ok]
  if (isTRUE(GOWER_REQUIRE_SUBSET_SUPPORT) && length(unsupported)) {
    w[unsupported] <- w_min
    allow_update[unsupported] <- FALSE
    if (verbose) {
      cat(sprintf(
        "[optim] %d variable(s) lack subset support and were put at w_min for this run.\n",
        length(unsupported)
      ))
    }
  }
  w[!allow_update] <- pmax(w_min, w[!allow_update])
  
  # Build cache
  cache <- make_NS_cache(Xs, type = typ)
  
  # Initialise state: current numerator/denominator and objective
  num_cur <- Reduce(`+`, Map(`*`, cache$N, as.list(w)))
  den_cur <- Reduce(`+`, Map(`*`, cache$S, as.list(w)))
  
  o0   <- calc_obj_fast(num_cur, den_cur, n_sub, w)
  id0  <- o0$id
  obj0 <- o0$obj
  
  hist <- data.frame(
    iter    = 0L,
    ID      = id0,
    Obj     = obj0,
    changed = NA_character_,
    note    = NA_character_,
    stringsAsFactors = FALSE
  )
  
  if (verbose) {
    cat(sprintf(
      "[optim] Start ID: %.3f | L1=%.3f L2=%.3f | N_sub=%d | Mode: STOCHASTIC\n",
      id0, lambda_l1, lambda_l2, n_sub
    ))
  }
  
  id      <- id0
  obj_cur <- obj0
  
  max_iter_eff <- if (is.null(max_iter) || !is.finite(max_iter)) 1000L else as.integer(max_iter)
  
  # Optimisation config
  N_SAMPLE_PER_ITER <- 50
  step_grid <- sort(step_grid, decreasing = FALSE)
  hot_vars  <- integer(0)  # momentum
  
  for (it in seq_len(max_iter_eff)) {
    can_all <- which(allow_update & (w > w_min + 1e-12))
    if (!length(can_all)) break
    
    # 1. Stochastic sampling: keep hot vars, fill rest with randoms
    n_rnd    <- max(10, N_SAMPLE_PER_ITER - length(hot_vars))
    rnd_vars <- .with_seed(
      as.integer(seed_iter + it),
      sample(can_all, min(length(can_all), n_rnd))
    )
    can_iter <- unique(c(hot_vars, rnd_vars))
    
    # 2. Evaluate candidates (greedy)
    best_res <- list(obj = Inf)
    
    for (j in can_iter) {
      w_base <- w[j]
      for (factor in step_grid) {
        w_new <- max(w_min, w_base * factor)
        if (w_new >= w_base - 1e-12) next
        
        delta   <- w_new - w_base
        num_try <- num_cur + delta * cache$N[[j]]
        den_try <- den_cur + delta * cache$S[[j]]
        
        w_try      <- w
        w_try[j]   <- w_new
        o_try      <- calc_obj_fast(num_try, den_try, n_sub, w_try)
        
        if (is.finite(o_try$obj) && o_try$obj < best_res$obj - 1e-9) {
          best_res <- list(
            obj   = o_try$obj,
            id    = o_try$id,
            j     = j,
            w_new = w_new,
            num   = num_try,
            den   = den_try
          )
        }
      }
    }
    
    changed <- FALSE
    if (is.finite(best_res$obj) && best_res$obj < obj_cur - 1e-6) {
      jbest <- best_res$j
      wbest <- best_res$w_new
      
      w[jbest] <- wbest
      num_cur  <- best_res$num
      den_cur  <- best_res$den
      id       <- best_res$id
      obj_cur  <- best_res$obj
      changed  <- TRUE
      
      hot_vars <- unique(c(jbest, hot_vars))
      if (length(hot_vars) > 10) hot_vars <- head(hot_vars, 10)
      
      hist <- rbind(
        hist,
        data.frame(
          iter    = it,
          ID      = id,
          Obj     = obj_cur,
          changed = vars[jbest],
          note    = sprintf("%.3f", wbest),
          stringsAsFactors = FALSE
        )
      )
      
      if (verbose) {
        cat(sprintf(
          "   iter %d: %s -> %.3f (ID: %.3f, Obj: %.3f)\n",
          it, vars[jbest], wbest, id, obj_cur
        ))
      }
    } else {
      # No improvement: flush momentum or stop
      if (length(hot_vars) > 0) {
        hot_vars <- integer(0)
        if (verbose) cat("   [optim] Momentum lost, flushing hot vars.\n")
      } else {
        if (verbose) cat("[optim] No improvement in random subset.\n")
        break
      }
    }
    
    if (!is.null(progress_fun) && (it == 1L || it %% 10L == 0L)) {
      progress_fun(list(
        iter = it,
        ID   = id,
        Obj  = obj_cur
      ))
    }
    
    # 3. Batch descent (every 5th iter), using penalised objective as well
    if (batch_k > 1 && (it %% 5 == 0)) {
      remain <- setdiff(can_all, best_res$j)
      remain <- .with_seed(
        as.integer(seed_iter + 10000L + it),
        sample(remain, min(length(remain), N_SAMPLE_PER_ITER))
      )
      
      if (length(remain) > 0) {
        scores <- numeric(length(remain))
        for (i in seq_along(remain)) {
          j <- remain[i]
          delta_b <- max(w_min, w[j] * batch_factor) - w[j]
          if (abs(delta_b) < 1e-12) {
            scores[i] <- Inf
          } else {
            num_try <- num_cur + delta_b * cache$N[[j]]
            den_try <- den_cur + delta_b * cache$S[[j]]
            w_try   <- w
            w_try[j] <- w[j] + delta_b
            scores[i] <- calc_obj_fast(num_try, den_try, n_sub, w_try)$obj
          }
        }
        
        ord       <- order(scores)
        take_idx  <- head(ord, min(batch_k, length(ord)))
        take_vars <- remain[take_idx]
        
        if (length(take_vars)) {
          num_b <- num_cur
          den_b <- den_cur
          w_b   <- w
          
          for (j in take_vars) {
            wn <- max(w_min, w_b[j] * batch_factor)
            d  <- wn - w_b[j]
            if (abs(d) < 1e-12) next
            num_b   <- num_b + d * cache$N[[j]]
            den_b   <- den_b + d * cache$S[[j]]
            w_b[j]  <- wn
          }
          
          o_b  <- calc_obj_fast(num_b, den_b, n_sub, w_b)
          id_b <- o_b$id
          obj_b <- o_b$obj
          
          if (is.finite(obj_b) && obj_b < obj_cur - 1e-6) {
            num_cur <- num_b
            den_cur <- den_b
            w       <- w_b
            id      <- id_b
            obj_cur <- obj_b
            changed <- TRUE
            
            hist <- rbind(
              hist,
              data.frame(
                iter    = it,
                ID      = id,
                Obj     = obj_cur,
                changed = "BATCH",
                note    = paste(length(take_vars), "vars"),
                stringsAsFactors = FALSE
              )
            )
            
            if (verbose) {
              cat(sprintf(
                "   iter %d: [BATCH] x%.2f on %d vars (ID: %.3f, Obj: %.3f)\n",
                it, batch_factor, length(take_vars), id, obj_cur
              ))
            }
          }
        }
      }
    }
  }
  
  if (plot_progress && requireNamespace("ggplot2", quietly = TRUE)) {
    gp <- ggplot2::ggplot(hist, ggplot2::aes(iter, ID)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::labs(
        title = sprintf(
          "Stochastic Gower Optim (N=%d, Batch=%d)",
          n_sub, N_SAMPLE_PER_ITER
        )
      ) +
      ggplot2::theme_minimal()
    print(gp)
  }
  
  list(
    weights    = w,
    history    = hist,
    final_ID   = id,
    final_obj  = obj_cur,
    idx_used   = ix_sub_from_reps,
    subset_support = subset_support
  )
}

# --- Post-Optimization Selection ---

# knee_triangle <- function(w) {
#   y <- sort(as.numeric(w), decreasing = TRUE)
#   n <- length(y)
#   
#   if (n < 3L) {
#     return(list(k = n, thr = if (n) y[n] else NA_real_, curve = data.frame(i = seq_len(n), w = y, d = rep(0, n))))
#   }
#   
#   x <- seq_len(n)
#   num <- abs((y[n] - y[1]) * x - (n - 1) * y + n * y[1] - y[n])
#   den <- sqrt((y[n] - y[1])^2 + (n - 1)^2)
#   d <- num / den
#   k <- which.max(d)
#   
#   list(k = k, thr = y[k], curve = data.frame(i = x, w = y, d = d))
# }

knee_satopaa <- function(w) {
  y <- sort(as.numeric(w), decreasing = TRUE)
  n <- length(y)
  x <- seq_len(n)
  
  if (n < 3L) return(list(k = n, thr = if (n) y[n] else 0.0))
  
  # Normalize to Unit Square
  y_norm <- (y - min(y)) / (max(y) - min(y))
  x_norm <- (x - min(x)) / (max(x) - min(x))
  
  # Calculate Distance from Diagonal 
  # The sensitivity line is y = 1 - x (for decreasing curves).
  # Use abs() to detect 'concave' knees (Toes) that dip below the line.
  d <- abs(y_norm - (1 - x_norm))
  
  # Find the point of maximum deviation
  k <- which.max(d)
  
  list(k = k, thr = y[k])
}

# --- One-sided Kneedle for decreasing curves (Satopaa-style) ---
# Works on a decreasing curve after normalisation to [0,1].
# For decreasing curves, reference line is y = 1 - x.
# We take *positive* deviation only (knee above the line), not abs().
kneedle_decreasing <- function(y) {
  y <- as.numeric(y)
  n <- length(y)
  if (n < 3L) return(list(k = n, thr = if (n) y[n] else 0))
  
  x <- seq_len(n)
  
  # normalise to unit square
  yr <- range(y, na.rm = TRUE)
  if (!is.finite(yr[1]) || !is.finite(yr[2]) || yr[2] <= yr[1]) {
    return(list(k = n, thr = y[n]))
  }
  y_norm <- (y - yr[1]) / (yr[2] - yr[1])
  x_norm <- (x - 1) / (n - 1)
  
  # one-sided deviation from decreasing diagonal
  # knee candidates are where curve lies ABOVE the diagonal
  d <- y_norm - (1 - x_norm)
  d[!is.finite(d)] <- -Inf
  
  k <- which.max(d)
  list(k = k, thr = y[k], d = d, y_norm = y_norm, x_norm = x_norm)
}

knee_transition_empirical_floor <- function(w_sorted,
                                            plateau_mult = 0.98,
                                            tail_frac = 0.35,
                                            floor_slack = 1.05) {
  y <- as.numeric(w_sorted)
  n <- length(y)
  if (n < 3L) return(list(k = n, thr = y[n], plateau_end = 1L, floor_thr = y[n]))
  
  ymax <- y[1]
  
  # Plateau end: allow slight sag under the max (your picture needs this)
  plateau_end <- max(which(y >= plateau_mult * ymax), 1L)
  
  # Empirical floor from the tail: treat the last tail_frac as "candidate floor region"
  t0 <- max(plateau_end + 1L, floor((1 - tail_frac) * n))
  tail <- y[t0:n]
  
  # Robust floor level: median of tail (insensitive to a few slightly-higher tail points)
  floor_level <- stats::median(tail, na.rm = TRUE)
  
  # Define "in floor" as being close to that level (slack factor)
  floor_thr <- floor_level * floor_slack
  
  # First entry into floor after plateau; knee is just before that
  j0 <- which((y <= floor_thr) & (seq_len(n) > plateau_end))
  if (!length(j0)) {
    # If you never enter floor, pick steepest drop after plateau as fallback
    dy <- diff(y)
    idx <- plateau_end:(n - 1L)
    k_step <- idx[which.min(dy[idx])]
    k <- min(n, k_step + 1L)
    return(list(k = k, thr = y[k], plateau_end = plateau_end, floor_thr = floor_thr,
                floor_level = floor_level, j_floor = NA_integer_))
  }
  
  j_floor <- min(j0)
  k <- max(plateau_end + 1L, j_floor - 1L)
  
  list(k = k, thr = y[k], plateau_end = plateau_end, floor_thr = floor_thr,
       floor_level = floor_level, j_floor = j_floor)
}

# Knee = first rank where weight < 0.01 (after sorting descending).
# Returns cut_idx = knee-1 (i.e., keep everything before it).
# If nothing drops below 0.01, keeps all.
survivors_from_weights <- function(w,
                                   thr = 0.01,
                                   w_min = W_MIN,
                                   make_plot = TRUE,
                                   plot_file = "FIG_weight_curve_knee.png") {
  
  w_names <- names(w)
  if (is.null(w_names)) w_names <- paste0("V", seq_along(w))
  w <- pmax(w_min, as.numeric(w))
  names(w) <- w_names
  p <- length(w)
  
  w_sorted <- sort(w, decreasing = TRUE)
  
  knee_idx <- which(w_sorted < thr)[1]
  cut_idx <- if (is.na(knee_idx)) p else max(1L, knee_idx - 1L)
  
  thr_final <- w_sorted[cut_idx]
  survivors <- names(w_sorted)[seq_len(cut_idx)]
  active_eps <- GOWER_ACTIVE_EPS
  if (!is.finite(active_eps) || active_eps <= 0) active_eps <- 1e-8
  active_eps <- max(active_eps, .Machine$double.eps)
  active_survivors <- survivors[w_sorted[seq_len(cut_idx)] > (w_min + active_eps)]
  
  cat(sprintf("[Selection] knee(first<%.4f)=%s | Final: %d/%d (thr=%.6f)\n",
              thr,
              if (is.na(knee_idx)) "NA" else as.integer(knee_idx),
              length(survivors), p, thr_final))
  
  if (isTRUE(make_plot) && requireNamespace("ggplot2", quietly = TRUE)) {
    curve_df <- data.frame(rank = seq_along(w_sorted), weight = w_sorted)
    pt_final <- curve_df[cut_idx, , drop = FALSE]
    
    pplt <- ggplot2::ggplot(curve_df, ggplot2::aes(rank, weight)) +
      ggplot2::geom_line() +
      ggplot2::geom_hline(yintercept = thr, linetype = "dashed", linewidth = 0.35) +
      ggplot2::geom_point(data = pt_final, ggplot2::aes(color = "Selected Cut"),
                          shape = 4, size = 4, stroke = 2) +
      ggplot2::scale_color_manual(values = c("Selected Cut" = "blue")) +
      ggplot2::labs(x = "Rank (descending)", y = "Weight", color = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(legend.position = c(0.8, 0.8))
    
    print(pplt)
    if (exists("save_plot_gg")) save_plot_gg(plot_file, pplt, width = 6, height = 4, dpi = 300)
  }
  
  list(
    survivors = survivors,
    active_survivors = active_survivors,
    active_eps = active_eps,
    thr_tail = thr_final,
    w_clamped = w
  )
}

# --- PCA and Residuals Helpers ---

design_with_map <- function(X, nominal_missing = PCA_NOMINAL_MISSING) {
  Xg <- as.data.frame(X, check.names = TRUE, stringsAsFactors = FALSE)
  if (!ncol(Xg)) stop("[design_with_map] input has 0 columns")

  keep <- vapply(Xg, function(v) length(unique(na.omit(v))) >= 2L, logical(1))
  if (!any(keep)) stop("[design_with_map] all columns are NA-only or constant")
  Xg <- Xg[, keep, drop = FALSE]

  nominal_missing <- match.arg(tolower(nominal_missing), c("mode", "level"))
  imputation_rows <- list()
  blocks <- vector("list", ncol(Xg))
  block_varmap <- vector("list", ncol(Xg))
  block_levelmap <- vector("list", ncol(Xg))

  for (j in seq_along(Xg)) {
    nm <- names(Xg)[j]
    v <- Xg[[j]]

    if (is.ordered(v) && isTRUE(TREAT_ORDINALS_AS_NOMINAL)) {
      v <- factor(as.character(v), levels = levels(v), ordered = FALSE, exclude = NA)
    } else if (is.ordered(v)) {
      v <- as.numeric(v)
    } else if (is.logical(v)) {
      v <- factor(v, levels = c(FALSE, TRUE), ordered = FALSE, exclude = NA)
    } else if (!is.factor(v) && !is.numeric(v) && !is.integer(v)) {
      v <- factor(as.character(v), ordered = FALSE, exclude = NA)
    }

    if (is.factor(v)) {
      counts <- table(v, useNA = "no")
      min_level_count <- max(
        as.integer(GOWER_SUBSET_MIN_LEVEL_N),
        ceiling(RARE_LEVEL_MIN_PROP * nrow(Xg))
      )
      rare_levels <- names(counts)[counts < min_level_count]
      if (length(rare_levels)) {
        values <- as.character(v)
        values[!is.na(values) & values %in% rare_levels] <- NA_character_
        v <- factor(
          values,
          levels = setdiff(levels(v), rare_levels),
          ordered = FALSE,
          exclude = NA
        )
      }
      n_missing <- sum(is.na(v))
      if (n_missing && nominal_missing == "mode") {
        counts <- table(v, useNA = "no")
        if (!length(counts)) stop("[design_with_map] no observed level for ", nm)
        fill <- names(counts)[which.max(counts)]
        values <- as.character(v)
        values[is.na(values)] <- fill
        v <- factor(values, levels = levels(v), ordered = FALSE, exclude = NA)
      } else if (n_missing && nominal_missing == "level") {
        values <- as.character(v)
        values[is.na(values)] <- "(Missing)"
        v <- factor(values, levels = c(levels(v), "(Missing)"), ordered = FALSE)
        fill <- "(Missing)"
      } else {
        fill <- NA_character_
      }

      lev <- levels(droplevels(v))
      if (length(lev) < 2L) next
      block <- vapply(lev, function(level) as.numeric(v == level), numeric(nrow(Xg)))
      if (is.null(dim(block))) block <- matrix(block, ncol = 1L)
      colnames(block) <- paste0(nm, "__level__", make.names(lev, unique = TRUE))
      blocks[[j]] <- block
      block_varmap[[j]] <- rep(nm, length(lev))
      block_levelmap[[j]] <- lev
      if (n_missing) {
        imputation_rows[[length(imputation_rows) + 1L]] <- data.frame(
          source_var = nm,
          n_missing = n_missing,
          method = nominal_missing,
          fill_value = fill,
          stringsAsFactors = FALSE
        )
      }
    } else {
      v <- as.numeric(v)
      n_missing <- sum(is.na(v))
      fill <- NA_real_
      if (n_missing) {
        fill <- stats::median(v, na.rm = TRUE)
        if (!is.finite(fill)) stop("[design_with_map] cannot impute numeric column ", nm)
        v[is.na(v)] <- fill
        imputation_rows[[length(imputation_rows) + 1L]] <- data.frame(
          source_var = nm,
          n_missing = n_missing,
          method = "median",
          fill_value = as.character(fill),
          stringsAsFactors = FALSE
        )
      }
      blocks[[j]] <- matrix(v, ncol = 1L, dimnames = list(NULL, nm))
      block_varmap[[j]] <- nm
      block_levelmap[[j]] <- NA_character_
    }
  }

  nonempty <- lengths(blocks) > 0L
  if (!any(nonempty)) stop("[design_with_map] no encodable columns remain")
  MM <- do.call(cbind, blocks[nonempty])
  storage.mode(MM) <- "double"
  rownames(MM) <- rownames(Xg)
  varmap <- unlist(block_varmap[nonempty], use.names = FALSE)
  levelmap <- unlist(block_levelmap[nonempty], use.names = FALSE)

  ok <- apply(MM, 2L, function(col) {
    vv <- stats::var(as.numeric(col), na.rm = TRUE)
    is.finite(vv) && vv > 1e-12
  })
  if (!any(ok)) stop("[design_with_map] all encoded columns were ~zero-variance")

  MM <- MM[, ok, drop = FALSE]
  attr(MM, "varmap") <- varmap[ok]
  attr(MM, "levelmap") <- levelmap[ok]
  attr(MM, "imputation") <- if (length(imputation_rows)) {
    do.call(rbind, imputation_rows)
  } else {
    data.frame(
      source_var = character(), n_missing = integer(), method = character(),
      fill_value = character(), stringsAsFactors = FALSE
    )
  }
  MM
}

residualise_foldsafe <- function(Xenc, Base, folds, k_gam = 6) {
  n <- nrow(Base)
  V <- colnames(Xenc)
  E <- matrix(NA_real_, n, length(V), dimnames = list(rownames(Base), V))
  sm_terms <- paste0("s(b", seq_len(ncol(Base)), ",k=", k_gam, ")")
  
  for (v in V) {
    z <- as.numeric(Xenc[, v])
    for (k in sort(unique(folds))) {
      tr <- which(folds != k)
      te <- which(folds == k)
      dftr <- data.frame(v = z[tr], Base[tr, , drop = FALSE])
      fml <- reformulate(sm_terms, response = "v")
      
      g <- try(mgcv::gam(fml, data = dftr, method = "REML"), silent = TRUE)
      if (inherits(g, "try-error")) next
      
      mu <- as.numeric(predict(g, newdata = data.frame(Base[te, , drop = FALSE]), type = "response"))
      E[te, v] <- z[te] - mu
    }
  }
  E[, colSums(is.finite(E)) > 0, drop = FALSE]
}

# --- Diagnostic Metrics Helpers ---

score_item_base <- function(nm, Z, varmap) {
  idx <- which(varmap == nm)
  if (!length(idx)) return(rep(NA_real_, nrow(Z)))
  
  if (length(idx) == 1L) return(as.numeric(Z[, idx]))
  
  # For multi-column items (factors), use first PC
  sc <- try(suppressWarnings(prcomp(Z[, idx, drop = FALSE], rank. = 1)$x[, 1]), silent = TRUE)
  if (inherits(sc, "try-error")) return(rep(NA_real_, nrow(Z)))
  
  as.numeric(sc)
}

e_from_E <- function(nm, E_scaled, Z, varmap) {
  idx_all <- which(varmap == nm)
  if (!length(idx_all)) return(rep(NA_real_, nrow(Z)))
  
  enc_names <- colnames(Z)
  cols_item <- enc_names[idx_all]
  
  cols_in_E <- intersect(cols_item, colnames(E_scaled))
  if (!length(cols_in_E)) return(rep(NA_real_, nrow(Z)))
  
  if (length(cols_in_E) == 1L) {
    return(as.numeric(E_scaled[, cols_in_E]))
  }
  
  # Use only the columns that exist in E_scaled for the item projection
  Zsub <- Z[, cols_in_E, drop = FALSE]
  Esub <- E_scaled[, cols_in_E, drop = FALSE]
  
  pc1 <- try(prcomp(Zsub, rank. = 1), silent = TRUE)
  if (inherits(pc1, "try-error")) return(rep(NA_real_, nrow(Z)))
  
  as.numeric(as.matrix(Esub) %*% pc1$rotation[, 1])
}

r2_base_linear <- function(Base, v, d = BASE_DIM) {
  if (nrow(Base) != length(v)) return(NA_real_)
  
  d_eff <- min(as.integer(d), ncol(Base))
  if (d_eff < 1L) return(NA_real_)
  
  Bsub <- as.data.frame(Base[, seq_len(d_eff), drop = FALSE])
  colnames(Bsub) <- paste0("b", seq_len(d_eff))
  
  df <- cbind.data.frame(v = as.numeric(v), Bsub)
  fit <- try(stats::lm(v ~ ., data = df), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  summary(fit)$r.squared
}

r2_residual_cv <- function(e, nb, folds = CV_FOLDS, seed = SEED_GLOBAL, key = NULL) {
  n <- length(e)
  if (n < 10 || var(e, na.rm = TRUE) == 0 || all(is.na(e))) return(NA_real_)
  
  # Accept either a scalar K or a precomputed fold_id vector
  if (length(folds) > 1L) {
    fold_id <- as.integer(folds)
    if (length(fold_id) != n) stop("fold_id length mismatch: ", length(fold_id), " vs n=", n)
    fold_id[!is.finite(fold_id)] <- 1L
  } else {
    K <- as.integer(folds[1])
    if (!is.finite(K) || K < 2L) K <- 2L
    if (K > n) K <- n
    
    s <- if (is.null(key)) as.integer(seed) else .seed_from_key(seed, paste0("r2_resid|", key))
    fold_id <- .with_seed(s, sample(rep_len(seq_len(K), n)))
  }
  
  pred <- rep(NA_real_, n)
  
  for (f in sort(unique(fold_id))) {
    te <- which(fold_id == f)
    if (!length(te)) next
    
    e_mask <- e
    e_mask[te] <- NA_real_
    
    nb_vals <- matrix(e_mask[nb], nrow = n)
    pr <- rowMeans(nb_vals, na.rm = TRUE)
    
    pr[is.na(pr)] <- mean(e[-te], na.rm = TRUE)
    pred[te] <- pr[te]
  }
  
  ve <- stats::var(e, na.rm = TRUE)
  mse <- mean((e - pred)^2, na.rm = TRUE)
  
  if (!is.finite(ve) || ve <= 1e-9) return(0)
  max(0, 1 - mse / ve)
}

# --- Parallel Worker for Interactions ---

interact_worker <- function(i, grid, Zmat, BaseDF, DxDF, map) {
  
  # Prevent over-subscription: Lock internal threads to 1 per worker.
  # This avoids the "exploding threads" issue when parallelizing high-level tasks.
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    try(RhpcBLASctl::blas_set_num_threads(1), silent = TRUE)
    try(RhpcBLASctl::omp_set_num_threads(1), silent = TRUE)
  }
  Sys.setenv(OMP_NUM_THREADS = "1")
  Sys.setenv(MKL_NUM_THREADS = "1")
  Sys.setenv(OPENBLAS_NUM_THREADS = "1")
  options(mc.cores = 1)
  
  dx_name <- grid$dx[i]
  nm      <- grid$var[i]
  
  y <- as.integer(DxDF[[dx_name]] > 0)
  v <- score_item_base(nm, Zmat, map) 
  
  if (all(is.na(v)) || stats::sd(v, na.rm = TRUE) < 1e-6) return(NULL)
  
  # Run OOF GAMs
  res <- oof_R2_two_gams(v, BaseDF, y, K_target = 5, k_gam = 10,
                         seed = SEED_GLOBAL,
                         seed_key = paste0(dx_name, "::", nm))
  
  if (is.finite(unname(res["dR2"]))) {
    return(data.frame(
      dx = dx_name, var = nm,
      R2_add = unname(res["R2_add"]), R2_int = unname(res["R2_int"]),
      dR2 = unname(res["dR2"]), p_like = unname(res["p_like"]),
      stringsAsFactors = FALSE
    ))
  }
  return(NULL)
}

# ==============================================================================
# 2. Data Ingest
# ==============================================================================

read_input_table <- function(path, ...) {
  first_line <- readLines(path, n = 1L, warn = FALSE)
  count_fixed <- function(x, pattern) {
    m <- gregexpr(pattern, x, fixed = TRUE)[[1L]]
    if (identical(m, -1L)) 0L else length(m)
  }
  
  n_semicolon <- count_fixed(first_line, ";")
  n_comma <- count_fixed(first_line, ",")
  delim <- if (n_semicolon > n_comma) ";" else ","
  decimal_mark <- if (identical(delim, ";")) "," else "."
  
  out <- readr::read_delim(
    path,
    delim = delim,
    locale = readr::locale(decimal_mark = decimal_mark),
    progress = FALSE,
    show_col_types = FALSE,
    ...
  )
  attr(out, "input_delim") <- delim
  attr(out, "input_decimal_mark") <- decimal_mark
  out
}

df <- read_input_table(PSY_CSV)
cat(sprintf(
  "[ingest] Read %s as delim='%s', decimal='%s': %d rows x %d columns\n",
  PSY_CSV, attr(df, "input_delim"), attr(df, "input_decimal_mark"),
  nrow(df), ncol(df)
))

id_col <- if ("participant_id" %in% names(df)) "participant_id" else names(df)[1]
ids_all <- as.character(df[[id_col]])

# Clean columns
drop_pattern <- "(?i)diagnosis|SCID|NODIAG"
cols_to_drop <- grep(drop_pattern, names(df), value = TRUE)

X <- dplyr::select(df, -dplyr::all_of(c(id_col, cols_to_drop))) |>
  as.data.frame(stringsAsFactors = FALSE)

# Apply explicit item metadata before heuristic type coercion. The psychometric
# CSV stores response codes as numbers, so R cannot otherwise distinguish an
# ordinal response scale from a genuinely continuous variable.
resolve_schema_path <- function(path, psych_csv) {
  if (is.null(path) || !nzchar(path)) return("")
  if (file.exists(path)) return(normalizePath(path, mustWork = TRUE))
  candidate <- file.path(dirname(psych_csv), basename(path))
  if (file.exists(candidate)) return(normalizePath(candidate, mustWork = TRUE))
  path
}

parse_schema_levels <- function(text) {
  if (is.na(text) || !nzchar(trimws(text))) return(numeric())
  values <- suppressWarnings(as.numeric(strsplit(text, "|", fixed = TRUE)[[1L]]))
  values[is.finite(values)]
}

schema_path <- resolve_schema_path(ITEM_RESPONSE_SCHEMA_CSV, PSY_CSV)
response_code_audit <- NULL
if (nzchar(schema_path) && file.exists(schema_path)) {
  response_schema <- read_input_table(schema_path, col_types = readr::cols(.default = readr::col_character()))
  required_schema_cols <- c("var", "kind", "valid_levels")
  missing_schema_cols <- setdiff(required_schema_cols, names(response_schema))
  if (length(missing_schema_cols)) {
    stop("[schema] Missing required columns: ", paste(missing_schema_cols, collapse = ", "))
  }
  if (anyDuplicated(response_schema$var)) stop("[schema] Duplicate variable names")
  missing_schema_vars <- setdiff(names(X), response_schema$var)
  if (length(missing_schema_vars) && isTRUE(REQUIRE_ITEM_RESPONSE_SCHEMA)) {
    stop("[schema] Variables absent from item_response_schema.csv: ",
         paste(head(missing_schema_vars, 20L), collapse = ", "))
  }

  audit_rows <- list()
  for (nm in intersect(names(X), response_schema$var)) {
    row <- response_schema[match(nm, response_schema$var), , drop = FALSE]
    kind <- tolower(trimws(row$kind[[1L]]))
    if (identical(kind, "nominal")) {
      valid <- parse_schema_levels(row$valid_levels[[1L]])
      if (length(valid) < 2L) stop("[schema] Fewer than two valid levels for ", nm)
      values <- suppressWarnings(as.numeric(X[[nm]]))
      invalid <- is.finite(values) & !values %in% valid
      invalid_values <- sort(unique(values[invalid]))
      if (any(invalid) && !identical(tolower(INVALID_RESPONSE_ACTION), "missing")) {
        stop("[schema] Invalid response code(s) in ", nm, ": ",
             paste(sort(unique(values[invalid])), collapse = ", "))
      }
      values[invalid] <- NA_real_
      X[[nm]] <- factor(values, levels = valid, ordered = FALSE, exclude = NA)
      audit_rows[[length(audit_rows) + 1L]] <- data.frame(
        var = nm,
        kind = kind,
        n_invalid_to_missing = sum(invalid),
        invalid_codes = paste(invalid_values, collapse = "|"),
        stringsAsFactors = FALSE
      )
    } else if (identical(kind, "continuous")) {
      X[[nm]] <- suppressWarnings(as.numeric(X[[nm]]))
      audit_rows[[length(audit_rows) + 1L]] <- data.frame(
        var = nm, kind = kind, n_invalid_to_missing = 0L,
        invalid_codes = "", stringsAsFactors = FALSE
      )
    } else {
      stop("[schema] Unknown kind for ", nm, ": ", kind)
    }
  }
  response_code_audit <- dplyr::bind_rows(audit_rows)
  cat(sprintf(
    "[schema] Applied %d item definitions (%d nominal, %d continuous); %d invalid code(s) set to missing.\n",
    nrow(response_code_audit),
    sum(response_code_audit$kind == "nominal"),
    sum(response_code_audit$kind == "continuous"),
    sum(response_code_audit$n_invalid_to_missing)
  ))
} else if (isTRUE(REQUIRE_ITEM_RESPONSE_SCHEMA)) {
  stop("[schema] Required item response schema not found: ", schema_path)
}

# Heuristic fallback for matrices without an explicit schema.
is_small_int_scale <- function(v) {
  vn <- suppressWarnings(as.numeric(v))
  if (all(is.na(vn))) return(FALSE)
  u <- sort(unique(na.omit(vn)))
  k <- length(u)
  k >= 3 && k <= as.integer(SMALL_INT_NOMINAL_MAX_UNIQUE) &&
    all(abs(u - round(u)) < 1e-8)
}
parse_numeric_text <- function(v) {
  ch <- trimws(as.character(v))
  is_missing_text <- is.na(ch) | ch %in% c("", ".", "NA", "NaN")
  n_observed <- sum(!is_missing_text)
  if (!n_observed) return(NULL)
  
  parse_with <- function(decimal_mark) {
    suppressWarnings(readr::parse_double(
      ch,
      locale = readr::locale(decimal_mark = decimal_mark),
      na = c("", ".", "NA", "NaN")
    ))
  }
  comma <- parse_with(",")
  dot <- parse_with(".")
  comma_ok <- sum(!is.na(comma[!is_missing_text]))
  dot_ok <- sum(!is.na(dot[!is_missing_text]))
  best <- if (comma_ok >= dot_ok) comma else dot
  best_ok <- max(comma_ok, dot_ok)
  
  if (best_ok == n_observed) best else NULL
}
coerce_small_int_numeric <- isTRUE(COERCE_SMALL_INT_NUMERIC_TO_FACTOR)
coerce_numeric_01 <- isTRUE(COERCE_NUMERIC_01_TO_FACTOR)
numeric_01_predictors <- character(0)

for (nm in names(X)) {
  v <- X[[nm]]
  if (is.factor(v)) next
  if (is.character(v)) {
    vn <- parse_numeric_text(v)
    if (!is.null(vn)) v <- vn
  }
  if (all(v %in% c(0, 1, NA))) {
    numeric_01_predictors <- union(numeric_01_predictors, nm)
    if (coerce_numeric_01 && isTRUE(TREAT_ORDINALS_AS_NOMINAL)) {
      X[[nm]] <- factor(v, ordered = FALSE, exclude = if (isTRUE(MISSING_AS_NOMINAL_LEVEL)) NULL else NA)
    } else if (coerce_numeric_01) {
      X[[nm]] <- factor(v, levels = c(0, 1), ordered = TRUE)
    } else {
      X[[nm]] <- as.numeric(v)
    }
  } else if (coerce_small_int_numeric && is_small_int_scale(v)) {
    if (isTRUE(TREAT_ORDINALS_AS_NOMINAL)) {
      X[[nm]] <- factor(
        as.integer(round(as.numeric(v))),
        ordered = FALSE,
        exclude = if (isTRUE(MISSING_AS_NOMINAL_LEVEL)) NULL else NA
      )
    } else {
      X[[nm]] <- factor(as.integer(round(as.numeric(v))), ordered = TRUE)
    }
  } else if (is.numeric(v)) {
    X[[nm]] <- as.numeric(v)
  } else {
    X[[nm]] <- factor(v)
  }
}

if (length(numeric_01_predictors)) {
  write_csv(
    data.frame(
      var = numeric_01_predictors,
      coerced_to_factor = isTRUE(coerce_numeric_01),
      stringsAsFactors = FALSE
    ),
    "numeric_01_predictors.csv"
  )
}

low_support_numeric_report <- function(X,
                                       max_unique = LOW_SUPPORT_NUMERIC_MAX_UNIQUE,
                                       min_prop = LOW_SUPPORT_NUMERIC_MIN_PROP,
                                       min_n = LOW_SUPPORT_NUMERIC_MIN_N) {
  n <- nrow(X)
  if (!n || !ncol(X)) return(data.frame())
  
  min_count <- max(as.integer(min_n), ceiling(min_prop * n))
  rows <- vector("list", ncol(X))
  
  for (i in seq_along(X)) {
    nm <- names(X)[[i]]
    v <- X[[i]]
    if (!is.numeric(v)) next
    
    vv <- v[is.finite(v)]
    u <- sort(unique(vv))
    if (length(u) < 2L || length(u) > max_unique) next
    
    tb <- table(vv)
    rare_count <- min(as.integer(tb))
    if (rare_count < min_count) {
      rare_values <- names(tb)[as.integer(tb) == rare_count]
      rows[[i]] <- data.frame(
        var = nm,
        n_unique = length(u),
        rare_count = rare_count,
        rare_prop = rare_count / n,
        min_count_required = min_count,
        rare_values = paste(head(rare_values, 12L), collapse = "|"),
        stringsAsFactors = FALSE
      )
    }
  }
  
  out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
  if (is.null(out)) data.frame() else out[order(out$rare_count, out$var), , drop = FALSE]
}

low_support_numeric <- low_support_numeric_report(X)
if (nrow(low_support_numeric)) {
  if (isTRUE(REPORT_LOW_SUPPORT_NUMERIC)) {
    write_csv(low_support_numeric, "low_support_numeric_predictors.csv")
  }
  if (isTRUE(DROP_LOW_SUPPORT_NUMERIC)) {
    X <- X[, setdiff(names(X), low_support_numeric$var), drop = FALSE]
    cat(sprintf(
      "[ingest] Dropped %d low-support numeric predictors; see low_support_numeric_predictors.csv\n",
      nrow(low_support_numeric)
    ))
  } else {
    cat(sprintf(
      "[ingest] Flagged %d low-support numeric predictors; DROP_LOW_SUPPORT_NUMERIC=FALSE, so none were removed.\n",
      nrow(low_support_numeric)
    ))
  }
}

predictor_missingness <- data.frame(
  var = names(X),
  n_missing = as.integer(colSums(is.na(X))),
  n_observed = as.integer(colSums(!is.na(X))),
  prop_missing = as.numeric(colMeans(is.na(X))),
  stringsAsFactors = FALSE
)
predictor_missingness <- predictor_missingness[order(
  -predictor_missingness$prop_missing,
  predictor_missingness$var
), , drop = FALSE]
if (isTRUE(REPORT_PREDICTOR_MISSINGNESS)) {
  write_csv(predictor_missingness, "predictor_missingness.csv")
}
if (!is.null(response_code_audit)) {
  write_csv(response_code_audit, "response_code_audit.csv")
}

max_missing_prop <- MAX_PREDICTOR_MISSING_PROP
high_missing_predictors <- predictor_missingness$var[
  predictor_missingness$prop_missing > max_missing_prop
]
if (length(high_missing_predictors) && isTRUE(DROP_HIGH_MISSING_PREDICTORS)) {
  X <- X[, setdiff(names(X), high_missing_predictors), drop = FALSE]
  cat(sprintf(
    "[ingest] Dropped %d predictors with > %.1f%% missing values; see predictor_missingness.csv\n",
    length(high_missing_predictors),
    100 * max_missing_prop
  ))
} else if (length(high_missing_predictors)) {
  cat(sprintf(
    "[ingest] Flagged %d predictors with > %.1f%% missing values; DROP_HIGH_MISSING_PREDICTORS=FALSE, so none were removed.\n",
    length(high_missing_predictors),
    100 * max_missing_prop
  ))
}

cat(sprintf("[debug] Antes de complete.cases: nrow(X)=%d, ncol(X)=%d\n", nrow(X), ncol(X)))
na_by_col <- colSums(is.na(X))
cat("[debug] Top 20 colunas com mais NA:\n")
print(sort(na_by_col, decreasing = TRUE)[1:20])

prop_obs_row <- rowMeans(!is.na(X))
cat("[debug] resumo proporção de itens observados por linha:\n")
print(summary(prop_obs_row))

row_missingness <- data.frame(
  participant_id = ids_all,
  n_missing = rowSums(is.na(X)),
  prop_missing = rowMeans(is.na(X)),
  stringsAsFactors = FALSE
)
write_csv(row_missingness, "participant_predictor_missingness.csv")

keep <- row_missingness$prop_missing <= MAX_PARTICIPANT_MISSING_PROP &
  !is.na(ids_all) & ids_all != ""
X <- X[keep, , drop = FALSE]
ids_all <- ids_all[keep]

if (nrow(X) == 0L) {
  stop(sprintf(
    "[ingest] No rows remained after applying MAX_PARTICIPANT_MISSING_PROP=%.3f (PSY_CSV=%s).",
    MAX_PARTICIPANT_MISSING_PROP,
    PSY_CSV
  ))
}
cat(sprintf(
  "[ingest] Retained %d/%d participants with predictor missingness <= %.1f%%; Gower uses pairwise available variables and PCA applies the configured deterministic imputation.\n",
  nrow(X), length(keep), 100 * MAX_PARTICIPANT_MISSING_PROP
))

rownames(X) <- make.unique(ids_all)
colnames(X) <- make.names(colnames(X), unique = TRUE)

if (isTRUE(VALIDATE_ENCODING_ONLY)) {
  nominal_names <- names(X)[vapply(X, function(v) is.factor(v) && !is.ordered(v), logical(1))]
  PX_check <- prep_X_for_gower(
    X, rare_prop = RARE_LEVEL_MIN_PROP, do_jitter = FALSE, seed = SEED_JITTER
  )
  incorrectly_binned <- intersect(nominal_names, PX_check$binned_continuous_cols)
  if (length(incorrectly_binned)) {
    stop("[validate] Nominal variables entered the continuous-binning branch: ",
         paste(head(incorrectly_binned, 20L), collapse = ", "))
  }
  if (any(vapply(PX_check$X[nominal_names], is.ordered, logical(1)))) {
    stop("[validate] At least one nominal variable remained ordered in the Gower input")
  }

  Xenc_check <- design_with_map(X, nominal_missing = PCA_NOMINAL_MISSING)
  varmap_check <- attr(Xenc_check, "varmap")
  levelmap_check <- attr(Xenc_check, "levelmap")
  encoded_counts_check <- table(varmap_check)
  if (any(encoded_counts_check[nominal_names] < 2L, na.rm = TRUE)) {
    stop("[validate] At least one nominal variable lacks response-level dummy columns")
  }
  if (any(is.na(levelmap_check[varmap_check %in% nominal_names]))) {
    stop("[validate] Missing response-level labels in nominal encoding")
  }

  validation_summary <- data.frame(
    n_participants = nrow(X),
    n_source_variables = ncol(X),
    n_nominal_variables = length(nominal_names),
    n_continuous_variables = ncol(X) - length(nominal_names),
    n_encoded_columns = ncol(Xenc_check),
    n_gower_binned_continuous = length(PX_check$binned_continuous_cols),
    n_invalid_codes_to_missing = if (is.null(response_code_audit)) 0L else
      sum(response_code_audit$n_invalid_to_missing),
    stringsAsFactors = FALSE
  )
  write_csv(validation_summary, "encoding_validation_summary.csv")
  cat("[validate] Nominal preprocessing and response-level encoding passed.\n")
  signalCondition(structure(
    list(message = "encoding validation complete"),
    class = c("encoding_validation_complete", "condition")
  ))
  stop("[validate] encoding_validation_complete condition was not handled")
}

degenerate_value <- constant_profile_value(X)
degenerate_mask <- !is.na(degenerate_value)
if (any(degenerate_mask)) {
  degenerate_constant_profiles <- data.frame(
    row = which(degenerate_mask),
    participant_id = rownames(X)[degenerate_mask],
    constant_value = degenerate_value[degenerate_mask],
    stringsAsFactors = FALSE
  )
  
  if (isTRUE(WRITE_DEGENERATE_CSV)) {
    write_csv(degenerate_constant_profiles, "degenerate_constant_profiles.csv")
  }
  
  deg_counts <- sort(table(degenerate_constant_profiles$constant_value), decreasing = TRUE)
  top_counts <- head(deg_counts, 5L)
  warning(sprintf(
    "[data] Found %d fully-constant profiles; kept in analysis%s. Top values: %s",
    nrow(degenerate_constant_profiles),
    if (isTRUE(WRITE_DEGENERATE_CSV)) " (see degenerate_constant_profiles.csv)" else "",
    paste(sprintf("%s=%d", names(top_counts), as.integer(top_counts)), collapse = ", ")
  ))
} else {
  degenerate_constant_profiles <- NULL
}

lock_vars_intersect <- intersect(LOCK_VARS, colnames(X))
if (length(lock_vars_intersect)) {
  message("[weights] lock_vars: ",
          paste(lock_vars_intersect, collapse = ", "))
}

# --- Diagnoses Table Ingest (OPTIONAL) ---
DX_AVAILABLE <- FALSE
diag_wide_full <- NULL

safe_ingest_dx <- function(path, ids_here, optional = TRUE) {
  if (!nzchar(path) || !file.exists(path)) {
    if (optional) {
      message("[dx] Diagnoses file missing; proceeding without diagnoses.")
      return(list(
        available = FALSE,
        diag_wide_full = data.frame(participant_id = ids_here, stringsAsFactors = FALSE)
      ))
    } else {
      stop("Diagnoses file not found: ", path)
    }
  }
  
  # Try to read and align
  dx <- try(read_input_table(path, col_types = readr::cols()), silent = TRUE)
  if (inherits(dx, "try-error")) {
    if (optional) {
      warning("[dx] Read failed; disabling diagnoses (optional).")
      return(list(
        available = FALSE,
        diag_wide_full = data.frame(participant_id = ids_here, stringsAsFactors = FALSE)
      ))
    } else {
      stop(dx)
    }
  }
  
  dx <- dx %>% dplyr::mutate(participant_id = as.character(participant_id))
  
  mm <- match(ids_here, as.character(dx$participant_id))
  # If the join fails for all rows, gracefully disable
  if (all(is.na(mm))) {
    if (optional) {
      warning("[dx] Join failed for 100% of rows; disabling diagnoses (optional).")
      return(list(
        available = FALSE,
        diag_wide_full = data.frame(participant_id = ids_here, stringsAsFactors = FALSE)
      ))
    } else {
      stop(sprintf("Diagnoses join failed for %d/%d rows.", sum(is.na(mm)), length(mm)))
    }
  }
  
  # Reorder and fill missing with 0
  dxA <- dx[mm, , drop = FALSE]
  dxA$participant_id <- ids_here
  if (ncol(dxA) > 1) {
    dxA <- dxA %>%
      dplyr::mutate(dplyr::across(-participant_id, ~ {
        x <- suppressWarnings(as.integer(.))
        x[is.na(x)] <- 0L
        pmin(pmax(x, 0L), 1L)
      }))
  }
  
  avail <- (ncol(dxA) > 1) && any(colSums(dxA[, -1, drop = FALSE] > 0, na.rm = TRUE) > 0)
  list(available = avail, diag_wide_full = dxA)
}

ids_here <- rownames(X)
dx_ing   <- safe_ingest_dx(DIAG_CSV, ids_here, optional = isTRUE(DX_OPTIONAL))

DX_AVAILABLE   <- isTRUE(dx_ing$available)
diag_wide_full <- dx_ing$diag_wide_full

cat(sprintf("[ingest] X rows=%d, cols=%d | DX_available=%s | DX_cols=%d\n",
            nrow(X), ncol(X), DX_AVAILABLE, max(1L, ncol(diag_wide_full))))

# Align to DX_wide matrix (or empty)
if (DX_AVAILABLE && ncol(diag_wide_full) > 1) {
  DX_wide <- as.data.frame(diag_wide_full[, -1, drop = FALSE])
  rownames(DX_wide) <- diag_wide_full$participant_id
} else {
  DX_wide <- NULL
}

# ==============================================================================
# 3. Deduplication and Coreset Selection
# ==============================================================================

dedup_mode_eff <- if (!isTRUE(DO_DEDUP)) "none" else DEDUP_MODE
if (!dedup_mode_eff %in% c("gower_complete", "hash_exact", "hash_round", "none")) {
  stop("Unknown DEDUP_MODE: ", dedup_mode_eff)
}

if (identical(dedup_mode_eff, "gower_complete")) {
  PX <- prep_X_for_gower(X, rare_prop = RARE_LEVEL_MIN_PROP, do_jitter = TRUE, seed = SEED_JITTER)
  X_for_id <- PX$X
  Dg <- gower_dist(X_for_id, type_list = PX$type, weights = PX$weights)
  
  # Auto-detect EPS_DEDUP using knee method if NA
  if (is.na(EPS_DEDUP)) {
    d1 <- first_nn_d1(Dg)
    qlo <- as.numeric(stats::quantile(d1[is.finite(d1)], probs = c(0.001, 0.01), na.rm = TRUE))
    eps_grid <- seq(from = max(0, min(qlo, na.rm = TRUE) - 0.02),
                    to = min(0.50, stats::quantile(d1, 0.30, na.rm = TRUE)),
                    by = 0.005)
    cc <- collapse_curve(Dg, eps_grid)
    dprop <- diff(cc$prop_retained) / diff(cc$eps)
    knee_i <- which.min(dprop)
    EPS_DEDUP <- mean(cc$eps[c(knee_i, knee_i + 1)])
  }
  
  gr_all <- complete_groups(Dg, EPS_DEDUP)
  
  # Identify medoids of complete-linkage groups
  Dm_g <- as.matrix(Dg)
  diag(Dm_g) <- 0
  split_idx <- split(seq_len(nrow(Dm_g)), gr_all)
  reps <- vapply(split_idx, function(ix) {
    ix[which.min(rowSums(Dm_g[ix, ix, drop = FALSE]))]
  }, integer(1))
  mult <- as.integer(lengths(split_idx))
  
  Dg_rep <- stats::as.dist(Dm_g[reps, reps, drop = FALSE])
  core_idx_rep <- twonn_core_by_slope(Dg_rep)
  
  cat(sprintf("[Dedup:gower_complete] eps=%.3f | reps=%d of %d | core_rep=%d\n",
              EPS_DEDUP, length(reps), nrow(X), length(core_idx_rep)))
  
  if (isTRUE(WRITE_DEDUP_CSV)) {
    mult_df <- data.frame(
      rep_row = reps,
      representative_id = rownames(X_for_id)[reps],
      multiplicity = mult
    )
    write_csv(mult_df, sprintf("near_duplicate_groups_complete_eps%g.csv", EPS_DEDUP))
  }
  
} else if (dedup_mode_eff %in% c("hash_exact", "hash_round")) {
  # Hash dedup disables jitter so duplicate keys stay deterministic.
  PX <- prep_X_for_gower(X, rare_prop = RARE_LEVEL_MIN_PROP, do_jitter = FALSE, seed = SEED_JITTER)
  X_for_id <- PX$X
  
  split_idx <- dedup_groups_from_hash(
    X_for_id,
    mode = dedup_mode_eff,
    digits = DEDUP_HASH_DIGITS,
    na_token = DEDUP_HASH_NA
  )
  
  reps <- vapply(split_idx, `[`, integer(1), 1L)
  mult <- as.integer(lengths(split_idx))
  core_idx_rep <- seq_along(reps)
  
  cat(sprintf("[Dedup:%s] reps=%d of %d | core_rep=%d | hash_digits=%d\n",
              dedup_mode_eff, length(reps), nrow(X), length(core_idx_rep), DEDUP_HASH_DIGITS))
  
  if (isTRUE(WRITE_DEDUP_CSV)) {
    mult_df <- data.frame(
      rep_row = reps,
      representative_id = rownames(X_for_id)[reps],
      multiplicity = mult
    )
    suffix <- if (identical(dedup_mode_eff, "hash_round")) {
      paste0("hash_round_d", DEDUP_HASH_DIGITS)
    } else {
      "hash_exact"
    }
    write_csv(mult_df, sprintf("near_duplicate_groups_%s.csv", suffix))
  }
  
} else {
  reps <- seq_len(nrow(X))
  mult <- rep(1L, nrow(X))
  core_idx_rep <- seq_along(reps)
  
  cat(sprintf("[Dedup:none] reps=%d of %d | core_rep=%d\n",
              length(reps), nrow(X), length(core_idx_rep)))
  
  if (isTRUE(WRITE_DEDUP_CSV)) {
    mult_df <- data.frame(
      rep_row = reps,
      representative_id = rownames(X)[reps],
      multiplicity = mult
    )
    write_csv(mult_df, "near_duplicate_groups_none.csv")
  }
}

# ==============================================================================
# 4. Gower Weight Optimization
# ==============================================================================

# Drop constant columns and optionally trim high-correlation pairs
drop_constant_cols <- function(X) {
  keep <- vapply(X, function(v) length(unique(na.omit(v))) >= 2L, logical(1))
  X[, keep, drop = FALSE]
}
X_pred <- drop_constant_cols(X)

if (DO_CORR_TRIM) {
  corr_trim <- function(X, thr = 0.95, lock_vars = character(0)) {
    num <- X[, vapply(X, is.numeric, logical(1)), drop = FALSE]
    if (!ncol(num)) return(X)
    
    C <- suppressWarnings(cor(num, use = "pairwise.complete.obs"))
    drop <- character(0)
    
    for (i in seq_len(ncol(C) - 1)) {
      vi <- colnames(C)[i]
      if (vi %in% drop) next
      
      j_cand <- which(abs(C[i, (i + 1):ncol(C)]) >= thr) + i
      if (!length(j_cand)) next
      
      for (j in j_cand) {
        vj <- colnames(C)[j]
        if (vj %in% drop) next
        
        if (vi %in% lock_vars && !(vj %in% lock_vars)) {
          # keep vi, drop vj
          drop <- union(drop, vj)
        } else if (vj %in% lock_vars && !(vi %in% lock_vars)) {
          # keep vj, drop vi
          drop <- union(drop, vi)
        } else if (!(vi %in% lock_vars) && !(vj %in% lock_vars)) {
          # neither locked: default behaviour, drop j
          drop <- union(drop, vj)
        } else {
          # both locked: drop neither
        }
      }
    }
    
    keep <- setdiff(colnames(X), drop)
    X[, keep, drop = FALSE]
  }
  X_pred <- drop_constant_cols(X)
  
  if (DO_CORR_TRIM) {
    X_pred <- corr_trim(X_pred, CORR_THRESH, lock_vars = lock_vars_intersect)
  }
}
cat(sprintf("Start: X_pred has %d columns after constant-drop%s.\n\n",
            ncol(X_pred), if (DO_CORR_TRIM) " + corr-trim" else ""))

WEIGHTING_MODE <- tolower(WEIGHTING_MODE)
if (!WEIGHTING_MODE %in% c("id_guided", "uniform")) {
  stop("Unknown WEIGHTING_MODE: ", WEIGHTING_MODE)
}

w_init <- setNames(rep(W_MIN, ncol(X)), colnames(X))
w_init[colnames(X_pred)] <- 1
allow <- setNames(rep(FALSE, ncol(X)), colnames(X))
allow[colnames(X_pred)] <- TRUE
# locked variables are off-limits to optimisation
if (length(lock_vars_intersect)) {
  allow[lock_vars_intersect] <- FALSE
}

# --- Gower optimisation: single reference run + optional parallel replicates ---
if (identical(WEIGHTING_MODE, "id_guided") && GOWER_MULTI_ENABLE && GOWER_MULTI_RUNS > 1L) {
  cat(sprintf("[Gower-multi] Running %d optimisation replicates (reference run = 1, base seed = %d).\n",
              GOWER_MULTI_RUNS, SEED_GLOBAL))
  
  run_ids <- seq_len(GOWER_MULTI_RUNS)
  
  wopt_list <- progressr::with_progress({
    p <- progressr::progressor(steps = length(run_ids))
    
    FUTURE_LAPPLY(
      run_ids,
      function(r) {
        # reference run uses the exact global seeds; others are keyed perturbations
        base_seed_r <- if (r == 1L) SEED_GLOBAL else .seed_from_key(SEED_GLOBAL, paste0("gower_run_", r))
        jitter_r    <- if (r == 1L) SEED_JITTER else .seed_from_key(SEED_JITTER, paste0("gower_run_", r))
        verbose_r   <- (r == 1L)
        plot_r      <- (r == 1L)
        
        progress_cb <- if (r == 1L) {
          function(info) {
            cat(sprintf("[Gower ref] iter %d: ID=%.3f Obj=%.3f\n",
                        info$iter, info$ID, info$Obj))
          }
        } else NULL
        
        res <- optimise_gower_weights_constrained(
          X,
          init_weights = w_init,
          allow_update = allow,
          objective    = "TwoNN_all",
          w_min        = W_MIN,
          step_grid    = W_STEP_GRID,
          batch_k      = W_BATCH_K,
          batch_factor = W_BATCH_FACTOR,
          max_iter     = W_MAX_ITERS,
          n_rows_sub   = N_ROWS_SUB,
          ncores       = NCORES_PAR,
          seed_jitter  = jitter_r,
          reps_idx     = reps,
          core_idx_rep = core_idx_rep,
          lambda_l1    = GOWER_L1,
          lambda_l2    = GOWER_L2,
          verbose      = verbose_r,
          plot_progress = plot_r,
          base_seed    = base_seed_r,
          progress_fun = progress_cb
        )
        
        p(sprintf("Gower run %d/%d done", r, length(run_ids)))
        res
      },
      future.seed = TRUE
    )
  })
  
  # Use run 1 (global seed) as reference for all downstream geometry
  wopt <- wopt_list[[1L]]
} else if (identical(WEIGHTING_MODE, "id_guided")) {
  wopt_list <- NULL
  wopt <- optimise_gower_weights_constrained(
    X,
    init_weights = w_init,
    allow_update = allow,
    objective    = "TwoNN_all",
    w_min        = W_MIN,
    step_grid    = W_STEP_GRID,
    batch_k      = W_BATCH_K,
    batch_factor = W_BATCH_FACTOR,
    max_iter     = W_MAX_ITERS,
    n_rows_sub   = N_ROWS_SUB,
    ncores       = NCORES_PAR,
    seed_jitter  = SEED_JITTER,
    reps_idx     = reps,
    core_idx_rep = core_idx_rep,
    lambda_l1    = GOWER_L1,
    lambda_l2    = GOWER_L2,
    verbose      = TRUE,
    plot_progress = TRUE,
    base_seed    = SEED_GLOBAL
  )
} else {
  survivors <- colnames(X_pred)
  if (!length(survivors)) {
    stop("[weights] uniform mode retained zero predictors after preprocessing.")
  }
  
  wopt_list <- NULL
  w_full <- setNames(rep(0, ncol(X)), colnames(X))
  w_full[survivors] <- 1
  wopt <- list(
    idx_used = reps,
    history = data.frame(),
    final_ID = NA_real_,
    weights = w_full
  )
  
  cat(sprintf(
    "[weights] Uniform mode: skipping Gower optimisation; retaining %d/%d predictors with unit weights.\n",
    length(survivors), ncol(X)
  ))
}

write_gower_multi_run_coverage <- function(wopt_runs, X, reps, allow, mult = NULL) {
  if (is.null(wopt_runs) || !length(wopt_runs)) return(invisible(NULL))
  
  n_runs <- length(wopt_runs)
  row_ids <- rownames(X)
  if (is.null(row_ids) || length(row_ids) != nrow(X)) row_ids <- as.character(seq_len(nrow(X)))
  
  run_rows <- do.call(rbind, lapply(seq_along(wopt_runs), function(r) {
    idx <- as.integer(wopt_runs[[r]]$idx_used)
    idx <- idx[is.finite(idx) & idx >= 1L & idx <= nrow(X)]
    data.frame(
      run = as.integer(r),
      row_index = idx,
      participant_id = row_ids[idx],
      stringsAsFactors = FALSE
    )
  }))
  
  if (is.null(run_rows) || !nrow(run_rows)) return(invisible(NULL))
  
  write_csv(run_rows, "gower_multi_run_run_rows.csv")
  
  support_by_run <- do.call(rbind, lapply(seq_along(wopt_runs), function(r) {
    ss <- wopt_runs[[r]]$subset_support
    if (!is.data.frame(ss) || !nrow(ss)) return(NULL)
    data.frame(run = as.integer(r), ss, check.names = FALSE, stringsAsFactors = FALSE)
  }))
  if (!is.null(support_by_run) && nrow(support_by_run)) {
    write_csv(support_by_run, "gower_subset_variable_support_by_run.csv")
  }
  
  run_count <- tabulate(run_rows$row_index, nbins = nrow(X))
  run_list <- split(run_rows$run, run_rows$row_index)
  run_list_chr <- rep("", nrow(X))
  run_list_chr[as.integer(names(run_list))] <- vapply(run_list, function(z) {
    paste(sort(unique(z)), collapse = "|")
  }, character(1))
  
  mult_by_row <- rep(NA_integer_, nrow(X))
  if (!is.null(mult) && length(mult) == length(reps)) {
    mult_by_row[as.integer(reps)] <- as.integer(mult)
  }
  
  row_cov <- data.frame(
    row_index = seq_len(nrow(X)),
    participant_id = row_ids,
    in_optimiser_pool = seq_len(nrow(X)) %in% as.integer(reps),
    dedup_multiplicity = mult_by_row,
    n_runs_used = as.integer(run_count),
    prop_runs_used = as.numeric(run_count / n_runs),
    covered = run_count > 0L,
    run_list = run_list_chr,
    stringsAsFactors = FALSE
  )
  write_csv(row_cov, "gower_multi_run_row_coverage.csv")
  
  rows_per_run <- as.integer(table(factor(run_rows$run, levels = seq_len(n_runs))))
  pool_n <- length(unique(as.integer(reps)))
  requested_n <- suppressWarnings(as.integer(N_ROWS_SUB))
  if (!is.finite(requested_n) || requested_n <= 0L) requested_n <- pool_n
  requested_n <- min(requested_n, pool_n)
  expected_unique_independent <- if (pool_n > 0L) {
    pool_n * (1 - (1 - requested_n / pool_n)^n_runs)
  } else {
    NA_real_
  }
  
  pool_rows <- as.integer(reps)
  pool_counts <- run_count[pool_rows]
  summary_tbl <- data.frame(
    n_runs = n_runs,
    n_rows_total = nrow(X),
    n_rows_optimiser_pool = pool_n,
    n_rows_sub_requested = N_ROWS_SUB,
    rows_per_run_min = min(rows_per_run),
    rows_per_run_mean = mean(rows_per_run),
    rows_per_run_max = max(rows_per_run),
    fix_rep_subset = isTRUE(FIX_REP_SUBSET),
    unique_pool_rows_used = sum(pool_counts > 0L),
    prop_pool_rows_used = if (pool_n > 0L) mean(pool_counts > 0L) else NA_real_,
    uncovered_pool_rows = sum(pool_counts == 0L),
    all_pool_rows_used = all(pool_counts > 0L),
    unique_total_rows_used = sum(run_count > 0L),
    prop_total_rows_used = mean(run_count > 0L),
    all_total_rows_used = all(run_count > 0L),
    min_runs_per_covered_row = if (any(run_count > 0L)) min(run_count[run_count > 0L]) else NA_integer_,
    max_runs_per_covered_row = if (any(run_count > 0L)) max(run_count[run_count > 0L]) else NA_integer_,
    expected_unique_pool_rows_if_independent_sampling = expected_unique_independent,
    expected_prop_pool_rows_if_independent_sampling = if (pool_n > 0L) expected_unique_independent / pool_n else NA_real_,
    stringsAsFactors = FALSE
  )
  write_csv(summary_tbl, "gower_multi_run_coverage_summary.csv")
  
  all_vars <- names(wopt_runs[[1L]]$weights)
  var_cov <- do.call(rbind, lapply(all_vars, function(v) {
    weights_v <- vapply(wopt_runs, function(run) {
      ww <- run$weights
      if (v %in% names(ww)) as.numeric(ww[[v]]) else NA_real_
    }, numeric(1))
    changed_v <- vapply(wopt_runs, function(run) {
      hist <- run$history
      is.data.frame(hist) && "changed" %in% names(hist) && any(hist$changed == v, na.rm = TRUE)
    }, logical(1))
    support_v <- vapply(wopt_runs, function(run) {
      ss <- run$subset_support
      if (!is.data.frame(ss) || !"var" %in% names(ss) || !"support_ok" %in% names(ss)) return(NA)
      ix <- match(v, ss$var)
      if (is.na(ix)) NA else isTRUE(ss$support_ok[[ix]])
    }, logical(1))
    data.frame(
      var = v,
      allowed_update = isTRUE(allow[[v]]),
      n_runs_with_subset_support = sum(support_v, na.rm = TRUE),
      prop_runs_with_subset_support = mean(support_v, na.rm = TRUE),
      n_runs_with_weight = sum(is.finite(weights_v)),
      prop_runs_with_weight = mean(is.finite(weights_v)),
      n_runs_changed = sum(changed_v),
      prop_runs_changed = mean(changed_v),
      weight_min_across_runs = suppressWarnings(min(weights_v, na.rm = TRUE)),
      weight_max_across_runs = suppressWarnings(max(weights_v, na.rm = TRUE)),
      stringsAsFactors = FALSE
    )
  }))
  write_csv(var_cov, "gower_multi_run_variable_coverage.csv")
  
  cat(sprintf(
    "[Gower coverage] wrote row/run/variable coverage exports; unique optimiser-pool rows used: %d/%d across %d run(s).\n",
    summary_tbl$unique_pool_rows_used, summary_tbl$n_rows_optimiser_pool, n_runs
  ))
  
  invisible(list(row_coverage = row_cov, summary = summary_tbl, variable_coverage = var_cov))
}

if (identical(WEIGHTING_MODE, "id_guided")) {
  coverage_runs <- if (!is.null(wopt_list) && length(wopt_list)) wopt_list else list(wopt)
  write_gower_multi_run_coverage(coverage_runs, X = X, reps = reps, allow = allow, mult = mult)
}

# --- DIAGNOSTIC: pre-selection geometry on optimiser subset ---
ix <- wopt$idx_used

PX_pre <- prep_X_for_gower(X[ix, , drop=FALSE],
                           rare_prop = RARE_LEVEL_MIN_PROP,
                           do_jitter = TRUE)  # this matches optimiser call

Xsub_pre <- PX_pre$X
typ_pre  <- PX_pre$type

D_un_pre <- gower_dist(Xsub_pre, type_list = typ_pre, weights = rep(1, ncol(Xsub_pre)))
w_al_pre <- wopt$weights[colnames(Xsub_pre)]
D_op_pre <- gower_dist(Xsub_pre, type_list = typ_pre, weights = w_al_pre)

dist_health_console(D_un_pre, "unweighted PRE")
dist_health_console(D_op_pre, "optimised PRE")

# Select survivor vars via knee on weight tail, with optional stability filtering
if (identical(WEIGHTING_MODE, "id_guided") && !is.null(wopt_list) && length(wopt_list) > 1L) {
  var_names <- names(wopt$weights)
  n_runs <- length(wopt_list)
  
  surv_mat      <- matrix(0L, nrow = length(var_names), ncol = n_runs,
                          dimnames = list(var_names, paste0("run", seq_len(n_runs))))
  sel_mat       <- matrix(0L, nrow = length(var_names), ncol = n_runs,
                          dimnames = list(var_names, paste0("run", seq_len(n_runs))))
  thr_tail_vec  <- rep(NA_real_, n_runs)
  final_ID_vec  <- rep(NA_real_, n_runs)
  
  # Reference run: this is the one that drives all downstream geometry and plots
  sel_ref <- survivors_from_weights(
    w          = wopt$weights,
    w_min      = W_MIN,
    thr        = KNEE_THR,
    make_plot  = TRUE,
    plot_file  = "FIG_weight_curve_knee.png"
  )
  sel_mat[sel_ref$survivors, 1L] <- 1L
  surv_mat[sel_ref$active_survivors, 1L] <- 1L
  thr_tail_vec[1L] <- sel_ref$thr_tail
  final_ID_vec[1L] <- wopt$final_ID
  
  if (n_runs > 1L) {
    for (r in 2:n_runs) {
      wr <- wopt_list[[r]]$weights
      sel_r <- survivors_from_weights(
        w         = wr,
        w_min     = W_MIN,
        thr       = KNEE_THR,
        make_plot = FALSE
      )
      sel_mat[sel_r$survivors, r] <- 1L
      surv_mat[sel_r$active_survivors, r] <- 1L
      thr_tail_vec[r] <- sel_r$thr_tail
      final_ID_vec[r] <- wopt_list[[r]]$final_ID
    }
  }
  
  surv_count <- rowSums(surv_mat)
  surv_prop  <- surv_count / n_runs
  sel_count <- rowSums(sel_mat)
  sel_prop <- sel_count / n_runs
  
  # build var x run weight matrix (all runs, all vars)
  Wmat <- vapply(seq_len(n_runs), function(r) {
    wr <- wopt_list[[r]]$weights
    # ensure consistent ordering, fill with NA for any missing names (should not happen)
    wr[var_names]
  }, numeric(length(var_names)))
  # vapply gives matrix with rows=vars, cols=runs
  dimnames(Wmat) <- list(var_names, paste0("run", seq_len(n_runs)))
  
  # basic summary stats of weight across runs
  weight_mean <- rowMeans(Wmat, na.rm = TRUE)
  weight_sd   <- apply(Wmat, 1L, stats::sd, na.rm = TRUE)
  weight_min  <- apply(Wmat, 1L, min, na.rm = TRUE)
  weight_max  <- apply(Wmat, 1L, max, na.rm = TRUE)
  
  stability_tbl <- data.frame(
    var         = var_names,
    count       = as.integer(surv_count),
    prop        = as.numeric(surv_prop),
    count_selected = as.integer(sel_count),
    prop_selected = as.numeric(sel_prop),
    selected_ref = var_names %in% sel_ref$survivors,
    in_ref      = var_names %in% sel_ref$active_survivors,
    weight_ref  = as.numeric(wopt$weights[var_names]),
    weight_mean = as.numeric(weight_mean),
    weight_sd   = as.numeric(weight_sd),
    weight_min  = as.numeric(weight_min),
    weight_max  = as.numeric(weight_max),
    stringsAsFactors = FALSE
  )
  stability_tbl <- stability_tbl[order(-stability_tbl$prop, -stability_tbl$weight_mean), ]
  
  write_csv(stability_tbl, "gower_survivor_stability.csv")
  
  w_wide <- data.frame(
    var = rownames(Wmat),
    Wmat,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  write_csv(w_wide, "gower_weights_by_run_wide.csv")
  
  # Save a compact summary of all runs for later inspection
  saveRDS(list(
    n_runs   = n_runs,
    ref_run  = 1L,
    min_prop = GOWER_MULTI_MIN_PROP,
    final_ID = final_ID_vec,
    thr_tail = thr_tail_vec,
    weights  = lapply(wopt_list, function(z) z$weights)
  ), file = "gower_multi_runs_summary.rds")
  
  # Survivors purely by multi-run stability
  active_floor <- W_MIN + sel_ref$active_eps
  robust_mask <- stability_tbl$prop >= GOWER_MULTI_MIN_PROP & stability_tbl$weight_max > active_floor
  if (!any(robust_mask)) {
    warning(sprintf(
      "[Gower-multi] No active variable met prop >= %.2f; falling back to active reference survivors.",
      GOWER_MULTI_MIN_PROP
    ))
    robust_set <- sel_ref$active_survivors
    if (!length(robust_set)) robust_set <- sel_ref$survivors
    robust_mask <- stability_tbl$var %in% robust_set
  } else {
    robust_set <- stability_tbl$var[robust_mask]
  }
  
  cat(sprintf("[Gower-multi] Reference selected = %d | reference active = %d | robust survivors = %d (min_prop=%.2f, runs=%d).\n",
              length(sel_ref$survivors),
              length(sel_ref$active_survivors),
              length(robust_set),
              GOWER_MULTI_MIN_PROP,
              n_runs))
  
  # Aggregate weights across runs (expected weight)
  w_mean <- stability_tbl$weight_mean
  names(w_mean) <- stability_tbl$var
  
  # Multi-run aggregate weights: prop × mean-active logic
  # w_multi_raw_i = weight_mean_i - W_MIN * (1 - prop_i)
  w_multi_raw <- with(stability_tbl,
                      weight_mean - W_MIN * (1 - prop))
  
  # Guard against numerical junk
  w_multi_raw[!is.finite(w_multi_raw)] <- 0
  w_multi_raw <- pmax(0, w_multi_raw)
  
  # Renormalise so the top robust variable has weight 1
  scale_fac <- max(w_multi_raw[robust_mask], na.rm = TRUE)
  if (!is.finite(scale_fac) || scale_fac <= 0) scale_fac <- 1
  
  w_multi <- w_multi_raw / scale_fac
  
  # Full weight vector for all vars (useful for exporting diagnostics)
  w_full <- setNames(as.numeric(w_multi), stability_tbl$var)
  
  # Final survivor set is robustness-based
  survivors <- robust_set
} else if (identical(WEIGHTING_MODE, "id_guided")) {
  # Single optimiser run: original behaviour
  sel <- survivors_from_weights(
    w          = wopt$weights,
    w_min      = W_MIN,
    thr        = KNEE_THR,
    make_plot  = TRUE,
    plot_file  = "FIG_weight_curve_knee.png"
  )
  w_full    <- sel$w_clamped
  survivors <- sel$survivors
}

X    <- X[, survivors, drop = FALSE]
w_all <- w_full[survivors]
cat(sprintf("[weights] Mode=%s | retained predictors: p=%d\n", WEIGHTING_MODE, ncol(X)))

# Recalculate final weighted distances on the requested diagnostic subset.
PX <- prep_X_for_gower(X, rare_prop = RARE_LEVEL_MIN_PROP, do_jitter = TRUE, seed = SEED_JITTER)
Xg <- PX$X
type_g <- PX$type
w_use <- w_all[colnames(Xg)]

if (!FINAL_DIAG_MODE %in% c("full", "reps", "sample_reps", "sample_all", "none")) {
  stop("Unknown FINAL_DIAG_MODE: ", FINAL_DIAG_MODE)
}

diag_cap <- as.integer(DIAG_N_MAX %||% 4000L)
if (!is.finite(diag_cap) || diag_cap <= 0L) diag_cap <- nrow(Xg)

diag_idx <- switch(
  FINAL_DIAG_MODE,
  "none" = integer(0),
  "full" = seq_len(nrow(Xg)),
  "reps" = {
    if (exists("reps") && length(reps) > 0L) reps else seq_len(nrow(Xg))
  },
  "sample_reps" = {
    pool <- if (exists("reps") && length(reps) > 0L) reps else seq_len(nrow(Xg))
    n_take <- min(length(pool), diag_cap)
    if (length(pool) > n_take) {
      warn_diag_subsample("Final diag", n_take, length(pool), diag_cap, "representatives")
      .with_seed(.seed_from_key(SEED_GLOBAL, "final_diag_sample_reps"), sample(pool, n_take))
    } else {
      pool
    }
  },
  "sample_all" = {
    pool <- seq_len(nrow(Xg))
    n_take <- min(length(pool), diag_cap)
    if (length(pool) > n_take) {
      warn_diag_subsample("Final diag", n_take, length(pool), diag_cap, "rows")
      .with_seed(.seed_from_key(SEED_GLOBAL, "final_diag_sample_all"), sample(pool, n_take))
    } else {
      pool
    }
  }
)

if (length(diag_idx) >= 3L) {
  Xg_diag <- Xg[diag_idx, , drop = FALSE]
  w_use_diag <- w_use[colnames(Xg_diag)]
  stopifnot(length(w_use_diag) == ncol(Xg_diag), all(is.finite(w_use_diag)))
  
  D_final <- cluster::daisy(Xg_diag, metric = "gower", type = type_g, weights = w_use_diag)
  D_un_final <- cluster::daisy(Xg_diag, metric = "gower", type = type_g, weights = rep(1, ncol(Xg_diag)))
  ID_all <- twonn_id_from_dist(D_final)
  core_ix <- twonn_core_by_slope(D_final)
  
  DmF <- as.matrix(D_final)
  diag(DmF) <- Inf
  
  ID_core <- if (length(core_ix) >= 3L) {
    twonn_id_from_dist(as.dist(DmF[core_ix, core_ix, drop = FALSE]))
  } else {
    NA_real_
  }
  ID_LB <- if (length(core_ix) >= 15L) {
    lb_mle_id(DmF[core_ix, core_ix, drop = FALSE], 5, 15)
  } else {
    NA_real_
  }
  
  cat(sprintf(
    "[Constrained %s | diag=%s | n=%d] TwoNN_all=%.2f | TwoNN_core=%.2f | LB_core=%.2f (n_core=%d, p_active=%d)\n",
    toupper(PREF_TARGET), FINAL_DIAG_MODE, nrow(Xg_diag),
    ID_all, ID_core, ID_LB, length(core_ix), ncol(Xg)
  ))
  
  dist_health_console(D_un_final, sprintf("unweighted FINAL (%s)", FINAL_DIAG_MODE))
  dist_health_console(D_final, sprintf("weighted FINAL (%s)", FINAL_DIAG_MODE))
  
  for (k in c(10, 15, 30)) {
    nn_un <- knn_from_dist(D_un_final, k = k)
    nn_op <- knn_from_dist(D_final, k = k)
    cat(sprintf("[kNN overlap FINAL] k=%d | mean Jaccard = %.3f\n", k, mean_jaccard_knn(nn_un, nn_op)))
  }
  
  Zp <- as.matrix(Xg_diag)
  Zr <- apply(Zp, 2, function(v) round(as.numeric(v), 6))
  key <- apply(Zr, 1, paste, collapse = "|")
  tb <- sort(table(key), decreasing = TRUE)
  
  cat(sprintf("[FINAL quant] unique patterns = %d / %d\n", length(tb), length(key)))
  cat(sprintf("[FINAL quant] max multiplicity = %d | mean multiplicity = %.2f\n",
              as.integer(tb[1]), mean(as.numeric(tb))))
  cat(sprintf("[FINAL quant] prop in non-singletons = %.3f\n",
              mean(key %in% names(tb)[tb > 1])))
} else {
  cat(sprintf(
    "[Constrained %s | diag=%s] Final distance-based diagnostics skipped.\n",
    toupper(PREF_TARGET), FINAL_DIAG_MODE
  ))
}

# ==============================================================================
# 5. PCA, Whitening, and Residuals
# ==============================================================================

# Encode
Xenc <- design_with_map(X, nominal_missing = PCA_NOMINAL_MISSING)
varmap_full <- attr(Xenc, "varmap")
names(varmap_full) <- colnames(Xenc)
levelmap_full <- attr(Xenc, "levelmap")
names(levelmap_full) <- colnames(Xenc)
encoding_imputation <- attr(Xenc, "imputation")
if (!is.null(encoding_imputation) && nrow(encoding_imputation)) {
  write_csv(encoding_imputation, "encoding_imputation.csv")
}

# Standardise + drop degenerate encoded cols
Z0 <- scale(Xenc, center = TRUE, scale = TRUE)
sds0 <- apply(Z0, 2, sd, na.rm = TRUE)
keep_cols <- is.finite(sds0) & sds0 > 1e-12
Z0 <- Z0[, keep_cols, drop = FALSE]
Z0[!is.finite(Z0)] <- 0
if (!ncol(Z0)) {
  stop("[PCA] No encoded columns remain after standardisation/variance filtering.")
}

# Map for surviving encoded columns
varmap <- varmap_full[colnames(Z0)]
levelmap <- levelmap_full[colnames(Z0)]
nominal_vars_encoded <- intersect(names(X)[vapply(X, is.factor, logical(1))], unique(varmap))
encoded_counts <- table(varmap)
bad_nominal_encoding <- nominal_vars_encoded[encoded_counts[nominal_vars_encoded] < 2L]
if (length(bad_nominal_encoding)) {
  stop("[PCA] Nominal variables were not expanded to at least two response-level columns: ",
       paste(head(bad_nominal_encoding, 20L), collapse = ", "))
}
cat(sprintf(
  "[PCA] Encoded %d source variables as %d columns; nominal sources contribute one column per observed response level.\n",
  length(unique(varmap)), ncol(Z0)
))

# Allocate weights AFTER filtering so per-item totals are preserved
w_enc <- setNames(numeric(ncol(Z0)), colnames(Z0))
for (nm in unique(varmap)) {
  idx <- which(varmap == nm)
  wj <- w_all[nm]
  if (!is.finite(wj)) wj <- 1
  w_enc[idx] <- wj / length(idx)
}

# Apply sqrt-weights
Xenc_w <- sweep(Z0, 2, sqrt(pmax(w_enc, 0)), "*")

tot_by_item <- tapply(w_enc, varmap, sum)
common <- intersect(names(tot_by_item), names(w_all))
weight_error <- max(abs(tot_by_item[common] - w_all[common]), na.rm = TRUE)
if (!is.finite(weight_error) || weight_error > 1e-10) {
  stop(sprintf("[PCA] Encoded weight allocation drifted from item weights (max abs error %.3g).", weight_error))
}

# PCA input: do not call scale() on Xenc_w. A second standardisation would erase
# the sqrt-weighting just applied above.
Z <- Xenc_w
if (!isTRUE(all.equal(unclass(Z), unclass(Xenc_w), check.attributes = FALSE))) {
  stop("[PCA] Weighted PCA input was unexpectedly modified before decomposition.")
}

BASE_DECOMP_METHOD <- tolower(as.character(BASE_DECOMP_METHOD))
if (!BASE_DECOMP_METHOD %in% c("robust_pca", "pca")) {
  stop("Unknown BASE_DECOMP_METHOD: ", BASE_DECOMP_METHOD)
}

max_k <- min(nrow(Z) - 1L, ncol(Z))

cat(sprintf("[PCA dbg] nrow(Z)=%d ncol(Z)=%d max_k=%d | requested BASE_DIM=%d M_STAR_FIXED=%d\n",
            nrow(Z), ncol(Z), max_k, as.integer(BASE_DIM), as.integer(M_STAR_FIXED)))

if (max_k < BASE_DIM) {
  stop(sprintf(
    "[PCA] Cannot realise BASE_DIM=%d: only %d feasible components after encoding/standardisation/weighting/filters (nrow=%d, ncol=%d).",
    as.integer(BASE_DIM), as.integer(max_k), nrow(Z), ncol(Z)
  ))
}

m_star <- max(as.integer(M_STAR_FIXED), as.integer(BASE_DIM), 2L)
k_eff  <- min(max_k, m_star)

if (identical(BASE_DECOMP_METHOD, "robust_pca")) {
  mcd_ctrl <- rrcov::CovControlMcd(
    nsamp = "deterministic",
    seed  = SEED_GLOBAL
  )
  
  rp <- rrcov::PcaHubert(
    x           = Z,
    k           = k_eff,
    kmax        = k_eff,
    scale       = FALSE,
    signflip    = TRUE,
    cov.control = mcd_ctrl
  )
  
  Base <- rp@scores[, seq_len(k_eff), drop = FALSE]
  base_loadings <- rp@loadings[, seq_len(k_eff), drop = FALSE]
  base_spectrum <- rp@eigenvalues[seq_len(k_eff)]
} else {
  pc <- stats::prcomp(Z, center = FALSE, scale. = FALSE, rank. = k_eff)
  Base <- pc$x[, seq_len(k_eff), drop = FALSE]
  base_loadings <- pc$rotation[, seq_len(k_eff), drop = FALSE]
  base_spectrum <- pc$sdev[seq_len(k_eff)]^2
}

base_total_var <- sum(apply(Z, 2, stats::var), na.rm = TRUE)
base_explained_variance_ratio <- base_spectrum / pmax(base_total_var, 1e-12)
colnames(Base) <- paste0("b", seq_len(ncol(Base)))
colnames(base_loadings) <- paste0("b", seq_len(ncol(base_loadings)))
cat(sprintf("[Base] Decomposition=%s. Proceeding to Whitening...\n", BASE_DECOMP_METHOD))

# ---- 2D/3D snapshot depending on BASE_DIM ----
if (ncol(Base) >= 3L) {
  Base_A <- as.data.frame(Base[, 1:3, drop = FALSE])
  names(Base_A) <- c("b1", "b2", "b3")
} else {
  Base_A <- as.data.frame(Base[, 1:2, drop = FALSE])
  names(Base_A) <- c("b1", "b2")
}

if (BASE_DIM == 3L && ncol(Base) < 3L) {
  stop("[sanity] 3D requested, but Base is only ", ncol(Base), "D.")
}

# --- Whitening for Neighbour Search ---
S <- stats::cov(Base)
U <- try(chol(S + diag(1e-8, ncol(Base))), silent = TRUE)
if (inherits(U, "try-error")) {
  eig <- eigen(S, symmetric = TRUE)
  U <- t(eig$vectors %*% diag(sqrt(pmax(eig$values, 1e-8))) %*% t(eig$vectors))
}
Base_w <- Base %*% solve(U)

KS_RESIDUAL <- c(6, 8, 10, 12, 16, 20)
nb_list <- setNames(lapply(KS_RESIDUAL, function(k) {
  RANN::nn2(Base_w, Base_w, k = pmin(k + 1L, nrow(Base_w)))$nn.idx[, -1L, drop = FALSE]
}), paste0("k", KS_RESIDUAL))

# --- Residual Extraction  ---
make_strat_folds <- function(y, K, seed = SEED_GLOBAL, key = NULL) {
  y <- as.integer(y)
  n <- length(y)
  K <- as.integer(K[1])
  if (!is.finite(K) || K < 2L) K <- 2L
  if (K > n) K <- n
  
  s <- if (is.null(key)) as.integer(seed) else .seed_from_key(seed, paste0("strat|", key))
  
  .with_seed(s, {
    folds <- integer(n)
    idx0 <- which(y == 0)
    idx1 <- which(y == 1)
    if (length(idx0)) folds[idx0] <- sample(rep(seq_len(K), length.out = length(idx0)))
    if (length(idx1)) folds[idx1] <- sample(rep(seq_len(K), length.out = length(idx1)))
    folds
  })
}

choose_K <- function(y, K_target = CV_FOLDS, min_per_class = 8) {
  y <- as.integer(y)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  max(2, min(K_target, floor(n1 / min_per_class), floor(n0 / min_per_class)))
}

# --- Unsupervised folds (no diagnoses required) ---
n <- nrow(X)
K_fold <- max(2L, min(CV_FOLDS, n))
fold_id <- .with_seed(SEED_GLOBAL, sample(rep(1:K_fold, length.out = n)))

knn_graph <- function(X, k = 15) {
  idx <- RANN::nn2(X, X, k = pmin(k + 1L, nrow(X)))$nn.idx[, -1L, drop = FALSE]
  i <- rep(seq_len(nrow(X)), each = ncol(idx))
  j <- as.vector(idx)
  g <- igraph::graph_from_edgelist(cbind(i, j), directed = FALSE)
  igraph::simplify(g)
}

if (isTRUE(DO_RESIDUALISATION)) {
E <- residualise_foldsafe(Xenc_w, Base, folds = fold_id, k_gam = 6)
E_scaled <- scale(E, center = TRUE, scale = TRUE)

XR <- E_scaled
XR[!is.finite(XR)] <- 0

cat(sprintf("[Residuals] XR matrix: %d rows × %d columns (post-OOF, scaled)\n", nrow(XR), ncol(XR)))

# ==============================================================================
# 6. Item Diagnostics (Base vs Resid Roles)
# ==============================================================================

Z0_std <- Z0
vars_diag <- unique(varmap)

# Helper for KNN K selection
choose_k_nb <- function(e, nb_list, folds = CV_FOLDS, seed = SEED_GLOBAL, key = NULL) {
  r2s <- vapply(names(nb_list), function(nm_nb) {
    r2_residual_cv(e, nb_list[[nm_nb]], folds = folds, seed = seed, key = paste0(key, "|", nm_nb))
  }, numeric(1))
  ix <- which.max(r2s)
  list(k = as.integer(sub("^k", "", names(nb_list)[ix])), R2_cv = r2s[ix])
}

RESIDUAL_PERM_B <- N_PERM
MIN_SD_ITEM <- 1e-6

if (isTRUE(RUN_ITEM_DIAGNOSTICS)) {
cat(sprintf("[Item Roles] Starting diagnostics on %d items...\n", length(vars_diag)))

roles_rows <- lapply(vars_diag, function(nm) {
  # Score item (Original)
  v <- score_item_base(nm, Z0_std, varmap)
  if (all(is.na(v)) || stats::sd(v, na.rm = TRUE) < MIN_SD_ITEM) return(NULL)
  
  # Score residual
  e_item <- e_from_E(nm, E_scaled, Z0_std, varmap)
  if (all(is.na(e_item)) || stats::sd(e_item, na.rm = TRUE) < MIN_SD_ITEM) return(NULL)
  
  # Metrics
  R2_base <- r2_base_linear(Base, v)
  sel <- choose_k_nb(e_item, nb_list, folds = CV_FOLDS, seed = SEED_GLOBAL)
  R2_resid <- sel$R2_cv
  
  data.frame(
    var = nm,
    R2_base = R2_base,
    R2_residual = R2_resid,
    k_residual = sel$k,
    stringsAsFactors = FALSE
  )
})

roles_df <- dplyr::bind_rows(roles_rows)

if (nrow(roles_df) > 0) {
  # Classify items
  ES_BASE <- 0.08      
  ES_RESIDUAL <- 0.05
  
  roles_df$role <- dplyr::case_when(
    roles_df$R2_base >= ES_BASE & roles_df$R2_residual >= ES_RESIDUAL ~ "mixed",
    roles_df$R2_base >= ES_BASE ~ "base-aligned",
    roles_df$R2_residual >= ES_RESIDUAL ~ "residual-structured",
    TRUE ~ "weak"
  )
  
  roles_df$role_final <- factor(roles_df$role, levels = c("base-aligned", "mixed", "residual-structured", "weak"))
  roles_df <- roles_df[order(-pmax(roles_df$R2_base, roles_df$R2_residual), roles_df$var), ]
  
  write_csv(roles_df, "predictive_item_roles_diagnostics.csv")
  cat(sprintf("[roles] wrote: %s  (p=%d items processed)\n", "predictive_item_roles_diagnostics.csv", nrow(roles_df)))
} else {
  warning("[roles] No valid items found for roles analysis.")
}
} else {
  message("[roles] Diagnostic-only parent-item PC1 summaries skipped (RUN_ITEM_DIAGNOSTICS=FALSE).")
}

# ==============================================================================
# 7. Resid-Only Decomposition & Clustering
# ==============================================================================

stopifnot(exists("E"), is.matrix(E), nrow(E) >= 3)

Ef <- E
Ef[!is.finite(Ef)] <- 0

n  <- nrow(Ef)
pE <- ncol(Ef)
if (pE < 2L || n < 4L) stop("[Resid-only] insufficient columns/rows in E.")

RESIDUAL_BASE_MAX <- min(6L, pE, n - 1L)

# PCA on covariance of weighted residuals (centre only, no extra scaling)
pc_f <- prcomp(Ef, center = TRUE, scale. = FALSE, rank. = max(2L, RESIDUAL_BASE_MAX))

Bprime_all <- pc_f$x[, 1:RESIDUAL_BASE_MAX, drop = FALSE]
colnames(Bprime_all) <- paste0("f", seq_len(ncol(Bprime_all)))
f3 <- Bprime_all[, 3]

Base_b1b2 <- Base[, 1:2, drop = FALSE]
Base_b1f3 <- cbind(b1 = Base[,1], b2 = as.numeric(f3))
Base_b2f3 <- cbind(b1 = Base[,2], b2 = as.numeric(f3))

oof_R2_surface_plus_f3 <- function(v, Base_b1b2, f3,
                                   K_target = 5,
                                   k2d = 10, k1d = 6,
                                   fold_id = NULL,
                                   seed = SEED_GLOBAL) {
  
  n <- length(v)
  if (is.null(fold_id)) {
    fold_id <- .with_seed(seed, sample(rep(1:K_target, length.out = n)))
  } else {
    stopifnot(length(fold_id) == n)
  }
  
  pred_add <- pred_full <- rep(NA_real_, n)
  ctrl <- mgcv::gam.control(maxit = 100, nthreads = 1)
  
  for (k in sort(unique(fold_id))) {
    tr <- which(fold_id != k); te <- which(fold_id == k)
    dtr <- data.frame(v = v[tr], b1 = Base_b1b2[tr,1], b2 = Base_b1b2[tr,2], f3 = f3[tr])
    dte <- data.frame(b1 = Base_b1b2[te,1], b2 = Base_b1b2[te,2], f3 = f3[te])
    
    g_add  <- try(mgcv::gam(v ~ s(b1,b2,k=k2d), data=dtr, method="REML", control=ctrl), silent=TRUE)
    g_full <- try(mgcv::gam(v ~ s(b1,b2,k=k2d) + s(f3,k=k1d), data=dtr, method="REML", control=ctrl), silent=TRUE)
    
    mu <- mean(dtr$v, na.rm=TRUE)
    if (!inherits(g_add,"try-error"))  { pa <- try(predict(g_add,  dte), silent=TRUE); if(!inherits(pa,"try-error")) pred_add[te]  <- pa }
    if (!inherits(g_full,"try-error")) { pf <- try(predict(g_full, dte), silent=TRUE); if(!inherits(pf,"try-error")) pred_full[te] <- pf }
    
    pred_add[te][!is.finite(pred_add[te])]   <- mu
    pred_full[te][!is.finite(pred_full[te])] <- mu
  }
  
  ve <- var(v, na.rm=TRUE)
  if (!is.finite(ve) || ve <= 0) return(c(R2_add=NA, R2_full=NA, dR2=NA))
  R2_add  <- max(0, 1 - mean((v - pred_add)^2,  na.rm=TRUE) / ve)
  R2_full <- max(0, 1 - mean((v - pred_full)^2, na.rm=TRUE) / ve)
  c(R2_add = R2_add, R2_full = R2_full, dR2 = R2_full - R2_add)
}

pick_m_via_tc <- function(Xhigh, Xlow_all, ks = 10:30, mmax = ncol(Xlow_all), lambda = 0.02) {
  trust_cont <- function(high, low, ks = 10:30) {
    high <- as.matrix(high)
    low <- as.matrix(low)
    stopifnot(nrow(high) == nrow(low))
    n <- nrow(high)
    Dh <- as.matrix(stats::dist(high))
    diag(Dh) <- Inf
    Dl <- as.matrix(stats::dist(low))
    diag(Dl) <- Inf
    
    rf <- function(D) {
      R <- matrix(0L, n, n)
      for (i in 1:n) {
        r <- D[i, ]
        ord <- order(r)
        R[i, ord] <- seq_len(n)
      }
      R
    }
    Rh <- rf(Dh)
    Rl <- rf(Dl)
    
    res <- lapply(ks, function(k) {
      H <- t(apply(Rh, 1, function(r) order(r)[1:k]))
      L <- t(apply(Rl, 1, function(r) order(r)[1:k]))
      Tsum <- 0
      Csum <- 0
      for (i in 1:n) {
        U <- setdiff(L[i, ], H[i, ])
        if (length(U)) Tsum <- Tsum + sum(pmax(Rh[i, U] - k, 0))
        V <- setdiff(H[i, ], L[i, ])
        if (length(V)) Csum <- Csum + sum(pmax(Rl[i, V] - k, 0))
      }
      denom <- n * k * (2 * n - 3 * k - 1)
      data.frame(k = k, Trust = 1 - (2 / denom) * Tsum, Continuity = 1 - (2 / denom) * Csum)
    })
    do.call(rbind, res)
  }
  
  vals <- lapply(2:mmax, function(m) {
    tc <- trust_cont(Xhigh, Xlow_all[, 1:m, drop = FALSE], ks)
    data.frame(m = m, T = mean(tc$Trust), C = mean(tc$Continuity))
  }) |> dplyr::bind_rows()
  
  vals$score <- with(vals, (T + C) - lambda * (m - min(m)))
  vals$m[which.max(vals$score)]
}

m_f <- pick_m_via_tc(Ef, Bprime_all, ks = 10:30, mmax = ncol(Bprime_all), lambda = 0.02)
Bprime <- Bprime_all[, 1:m_f, drop = FALSE]
colnames(Bprime) <- paste0("f", seq_len(ncol(Bprime)))

# OOF residuals of E on B' -> F' (linear)
Kf <- max(2L, min(CV_FOLDS, n))
folds_f <- .with_seed(SEED_PRED, sample(rep(1:Kf, length.out = n)))

residualise_linear_oof <- function(Y, X, folds) {
  n <- nrow(Y)
  p <- ncol(Y)
  R <- matrix(NA_real_, n, p, dimnames = list(rownames(Y), colnames(Y)))
  Xdf <- as.data.frame(X)
  
  for (j in seq_len(p)) {
    y <- Y[, j]
    if (all(!is.finite(y))) next
    for (k in sort(unique(folds))) {
      tr <- which(folds != k)
      te <- which(folds == k)
      fit <- try(stats::lm(y ~ ., data = cbind.data.frame(y = y[tr], Xdf[tr, , drop = FALSE])), silent = TRUE)
      if (inherits(fit, "try-error")) next
      mu <- as.numeric(predict(fit, newdata = Xdf[te, , drop = FALSE]))
      R[te, j] <- y[te] - mu
    }
  }
  R[, colSums(is.finite(R)) > 0, drop = FALSE]
}
Fprime <- residualise_linear_oof(Ef, Bprime, folds_f)

# ID diagnostics in residual spaces
D_Bprime <- stats::dist(Bprime)
MB <- as.matrix(D_Bprime)
diag(MB) <- Inf
core_B <- core_band_idx(D_Bprime, k = CORE_KNN_K, band = CORE_BAND)
ID_B_all <- twonn_id_from_dist(D_Bprime)
ID_B_core <- if (length(core_B) >= 3) twonn_id_from_dist(stats::as.dist(MB[core_B, core_B])) else NA_real_
ID_B_LB <- if (length(core_B) >= 3) lb_mle_id(MB[core_B, core_B, drop = FALSE], 5, 15) else NA_real_

if (ncol(Fprime) >= 2) {
  D_Fprime <- stats::dist(Fprime)
  MFp <- as.matrix(D_Fprime)
  diag(MFp) <- Inf
  core_Fp <- core_band_idx(D_Fprime, k = CORE_KNN_K, band = CORE_BAND)
  ID_Fp_all <- twonn_id_from_dist(D_Fprime)
  ID_Fp_core <- if (length(core_Fp) >= 3) twonn_id_from_dist(stats::as.dist(MFp[core_Fp, core_Fp])) else NA_real_
  ID_Fp_LB <- if (length(core_Fp) >= 3) lb_mle_id(MFp[core_Fp, core_Fp, drop = FALSE], 5, 15) else NA_real_
} else {
  ID_Fp_all <- ID_Fp_core <- ID_Fp_LB <- NA_real_
}

cat(sprintf("[Resid-only] ID(B'): TwoNN_all=%.2f | TwoNN_core=%.2f | LB_core=%.2f (n_core=%d)\n",
            ID_B_all, ID_B_core, ID_B_LB, length(core_B)))
cat(sprintf("[Resid-only] ID(F'): TwoNN_all=%.2f | TwoNN_core=%.2f | LB_core=%.2f (n_core=%s)\n",
            ID_Fp_all, ID_Fp_core, ID_Fp_LB, if (exists("core_Fp")) length(core_Fp) else "NA"))

# Clustering in B' (Louvain on kNN graph)
gF <- knn_graph(Bprime, k = 15)
clF <- igraph::cluster_louvain(gF)$membership

saveRDS(list(
  Bprime = Bprime,
  Fprime = Fprime,
  m_f = m_f,
  ID_Bprime = c(all = ID_B_all, core = ID_B_core, LB = ID_B_LB),
  ID_Fprime = c(all = ID_Fp_all, core = ID_Fp_core, LB = ID_Fp_LB),
  clusters = clF
), file = file.path(OUTPUTS_DIR, "residual_only_summary.rds"))

cat("[Resid-only] wrote:", file.path(OUTPUTS_DIR, "residual_only_summary.rds"), "\n")
} else {
  message("[Residuals] DO_RESIDUALISATION=FALSE; skipping E/E_scaled, item residual diagnostics, Bprime, and Fprime.")
  for (.res_obj in c("E", "E_scaled", "XR", "Bprime", "Bprime_all", "Fprime", "emb_residual")) {
    if (exists(.res_obj, inherits = FALSE)) rm(list = .res_obj)
  }
  gF <- knn_graph(Base, k = 15)
  clF <- igraph::cluster_louvain(gF)$membership
}

# ==============================================================================
# 8. Predictive Diagnostics (OOF & Interactions)
# ==============================================================================
if (DX_AVAILABLE) {
  # Align dx labels to Base
  ids_base <- rownames(Base)
  stopifnot(length(ids_base) == nrow(Base), all(nzchar(ids_base)))
  
  DxW <- diag_wide_full %>%
    dplyr::transmute(participant_id = 
                       trimws(as.character(participant_id)),
                     dplyr::across(-participant_id, ~ as.integer(.x))) %>%
    dplyr::right_join(tibble::tibble(participant_id = ids_base, .row = seq_along(ids_base)), by = "participant_id") %>%
    dplyr::arrange(.row) %>%
    dplyr::select(-participant_id, -.row)
  DxW[is.na(DxW)] <- 0L
  
  prev <- colMeans(DxW > 0, na.rm = TRUE)
  cases <- colSums(DxW > 0, na.rm = TRUE)
  keep_dx <- names(prev)[(prev >= DX_PREV_MIN) & (cases >= DX_CASES_MIN)]
  
  if (!length(keep_dx)) {
    warning("[dx] No diagnosis passes thresholds; skipping predictive diagnostics.")
    keep_dx <- intersect(names(prev), names(prev)[cases > 0 & (1 - prev) * nrow(DxW) > 0])
  }
  DxW_A <- as.data.frame(DxW)[, keep_dx, drop = FALSE]
  rownames(DxW_A) <- ids_base
  
  # --- Diagnosis probability fields over Base ---
  dir.create(paste0(OUTPUTS_DIR, "/base_prob_grids"), showWarnings = FALSE)
  
  grid_from_base <- function(Base, nx = 140, ny = 140, pad = 0.05, q = c(0.01, 0.99)) {
    rx <- quantile(Base[, 1], q, na.rm = TRUE)
    ry <- quantile(Base[, 2], q, na.rm = TRUE)
    wx <- diff(rx); wy <- diff(ry)
    xs <- seq(rx[1] - pad * wx, rx[2] + pad * wx, length.out = nx)
    ys <- seq(ry[1] - pad * wy, ry[2] + pad * wy, length.out = ny)
    as.matrix(expand.grid(b1 = xs, b2 = ys))
  }
  
  predict_dx_surface <- function(y, Base, gridXY, k_gam = 30) {
    df <- data.frame(y = as.integer(y), b1 = Base[, 1], b2 = Base[, 2])
    fit <- try(mgcv::gam(y ~ s(b1, b2, k = k_gam), data = df, family = binomial(), method = "REML"), silent = TRUE)
    if (inherits(fit, "try-error")) return(NULL)
    
    pg <- try(as.numeric(predict(fit,
                                 newdata = data.frame(
                                   b1 = gridXY[, 1],
                                   b2 = gridXY[, 2]),
                                 type = "response")),
              silent = TRUE)
    
    if (inherits(pg, "try-error")) return(NULL)
    pmin(pmax(pg, 1e-6), 1 - 1e-6)
  }
  
  make_fname_safe <- function(x, max_len = 80) {
    s <- iconv(as.character(x), to = "ASCII//TRANSLIT")
    s <- tolower(s)
    s <- gsub("[^a-z0-9]+", "_", s)
    s <- gsub("^_+|_+$", "", s)
    s <- substr(s, 1, max_len)
    if (!nzchar(s)) s <- "dx"
    s
  }
  
  gridXY <- grid_from_base(Base, nx = as.integer(BF_NGR_UV), ny = as.integer(BF_NGR_UV))
  dx_surface_files <- c()
  
  for (dx in names(DxW_A)) {
    y <- as.integer(DxW_A[[dx]] > 0)
    pg <- predict_dx_surface(y, Base, gridXY, k_gam = 30)
    if (is.null(pg)) {
      warning(sprintf("[base-field] skip dx=%s (fit failed).", dx))
      next
    }
    out <- data.frame(b1 = gridXY[, 1], b2 = gridXY[, 2], p = pg, dx = dx)
    fp <- file.path(paste0(OUTPUTS_DIR, "/base_prob_grids"), paste0("probgrid_", make_fname_safe(dx), ".csv"))
    write_csv(out, fp)
    dx_surface_files <- c(dx_surface_files, fp)
  }
  cat(sprintf("[base-field] wrote %d grid(s) to %s\n", length(dx_surface_files), paste0(OUTPUTS_DIR, "/base_prob_grids/")))
  
  # --- Item x Diagnosis Interactions (Parallelized OOF GAM) ---

  if (isTRUE(RUN_ITEM_DIAGNOSTICS)) {
  vars_item_interact <- unique(varmap)
  n_probe_count <- if (exists("N_TOP_PER_DX")) N_TOP_PER_DX else 50
  vars_probe <- head(vars_item_interact, min(n_probe_count, length(vars_item_interact)))
  
  # Reuse the filtered, unweighted standardized encoding so varmap and columns stay aligned.
  if (!exists("Z0_std")) Z0_std <- Z0
  
  task_grid <- expand.grid(
    dx = names(DxW_A), 
    var = vars_probe, 
    stringsAsFactors = FALSE
  )
  
  # Filter rare diagnoses for stability
  dx_counts <- colSums(DxW_A > 0, na.rm = TRUE)
  valid_dx <- names(dx_counts)[dx_counts >= 6]
  task_grid <- task_grid[task_grid$dx %in% valid_dx, , drop = FALSE]
  
  cat(sprintf("[Interactions] Parallelizing %d tasks (%d items x %d diagnoses)...\n", 
              nrow(task_grid), length(unique(task_grid$var)), length(unique(task_grid$dx))))
  
  # Helper Helpers
  choose_K_dx <- function(y, K_target = 5L, min_per_class = 6L) {
    y <- as.integer(y > 0)
    n1 <- sum(y == 1); n0 <- sum(y == 0)
    max(2L, min(K_target, floor(n1 / min_per_class), floor(n0 / min_per_class)))
  }
  
  oof_R2_two_gams <- function(v, Base, dx,
                              K_target = 5,
                              k_gam = 10,
                              fold_id = NULL,
                              seed = SEED_GLOBAL,
                              seed_key = NULL) {
    
    n <- length(v)
    if (n != nrow(Base) || n != length(dx)) {
      return(c(R2_add = NA, R2_int = NA, p_like = NA, dR2 = NA))
    }
    
    K <- choose_K_dx(dx, K_target = K_target, min_per_class = 6L)
    if (K < 2) return(c(R2_add = NA, R2_int = NA, p_like = NA, dR2 = NA))
    
    if (is.null(fold_id)) {
      seed_eff <- if (!is.null(seed_key)) .seed_from_key(seed, seed_key) else seed
      ybin <- as.integer(dx > 0)
      fold_id <- make_strat_folds(ybin, K, seed = seed_eff)
    } else {
      stopifnot(length(fold_id) == n)
    }
    
    lev_all <- levels(factor(dx))
    y_add <- rep(NA_real_, n)
    y_int <- rep(NA_real_, n)
    
    ctrl <- mgcv::gam.control(maxit = 100, nthreads = 1)
    
    for (k in sort(unique(fold_id))) {
      tr <- which(fold_id != k)
      te <- which(fold_id == k)
      if (length(tr) < 10 || length(te) == 0) next
      
      dtr <- data.frame(v = v[tr], b1 = Base[tr, 1], b2 = Base[tr, 2], dx = factor(dx[tr], levels = lev_all))
      dte <- data.frame(b1 = Base[te, 1], b2 = Base[te, 2], dx = factor(dx[te], levels = lev_all))
      
      f_add <- try(mgcv::gam(v ~ s(b1, b2, k = k_gam, bs = "tp", m = 2),
                             data = dtr, method = "REML", gamma = 1.4, control = ctrl), silent = TRUE)
      
      f_int <- try(mgcv::gam(v ~ s(b1, b2, k = k_gam, bs = "tp", m = 2) +
                               dx +
                               s(b1, b2, by = dx, k = k_gam, bs = "tp", m = 2),
                             data = dtr, method = "REML", select = TRUE, gamma = 1.4, control = ctrl), silent = TRUE)
      
      mu <- mean(dtr$v, na.rm = TRUE)
      
      if (!inherits(f_add, "try-error")) {
        pa <- try(predict(f_add, newdata = dte, type = "response"), silent = TRUE)
        if (!inherits(pa, "try-error")) y_add[te] <- as.numeric(pa)
      }
      
      if (!inherits(f_int, "try-error")) {
        pi <- try(predict(f_int, newdata = dte, type = "response"), silent = TRUE)
        if (!inherits(pi, "try-error")) y_int[te] <- as.numeric(pi)
      }
      
      y_add[te][!is.finite(y_add[te])] <- mu
      y_int[te][!is.finite(y_int[te])] <- mu
    }
    
    ve <- stats::var(v, na.rm = TRUE)
    if (!is.finite(ve) || ve <= 0) return(c(R2_add = NA, R2_int = NA, p_like = NA, dR2 = NA))
    
    R2_add <- max(0, 1 - mean((v - y_add)^2, na.rm = TRUE) / ve)
    R2_int <- max(0, 1 - mean((v - y_int)^2, na.rm = TRUE) / ve)
    
    d_sq <- (v - y_add)^2 - (v - y_int)^2
    d_sq <- d_sq[is.finite(d_sq)]
    
    p_like <- if (length(d_sq) < 10L || all(abs(d_sq) < .Machine$double.eps)) NA_real_ else
      tryCatch(stats::wilcox.test(d_sq, mu = 0, alternative = "greater", exact = FALSE)$p.value,
               error = function(e) NA_real_)
    
    c(R2_add = R2_add, R2_int = R2_int, p_like = p_like, dR2 = R2_int - R2_add)
  }
  
  # Execute Parallel Loop
  rows <- progressr::with_progress({
    p <- progressr::progressor(steps = nrow(task_grid))
    
    FUTURE_LAPPLY(
      seq_len(nrow(task_grid)),
      function(i) {
        out <- interact_worker(i, task_grid, Z0_std, Base, DxW_A, varmap)
        p() 
        out
      },
      future.globals = list(
        task_grid = task_grid, Z0_std = Z0_std, Base = Base, DxW_A = DxW_A, varmap = varmap,
        score_item_base = score_item_base, oof_R2_two_gams = oof_R2_two_gams,
        choose_K_dx = choose_K_dx, interact_worker = interact_worker, p = p,
        .seed_from_key = .seed_from_key, .with_seed = .with_seed, .hash32 = .hash32,
        make_strat_folds = make_strat_folds, SEED_GLOBAL = SEED_GLOBAL
      ),
      future.packages = c("mgcv", "stats"),
      future.seed = TRUE,
      future.scheduling = 1  # Force 1 item per worker at a time
    )
  })
  
  # Aggregate Results
  rows <- Filter(Negate(is.null), rows)
  
  if (length(rows) > 0) {
    int_tbl <- dplyr::bind_rows(rows) %>% dplyr::filter(is.finite(dR2), is.finite(p_like))
    
    if (nrow(int_tbl) > 0) {
      int_tbl <- int_tbl %>%
        dplyr::group_by(dx) %>%
        dplyr::mutate(q_like = p.adjust(p_like, method = "BH")) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(sig = q_like < SIG_Q) %>%
        dplyr::arrange(dplyr::desc(dR2))
      
      write_csv(int_tbl, "item_dx_interactions.csv")
      cat(sprintf("[interactions] wrote: %s (rows=%d)\n", "item_dx_interactions.csv", nrow(int_tbl)))
    } else {
      warning("[interactions] Table empty after filtering.")
    }
  } else {
    warning("[interactions] No valid interactions computed.")
  }
  } else {
    message("[interactions] Diagnostic-only parent-item PC1 interactions skipped (RUN_ITEM_DIAGNOSTICS=FALSE).")
  }
} else {
  message("[dx] Section 8 (predictive diagnostics & interactions) disabled: DX set as optional (setup.R)")
}
# ==============================================================================
# 9. Outputs and Session
# ==============================================================================
if (DX_AVAILABLE && !is.null(DX_wide) && ncol(DX_wide) > 0) {
  # Cluster Summaries
  kF <- length(unique(clF))
  cat(sprintf("[clusters] using %d clusters.\n", kF))
  
  clusters <- sort(unique(clF))
  cl_fac <- factor(clF, levels = clusters)
  
  mean_by_cluster <- function(vec, clf) {
    out <- tapply(vec, clf, function(z) mean(z, na.rm = TRUE))
    as.numeric(out)
  }
  
  n_by_cluster <- as.integer(tabulate(cl_fac))
  base_means <- cbind(b1_mean = mean_by_cluster(Base[, 1], cl_fac), b2_mean = mean_by_cluster(Base[, 2], cl_fac))
  
  f_means <- NULL
  if (exists("Bprime") && is.matrix(Bprime) && ncol(Bprime) > 0L) {
    f_means <- sapply(seq_len(ncol(Bprime)), function(j) mean_by_cluster(Bprime[, j], cl_fac))
    if (!is.null(dim(f_means))) {
      colnames(f_means) <- paste0("f", seq_len(ncol(Bprime)), "_mean")
    } else {
      f_means <- cbind(`f1_mean` = as.numeric(f_means))
    }
  }
  
  sum_tbl <- data.frame(cluster = clusters, n = n_by_cluster, base_means, check.names = FALSE, stringsAsFactors = FALSE)
  if (!is.null(f_means)) {
    f_means <- as.matrix(f_means)
    storage.mode(f_means) <- "double"
    sum_tbl <- cbind(sum_tbl, f_means, stringsAsFactors = FALSE)
  }
  write_csv(sum_tbl, "cluster_summary.csv")
  
  # Cluster Enrichment Analysis
  DxW2 <- DxW_A
  DxW2[is.na(DxW2)] <- 0L
  prev2 <- colSums(DxW2 > 0L, na.rm = TRUE)
  keepc <- prev2 >= max(DX_CASES_MIN, ceiling(0.01 * nrow(DxW2)))
  DxW2 <- DxW2[, keepc, drop = FALSE]
  
  if (ncol(DxW2) > 0) {
    tab <- sapply(colnames(DxW2), function(dn) {
      tapply(DxW2[[dn]] > 0L, factor(clF, levels = clusters), sum, na.rm = TRUE)
    })
    tab <- as.matrix(tab)
    storage.mode(tab) <- "double"
    rownames(tab) <- as.character(clusters)
    
    suppressWarnings({
      chi <- try(chisq.test(tab), silent = TRUE)
    })
    
    # Always compute expected counts ourselves from tab
    Eexp_used <- outer(
      rowSums(tab),
      colSums(tab),
      function(r, c) r * c / max(sum(tab), 1)
    )
    dimnames(Eexp_used) <- dimnames(tab)
    
    if (inherits(chi, "try-error")) {
      Z_enrich <- (tab - Eexp_used) / sqrt(pmax(Eexp_used, 1e-9))
    } else {
      Z_enrich <- chi$stdres
      # enforce matrix shape even if stdres came back as a vector
      if (is.null(dim(Z_enrich))) {
        Z_enrich <- matrix(
          Z_enrich,
          nrow = nrow(tab),
          ncol = ncol(tab),
          dimnames = dimnames(tab)
        )
      } else {
        dimnames(Z_enrich) <- dimnames(tab)
      }
    }
    
    # from here on, Z_enrich and Eexp_used are guaranteed matrices with same dims as tab
    pmat <- 2 * pnorm(-abs(Z_enrich))
    pmat <- matrix(
      pmat,
      nrow = nrow(tab),
      ncol = ncol(tab),
      dimnames = dimnames(tab)
    )
    
    # BH corrections
    q_by_dx <- apply(pmat, 2, function(p) p.adjust(p, method = "BH"))
    q_global <- p.adjust(as.vector(pmat), method = "BH")
    q_global_mat <- matrix(
      q_global,
      nrow = nrow(pmat), ncol = ncol(pmat),
      byrow = FALSE, dimnames = dimnames(pmat)
    )
    
    # Prevalences
    ncl <- as.integer(tabulate(factor(clF, levels = clusters)))
    names(ncl) <- as.character(clusters)
    prev_in_cluster <- sweep(tab, 1, ncl, "/")
    prev_overall <- colSums(DxW2 > 0L) / nrow(DxW2)
    
    # Export wide matrices
    write_csv(
      data.frame(cluster = rownames(Z_enrich), as.data.frame(Z_enrich), check.names = FALSE),
      "dx_cluster_enrichment_Z.csv"
    )
    write_csv(
      data.frame(cluster = rownames(pmat), as.data.frame(pmat), check.names = FALSE),
      "dx_cluster_enrichment_p.csv"
    )
    write_csv(
      data.frame(cluster = rownames(q_by_dx), as.data.frame(q_by_dx), check.names = FALSE),
      "dx_cluster_enrichment_q_by_dx.csv"
    )
    write_csv(
      data.frame(cluster = rownames(q_global_mat), as.data.frame(q_global_mat), check.names = FALSE),
      "dx_cluster_enrichment_q_global.csv"
    )
    write_csv(
      data.frame(cluster = rownames(prev_in_cluster), as.data.frame(prev_in_cluster), check.names = FALSE),
      "dx_cluster_prev_in_cluster.csv"
    )
    
    # Diagnostics on shapes, right before long table construction
    cat("[diag] dim(tab)         =", paste(dim(tab), collapse = " x "), "\n")
    cat("[diag] dim(Eexp_used)   =", paste(dim(Eexp_used), collapse = " x "), "\n")
    cat("[diag] dim(Z_enrich)    =", paste(dim(Z_enrich), collapse = " x "), "\n")
    cat("[diag] dim(pmat)        =", paste(dim(pmat), collapse = " x "), "\n")
    
    # Export Long table
    grid <- expand.grid(
      cluster = rownames(tab),
      dx = colnames(tab),
      stringsAsFactors = FALSE
    )
    idx_r <- match(grid$cluster, rownames(tab))
    idx_c <- match(grid$dx, colnames(tab))
    
    cat("[diag] nrow(grid)       =", nrow(grid), "\n")
    cat("[diag] length(idx_r)    =", length(idx_r), "\n")
    cat("[diag] length(idx_c)    =", length(idx_c), "\n")
    
    # sanity checks – if these fail, we want to know
    stopifnot(
      length(idx_r) == nrow(grid),
      length(idx_c) == nrow(grid),
      length(Eexp_used) == nrow(tab) * ncol(tab)
    )
    
    grid$n_cluster     <- ncl[grid$cluster]
    grid$cases_cluster <- tab[cbind(idx_r, idx_c)]
    grid$prev_cluster  <- prev_in_cluster[cbind(idx_r, idx_c)]
    grid$prev_overall  <- prev_overall[grid$dx]
    grid$Z             <- Z_enrich[cbind(idx_r, idx_c)]
    grid$p             <- pmat[cbind(idx_r, idx_c)]
    grid$q_by_dx       <- q_by_dx[cbind(idx_r, idx_c)]
    grid$q_global      <- q_global_mat[cbind(idx_r, idx_c)]
    grid$expected_cases <- Eexp_used[cbind(idx_r, idx_c)]
    grid$enriched      <- (grid$Z > 0) & (grid$q_by_dx < SIG_Q)
    
    grid <- grid[order(-grid$Z, grid$q_by_dx, grid$dx, as.integer(grid$cluster)), ]
    write_csv(grid, "dx_cluster_enrichment_long.csv")
    
    if (requireNamespace("ggplot2", quietly = TRUE)) {
      # Cap Z for display
      Zcap <- pmax(pmin(Z_enrich, 5), -5)
      z_long <- data.frame(
        cluster = rep(rownames(Zcap), times = ncol(Zcap)),
        dx = rep(colnames(Zcap), each = nrow(Zcap)),
        Z = as.vector(Zcap),
        stringsAsFactors = FALSE
      )
      gp <- ggplot2::ggplot(z_long, ggplot2::aes(x = dx, y = factor(cluster))) +
        ggplot2::geom_tile(ggplot2::aes(fill = Z)) +
        ggplot2::scale_fill_gradient2() +
        ggplot2::labs(
          x = "Diagnosis", y = "Cluster", fill = "Z (std. resid.)",
          title = "Diagnosis enrichment by cluster (std. residuals)"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
      print(gp)
      save_plot_gg("FIG_dx_cluster_enrichment_heatmap.png",
                   gp,
                   width = 10, height = 6, dpi = 150
      )
    }
  } else {
    warning("[clusters] No Dx columns after filtering for enrichment; skipped.")
  }
} else {
  message("[dx] Enrichment by diagnosis disabled: DX set as optional.")
}
# ------------------------------------------------------------------------------
# Encoding & geometry
# ------------------------------------------------------------------------------

as_base2_numeric <- function(Base_in) {
  A <- as.data.frame(Base_in[, 1:2, drop = FALSE])
  for (j in 1:2) {
    v <- suppressWarnings(as.numeric(A[[j]]))
    m <- mean(v, na.rm = TRUE)
    if (!is.finite(m)) m <- 0
    A[[j]] <- ifelse(is.finite(v), v, m)
  }
  B <- as.matrix(A)
  storage.mode(B) <- "double"
  colnames(B) <- c("b1", "b2")
  rownames(B) <- rownames(Base_in)
  B
}

# whiten base to unit disk (pure linear algebra)
standardise_to_circle <- function(Base2, cover = 0.95) {
  X <- as.matrix(Base2[, 1:2, drop = FALSE]); storage.mode(X) <- "double"
  mu <- colMeans(X); S <- stats::cov(X)
  if (!all(is.finite(S)) || det(S) <= 0) S <- S + diag(1e-8, 2)
  
  eig <- eigen(S, symmetric = TRUE)
  V   <- eig$vectors
  L   <- pmax(eig$values, 1e-12)
  S_half <- V %*% diag(sqrt(L), 2) %*% t(V)
  S_hi   <- V %*% diag(1/sqrt(L), 2) %*% t(V)
  
  Xc <- sweep(X, 2, mu, "-")
  U0 <- Xc %*% t(S_hi)
  r  <- sqrt(rowSums(U0^2))
  s  <- as.numeric(stats::quantile(r, probs = cover, na.rm = TRUE))
  if (!is.finite(s) || s <= 0) s <- max(r, na.rm = TRUE)
  
  fwd <- function(xb1, xb2) {
    xb1 <- as.numeric(xb1); xb2 <- as.numeric(xb2)
    Xnew <- cbind(xb1, xb2)
    Uq   <- (sweep(Xnew, 2, mu, "-") %*% t(S_hi)) / s
    colnames(Uq) <- c("u1", "u2"); Uq
  }
  inv <- function(u1, u2) {
    Unew <- cbind(as.numeric(u1), as.numeric(u2)) * s
    Xq   <- Unew %*% t(S_half) + matrix(mu, nrow(Unew), 2, byrow = TRUE)
    colnames(Xq) <- c("b1", "b2"); Xq
  }
  xstd_env <- list2env(
    list(mu = mu, S_half = S_half, S_hi = S_hi, s = s),
    parent = baseenv()
  )
  environment(fwd) <- xstd_env
  environment(inv) <- xstd_env
  
  list(mu = mu, S_half = S_half, S_half_inv = S_hi, s = s, fwd = fwd, inv = inv)
}

standardise_to_ball <- function(Base3, cover = 0.95, eps = 1e-8) {
  X <- as.matrix(Base3[, 1:3, drop = FALSE])
  storage.mode(X) <- "double"
  
  mu <- colMeans(X, na.rm = TRUE)
  Xc <- sweep(X, 2, mu, "-")
  
  S <- stats::cov(Xc, use = "pairwise.complete.obs")
  U <- chol(S + diag(eps, 3))              # upper-triangular; Xc %*% solve(U) whitens
  
  Y <- Xc %*% solve(U)
  r <- sqrt(rowSums(Y^2))
  s <- as.numeric(stats::quantile(r, probs = cover, na.rm = TRUE))
  if (!is.finite(s) || s <= 0) s <- max(r, na.rm = TRUE)
  if (!is.finite(s) || s <= 0) s <- 1
  
  inv <- function(u1, u2, u3) {
    u1 <- as.numeric(u1); u2 <- as.numeric(u2); u3 <- as.numeric(u3)
    n  <- max(length(u1), length(u2), length(u3))
    if (n == 0L) {
      out <- matrix(numeric(0), 0, 3)
      colnames(out) <- c("b1","b2","b3")
      return(out)
    }
    Unew <- cbind(rep_len(u1, n), rep_len(u2, n), rep_len(u3, n)) * s
    Xq   <- Unew %*% U + matrix(mu, nrow(Unew), 3, byrow = TRUE)
    colnames(Xq) <- c("b1","b2","b3")
    Xq
  }
  
  fwd <- function(xb1, xb2, xb3) {
    xb1 <- as.numeric(xb1); xb2 <- as.numeric(xb2); xb3 <- as.numeric(xb3)
    n   <- max(length(xb1), length(xb2), length(xb3))
    if (n == 0L) {
      out <- matrix(numeric(0), 0, 3)
      colnames(out) <- c("u1","u2","u3")
      return(out)
    }
    Xnew <- cbind(rep_len(xb1, n), rep_len(xb2, n), rep_len(xb3, n))
    Uq   <- (sweep(Xnew, 2, mu, "-") %*% solve(U)) / s
    colnames(Uq) <- c("u1","u2","u3")
    Uq
  }
  xstd_env <- list2env(
    list(mu = mu, U = U, s = s),
    parent = baseenv()
  )
  environment(fwd) <- xstd_env
  environment(inv) <- xstd_env
  
  list(mu = mu, U = U, s = s, fwd = fwd, inv = inv)
}

make_unitdisk_square <- function(nu) {
  u1 <- seq(-1, 1, length.out = nu)
  u2 <- seq(-1, 1, length.out = nu)
  g <- expand.grid(u1 = u1, u2 = u2)
  mask <- with(g, sqrt(u1^2 + u2^2) <= 1 + 1e-9)
  list(u1 = u1, u2 = u2, grid = g, mask_vec = mask)
}

make_unitball_cube <- function(nu) {
  u1 <- seq(-1, 1, length.out = nu)
  u2 <- seq(-1, 1, length.out = nu)
  u3 <- seq(-1, 1, length.out = nu)
  g  <- expand.grid(u1 = u1, u2 = u2, u3 = u3)
  mask <- with(g, sqrt(u1^2 + u2^2 + u3^2) <= 1 + 1e-9)
  list(u1 = u1, u2 = u2, u3 = u3, grid = g, mask_vec = mask)
}

inside_hull <- function(px, py, poly, eps = 1e-12) {
  n <- nrow(poly)
  j <- n
  inside <- rep(FALSE, length(px))
  on_edge <- rep(FALSE, length(px))
  for (i in seq_len(n)) {
    xi <- poly$b1[i]; yi <- poly$b2[i]
    xj <- poly$b1[j]; yj <- poly$b2[j]
    dx <- xj - xi; dy <- yj - yi
    seg_len2 <- dx * dx + dy * dy + eps
    t <- pmin(1, pmax(0, ((px - xi) * dx + (py - yi) * dy) / seg_len2))
    projx <- xi + t * dx; projy <- yi + t * dy
    on_edge <- on_edge | ((px - projx)^2 + (py - projy)^2 <= (1e-9)^2)
    cross <- ((yi > py) != (yj > py)) & (px < (xj - xi) * (py - yi) / (yj - yi + eps) + xi)
    inside <- xor(inside, cross)
    j <- i
  }
  inside | on_edge
}

build_geometry <- function(Base_A, GRID_N_B = 300, GRID_N_U = 300) {
  d <- ncol(Base_A)
  if (d < 2) stop("Base_A must have at least 2 columns")
  
  if (d >= 3) {
    B1 <- as.numeric(Base_A[,1]); B2 <- as.numeric(Base_A[,2]); B3 <- as.numeric(Base_A[,3])
    df_base <- data.frame(b1 = B1, b2 = B2, b3 = B3)
    rownames(df_base) <- rownames(Base_A)
    
    Xstd <- standardise_to_ball(Base_A[,1:3, drop = FALSE], cover = 0.95)
    U <- as.data.frame(Xstd$fwd(B1, B2, B3))
    rownames(U) <- rownames(Base_A)
    
    U$rin  <- sqrt(U$u1^2 + U$u2^2 + U$u3^2) <= 1 + 1e-9
    U$insq <- abs(U$u1) <= 1 + 1e-9 & abs(U$u2) <= 1 + 1e-9 & abs(U$u3) <= 1 + 1e-9
    
    # A full cube scales as GRID_N_B^3. With GRID_N_B=400 that is 64M
    # points, and contplot only needs 2D plane slices built later.
    return(list(
      Xstd = Xstd, U = U,
      df_base = df_base
    ))
  }
  
  # 2D path (your existing code), unchanged:
  B1 <- Base_A[, 1]; B2 <- Base_A[, 2]
  df_base <- data.frame(b1 = as.numeric(B1), b2 = as.numeric(B2))
  rownames(df_base) <- rownames(Base_A)
  
  Xstd <- standardise_to_circle(Base_A, cover = 0.95)
  U <- as.data.frame(Xstd$fwd(B1, B2))
  names(U) <- c("u1", "u2")
  rownames(U) <- rownames(Base_A)
  
  UD <- make_unitdisk_square(GRID_N_B)
  gridX_sq <- as.data.frame(Xstd$inv(UD$grid$u1, UD$grid$u2))
  mask_sq <- UD$mask_vec
  
  hidx <- grDevices::chull(df_base$b1, df_base$b2)
  hpoly <- df_base[hidx, , drop = FALSE]
  
  qx <- range(df_base$b1); qy <- range(df_base$b2)
  gridB_full <- expand.grid(
    b1 = seq(qx[1], qx[2], length.out = GRID_N_U),
    b2 = seq(qy[1], qy[2], length.out = GRID_N_U)
  )
  mask_hull <- inside_hull(gridB_full$b1, gridB_full$b2, hpoly)
  
  U$rin  <- sqrt(U$u1^2 + U$u2^2) <= 1 + 1e-9
  U$insq <- abs(U$u1) <= 1 + 1e-9 & abs(U$u2) <= 1 + 1e-9
  
  list(
    Xstd = Xstd, U = U, UD = UD,
    gridX_sq = gridX_sq, mask_sq = mask_sq,
    df_base = df_base, gridB_full = gridB_full, mask_hull = mask_hull
  )
}

# Basic objects
Z_A    <- Z

# Geometry & grids
geom <- build_geometry(Base_A, GRID_N_B = GRID_N_B, GRID_N_U = GRID_N_U)

if (ncol(Base_A) >= 3) {
  geom$Xstd3 <- standardise_to_ball(Base_A, cover = 0.95)
  U3 <- geom$Xstd3$fwd(Base_A[,1], Base_A[,2], Base_A[,3])
  geom$U3 <- as.data.frame(U3)
  rownames(geom$U3) <- rownames(Base_A)
}

saveRDS(
  list(
    Base_A = Base_A,
    geom = geom,
    fold_id = if (exists("fold_id", inherits = FALSE)) fold_id else NULL
  ),
  file = "contplot_state.rds"
)
message("[export] contplot_state.rds written for clean-session contplot runs.")

# ================= 9a) EXPORT EMBEDDINGS (robust, DX-agnostic) =================

# Ensure Base has rownames we can trust
if (is.null(rownames(Base)) || length(rownames(Base)) != nrow(Base)) {
  if (exists("Xenc_w") && !is.null(rownames(Xenc_w)) &&
      length(rownames(Xenc_w)) == nrow(Base)) {
    rownames(Base) <- rownames(Xenc_w)
  } else {
    rownames(Base) <- make.unique(rep("row", nrow(Base)))
  }
}

# Name cluster membership by the object it was computed on (usually Bprime)
if (exists("clF")) {
  rn_bprime <- if (exists("Bprime") && !is.null(rownames(Bprime)) &&
                   length(rownames(Bprime)) == length(clF)) rownames(Bprime) else rownames(Base)
  cl_named <- setNames(as.integer(clF), rn_bprime)
} else {
  warning("[export] cluster membership 'clF' not found; exporting clusters as NA.")
  cl_named <- setNames(rep(NA_integer_, nrow(Base)), rownames(Base))
}

# Align IDs between Base and clusters
id_base   <- rownames(Base)
id_common <- intersect(id_base, names(cl_named))
if (!length(id_common)) stop("[export] No overlapping IDs between Base and clusters.")

# Put everything in Base order
o <- match(id_base, id_common, nomatch = 0L)
keep <- which(o > 0L)
id_use <- id_base[keep]

Base_use    <- Base[id_use, , drop = FALSE]
cluster_use <- unname(cl_named[id_use])

# Write emb_base to OUTPUTS_DIR (this is where the longitudinal script reads it from)
dir.create(OUTPUTS_DIR, recursive = TRUE, showWarnings = FALSE)
emb_base <- data.frame(
  participant_id = id_use,
  b1 = Base_use[, 1],
  b2 = Base_use[, 2],
  cluster = as.integer(cluster_use),
  stringsAsFactors = FALSE
)
readr::write_csv2(emb_base, file.path(OUTPUTS_DIR, "embedding_base_b1b2.csv"))

if (BASE_DIM >= 3L && ncol(Base_use) >= 3L) {
  emb_base3 <- data.frame(
    participant_id = id_use,
    b1 = Base_use[, 1],
    b2 = Base_use[, 2],
    b3 = Base_use[, 3],
    cluster = as.integer(cluster_use),
    stringsAsFactors = FALSE
  )
  readr::write_csv2(emb_base3, file.path(OUTPUTS_DIR, "embedding_base_b1b2b3.csv"))
  message(sprintf(
    "[export] embedding_base_b1b2b3.csv rows=%d (BASE_DIM=%d)",
    nrow(emb_base3), BASE_DIM
  ))
}

# Refresh ids_base for downstream code to the aligned IDs we actually exported
ids_base <- id_use

# Export residual embedding (align rows by name; fall back to intersection)
if (exists("Bprime")) {
  if (is.null(rownames(Bprime)) || length(rownames(Bprime)) != nrow(Bprime)) {
    rownames(Bprime) <- ids_base[seq_len(min(length(ids_base), nrow(Bprime)))]
  }
  id_res <- intersect(ids_base, rownames(Bprime))
  Bprime_use <- Bprime[id_res, , drop = FALSE]
  cl_res_use <- unname(cl_named[id_res])
  emb_residual <- data.frame(
    participant_id = id_res,
    Bprime_use,
    cluster = as.integer(cl_res_use),
    stringsAsFactors = FALSE
  )
  readr::write_csv2(emb_residual, file.path(OUTPUTS_DIR, "embedding_residual_Bprime.csv"))
}

# Export XR/Fprime matrix with aligned participant_id
ids_for_mats <- ids_base
XR_use <- if (exists("E_scaled") && !is.null(rownames(E_scaled))) {
  ix <- match(ids_for_mats, rownames(E_scaled), nomatch = 0L); E_scaled[ix[ix > 0], , drop = FALSE]
} else if (exists("E_scaled")) {
  E_scaled
} else {
  NULL
}
Fprime_use <- if (exists("Fprime") && !is.null(rownames(Fprime))) {
  ix <- match(ids_for_mats, rownames(Fprime), nomatch = 0L); Fprime[ix[ix > 0], , drop = FALSE]
} else if (exists("Fprime")) {
  Fprime
} else {
  NULL
}

saveRDS(list(
  participant_id = ids_for_mats,
  XR = XR_use,
  Fprime = Fprime_use
), file = file.path(OUTPUTS_DIR, "Fprime_matrix.rds"))

nrow_or_na <- function(x) if (is.null(x)) "NA" else as.character(nrow(x))
message(sprintf("[export] embedding_base_b1b2.csv rows=%d; residual_export=%s; XR rows=%s; Fprime rows=%s",
                nrow(emb_base),
                if (exists("emb_residual")) nrow(emb_residual) else "NA",
                nrow_or_na(XR_use),
                nrow_or_na(Fprime_use)))

# Export Weights and Maps
w_tbl <- data.frame(
  var = names(w_full),
  weight = as.numeric(w_full),
  selected = names(w_full) %in% survivors,
  stringsAsFactors = FALSE
)
write_csv(w_tbl, "gower_weights_optimised.csv")

enc_map <- data.frame(
  mm_col = colnames(Xenc_w),
  source_var = as.character(varmap),
  response_level = as.character(levelmap),
  weight_share = as.numeric(w_enc[colnames(Xenc_w)]),
  stringsAsFactors = FALSE
)
write_csv(enc_map, "encoding_map_and_weight_share.csv")

# Session Info
setup_snapshot <- list(
  outputs_dir = OUTPUTS_DIR,
  base_dim = BASE_DIM,
  base_decomp_method = BASE_DECOMP_METHOD,
  weighting_mode = WEIGHTING_MODE,
  treat_ordinals_as_nominal = TREAT_ORDINALS_AS_NOMINAL,
  item_response_schema_csv = schema_path,
  invalid_response_action = INVALID_RESPONSE_ACTION,
  max_participant_missing_prop = MAX_PARTICIPANT_MISSING_PROP,
  pca_nominal_missing = PCA_NOMINAL_MISSING,
  missing_as_nominal_level = MISSING_AS_NOMINAL_LEVEL,
  run_item_diagnostics = RUN_ITEM_DIAGNOSTICS,
  n_rows_sub = N_ROWS_SUB,
  final_diag_mode = FINAL_DIAG_MODE,
  diag_n_max = DIAG_N_MAX,
  dedup_mode = DEDUP_MODE,
  dedup_hash_digits = DEDUP_HASH_DIGITS,
  gower_multi_enable = GOWER_MULTI_ENABLE,
  gower_multi_runs = GOWER_MULTI_RUNS,
  grid_n_b = GRID_N_B,
  grid_n_u = GRID_N_U,
  palette_engine = PALETTE_ENGINE,
  palette_name = PALETTE_NAME
)
capture.output(
  {
    print(sessionInfo())
    cat("\n\nSetup snapshot:\n")
    str(setup_snapshot, max.level = 1)
  },
  file = "sessionInfo.txt"
)

cat("[export] wrote embeddings, enrichment tables, weights, and session info.\n")

# Check alignment
cat("[Orientation] Base head:\n")
print(head(Base[, 1:2]))
cat("[Orientation] Base_A head:\n")
print(head(Base_A))

# ------------------------------------------------------------------------------
# Base-space Plots
# ------------------------------------------------------------------------------

density_plots <- function(geom) {
  stopifnot(all(c("U") %in% names(geom)))
  U_df <- geom$U
  if ("insq" %in% names(U_df)) U_df <- subset(U_df, insq)
  
  ggplot() +
    stat_density_2d(
      data = U_df,
      aes(u1, u2, colour = after_stat(level)),
      bins = 6, linewidth = 0.45, alpha = 0.9
    ) +
    { scico::scale_colour_scico(palette = PALETTE_NAME, direction = PALETTE_DIRECTION) } +
    geom_point(
      data = U_df, aes(u1, u2),
      shape = 16, size = 0.9, colour = scales::alpha("black", 0.35)
    ) +
    annotate("rect", xmin = -1, xmax = 1, ymin = -1, ymax = 1,
             fill = NA, colour = "black", linewidth = 0.35) +
    coord_equal(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
    labs(x = "u1 (whitened b1,b2)", y = "u2") +
    theme_pub(12)
}

density_plots_base <- function(Base_A, cover = 1, pad_frac = 0.10) {
  dfB <- data.frame(
    b1 = as.numeric(Base_A[, 1]),
    b2 = as.numeric(Base_A[, 2]),
    stringsAsFactors = FALSE
  )
  dfB <- dfB[is.finite(dfB$b1) & is.finite(dfB$b2), , drop = FALSE]
  if (!nrow(dfB)) return(NULL)
  
  rx <- range(dfB$b1, na.rm = TRUE)
  ry <- range(dfB$b2, na.rm = TRUE)
  wx <- diff(rx); wy <- diff(ry)
  if (!is.finite(wx) || wx <= 0) wx <- 1
  if (!is.finite(wy) || wy <= 0) wy <- 1
  
  xlim <- c(rx[1] - pad_frac * wx, rx[2] + pad_frac * wx)
  ylim <- c(ry[1] - pad_frac * wy, ry[2] + pad_frac * wy)
  
  ggplot2::ggplot() +
    ggplot2::stat_density_2d(
      data = dfB,
      ggplot2::aes(b1, b2, colour = after_stat(level)),
      bins = 6, linewidth = 0.45, alpha = 0.9
    ) +
    { scico::scale_colour_scico(palette = PALETTE_NAME, direction = PALETTE_DIRECTION) } +
    ggplot2::geom_point(
      data = dfB, ggplot2::aes(b1, b2),
      shape = 16, size = 0.9, colour = scales::alpha("black", 0.35)
    ) +
    ggplot2::scale_x_continuous(limits = xlim, expand = ggplot2::expansion(mult = 0)) +
    ggplot2::scale_y_continuous(limits = ylim, expand = ggplot2::expansion(mult = 0)) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(x = "b1", y = "b2") +
    theme_pub(12)
}

density_plots_base_pair <- function(Base, i = 1L, j = 2L, name = NULL, pad_frac = 0.10) {
  stopifnot(ncol(Base) >= max(i, j))
  dfB <- data.frame(
    x = as.numeric(Base[, i]),
    y = as.numeric(Base[, j]),
    stringsAsFactors = FALSE
  )
  dfB <- dfB[is.finite(dfB$x) & is.finite(dfB$y), , drop = FALSE]
  if (!nrow(dfB)) return(NULL)
  
  rx <- range(dfB$x, na.rm = TRUE)
  ry <- range(dfB$y, na.rm = TRUE)
  wx <- diff(rx); wy <- diff(ry)
  if (!is.finite(wx) || wx <= 0) wx <- 1
  if (!is.finite(wy) || wy <= 0) wy <- 1
  
  xlim <- c(rx[1] - pad_frac * wx, rx[2] + pad_frac * wx)
  ylim <- c(ry[1] - pad_frac * wy, ry[2] + pad_frac * wy)
  
  ggplot2::ggplot() +
    ggplot2::stat_density_2d(
      data = dfB,
      ggplot2::aes(x, y, colour = after_stat(level)),
      bins = 6, linewidth = 0.45, alpha = 0.9
    ) +
    { scico::scale_colour_scico(palette = PALETTE_NAME, direction = PALETTE_DIRECTION) } +
    ggplot2::geom_point(
      data = dfB,
      ggplot2::aes(x, y),
      shape = 16, size = 0.9,
      colour = scales::alpha("black", 0.35)
    ) +
    ggplot2::scale_x_continuous(limits = xlim, expand = ggplot2::expansion(mult = 0)) +
    ggplot2::scale_y_continuous(limits = ylim, expand = ggplot2::expansion(mult = 0)) +
    ggplot2::coord_equal(expand = FALSE) +
    ggplot2::labs(
      x = paste0("b", i),
      y = paste0("b", j),
      title = name
    ) +
    theme_pub(12)
}

plot_density_fog_3d <- function(X3,
                                cols = 1:3,
                                trim_q = c(0.01, 0.99),
                                ng = 70,
                                cut = 0.02,
                                surface_count = 10,
                                opacity = 0.18,
                                bandwidth = c("Hpi", "Hscv"),
                                h_mult = 1.8,              # << add: inflate bandwidth
                                isomin_q = 0.93,           # << add: render only top tail
                                trace = c("isosurface","volume"),
                                axis_titles = c("dim1", "dim2", "dim3"),
                                show_points = TRUE,
                                points_n = 1500,
                                points_size = 4) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Need package: plotly")
  if (!requireNamespace("ks", quietly = TRUE)) stop("Need package: ks")
  
  trace <- match.arg(trace)
  
  X <- as.matrix(X3)
  storage.mode(X) <- "double"
  X <- X[, cols, drop = FALSE]
  
  okf <- apply(X, 1, function(r) all(is.finite(r)))
  X <- X[okf, , drop = FALSE]
  if (nrow(X) < 10) stop("Not enough finite rows for a stable 3D KDE.")
  
  lims <- apply(X, 2, stats::quantile, probs = trim_q, na.rm = TRUE)
  keep <- X[,1] >= lims[1,1] & X[,1] <= lims[2,1] &
    X[,2] >= lims[1,2] & X[,2] <= lims[2,2] &
    X[,3] >= lims[1,3] & X[,3] <= lims[2,3]
  Xt <- X[keep, , drop = FALSE]
  if (nrow(Xt) < 10) Xt <- X
  
  # bandwidth
  if (is.matrix(bandwidth)) {
    H <- bandwidth
    if (!all(dim(H) == c(3,3))) stop("If bandwidth is a matrix, it must be 3x3.")
  } else {
    bw <- match.arg(bandwidth)
    H <- if (bw == "Hscv") ks::Hscv(Xt) else ks::Hpi(Xt)
  }
  H <- H * (h_mult^2)  # << key: stronger smoothing
  
  # KDE on a grid
  k3 <- ks::kde(
    x        = Xt,
    H        = H,
    binned   = TRUE,
    gridsize = rep(ng, 3),
    xmin     = lims[1, ],
    xmax     = lims[2, ]
  )
  
  gx <- k3$eval.points[[1]]
  gy <- k3$eval.points[[2]]
  gz <- k3$eval.points[[3]]
  D  <- k3$estimate
  
  mx <- max(D, na.rm = TRUE)
  if (!is.finite(mx) || mx <= 0) stop("KDE returned non-finite/degenerate density.")
  D <- D / mx
  
  G <- expand.grid(x = gx, y = gy, z = gz)
  G$val <- as.vector(D)
  
  vf <- G$val[is.finite(G$val)]
  isomin <- max(cut, as.numeric(stats::quantile(vf, probs = isomin_q, na.rm = TRUE)))
  isomax <- 1
  
  p <- plotly::plot_ly(
    data  = G,
    x     = ~x, y = ~y, z = ~z,
    value = ~val,
    type  = trace,
    isomin = isomin,
    isomax = isomax,
    opacity = opacity,
    surface = list(count = surface_count),
    caps = list(x = list(show = FALSE), y = list(show = FALSE), z = list(show = FALSE))
  ) %>%
    plotly::layout(
      scene = list(
        xaxis = list(title = axis_titles[1]),
        yaxis = list(title = axis_titles[2]),
        zaxis = list(title = axis_titles[3])
      )
    )
  
  if (isTRUE(show_points)) {
    n0 <- nrow(X)
    take <- if (n0 <= points_n) seq_len(n0) else .with_seed(SEED_GLOBAL, sample.int(n0, points_n))
    P <- as.data.frame(X[take, , drop = FALSE])
    names(P) <- c("x", "y", "z")
    p <- p %>% plotly::add_markers(
      data = P,
      x = ~x, y = ~y, z = ~z,
      marker = list(size = points_size, opacity = 0.35),
      inherit = FALSE,
      showlegend = FALSE
    )
  }
  
  p
}

direction_wheel_plot <- function(geom) {
  U_sq <- subset(geom$U, insq)
  if (!nrow(U_sq)) return(NULL)
  
  nu <- 500
  pad <- 0.75
  gx <- seq(-1 - pad, 1 + pad, length.out = nu)
  gy <- seq(-1 - pad, 1 + pad, length.out = nu)
  
  kd <- with(U_sq, MASS::kde2d(u1, u2, n = nu,
                               lims = c(-1 - pad, 1 + pad, -1 - pad, 1 + pad)))
  D <- kd$z
  D <- log1p(D / max(D, na.rm = TRUE))
  D <- D / quantile(D, 0.99, na.rm = TRUE)
  D[D > 1] <- 1; D[D < 0] <- 0
  ALPHA <- D^0.70
  
  G <- expand.grid(u1 = gx, u2 = gy)
  theta <- atan2(G$u2, G$u1)
  r <- sqrt(G$u1^2 + G$u2^2)
  
  H0 <- 170; L0 <- 60; Cmax <- 110; betaC <- 0.90
  H <- (H0 + theta * 180 / pi) %% 360
  C <- pmin(Cmax * (pmin(r, 1)^betaC), Cmax)
  L <- pmax(0, pmin(100, L0 - 6 * (pmin(r, 1)^1.1)))
  
  G$fill <- grDevices::hcl(H, C, L)
  G$alpha <- as.vector(ALPHA)
  
  feather1d <- function(x, lo = -1, hi = 1, w = 0.08) {
    tL <- pmin(pmax((x - lo) / w, 0), 1)
    tR <- pmin(pmax((hi - x) / w, 0), 1)
    fL <- (cos(tL * pi / 2))^2
    fR <- (cos(tR * pi / 2))^2
    pmin(fL, fR)
  }
  Fx <- feather1d(G$u1); Fy <- feather1d(G$u2)
  G$alpha <- G$alpha * Fx * Fy
  
  anchor <- data.frame(
    x = c(0.85, -0.85, 0.02, 0.02),
    y = c(0.02, 0.02, 0.85, -0.85),
    lab = c("+b1", "-b1", "+b2", "-b2"),
    col = grDevices::hcl((H0 + c(0, 180, 90, -90)) %% 360, Cmax, L0)
  )
  
  ggplot() +
    geom_raster(data = G, aes(u1, u2, fill = I(fill), alpha = alpha), interpolate = TRUE) +
    scale_alpha(range = c(0, 1), guide = "none") +
    geom_point(data = U_sq, aes(u1, u2),
               shape = 16, size = 0.8, colour = scales::alpha("black", 0.32)) +
    coord_equal(xlim = c(-1.5, 1.5), ylim = c(-1.5, 1.5),
                expand = FALSE, clip = "on") +
    geom_point(data = anchor, aes(x, y), shape = 15, size = 3, colour = anchor$col) +
    geom_text(data = anchor, aes(x, y, label = lab), nudge_x = 0.07, size = 3.2) +
    labs(x = "u1 (whitened b1,b2)", y = "u2") +
    theme_pub(12)
}

direction_wheel_plot_base <- function(Base_A,
                                      cover = 1,
                                      pad_frac = 0.18,
                                      bandwidth_mult = 1.20,
                                      alpha_power = 0.58,
                                      edge_fade_frac = 0.10) {
  dfB <- data.frame(
    b1 = as.numeric(Base_A[, 1]),
    b2 = as.numeric(Base_A[, 2]),
    stringsAsFactors = FALSE
  )
  dfB <- dfB[is.finite(dfB$b1) & is.finite(dfB$b2), , drop = FALSE]
  if (!nrow(dfB)) return(NULL)
  
  Xstd <- standardise_to_circle(Base_A, cover = cover)
  mu <- Xstd$mu
  S_hi <- Xstd$S_half_inv
  s <- Xstd$s
  
  nu <- 500
  
  rx <- range(dfB$b1, na.rm = TRUE)
  ry <- range(dfB$b2, na.rm = TRUE)
  wx <- diff(rx); wy <- diff(ry)
  if (!is.finite(wx) || wx <= 0) wx <- 1
  if (!is.finite(wy) || wy <= 0) wy <- 1
  
  lims <- c(
    rx[1] - pad_frac * wx, rx[2] + pad_frac * wx,
    ry[1] - pad_frac * wy, ry[2] + pad_frac * wy
  )
  
  h <- c(MASS::bandwidth.nrd(dfB$b1), MASS::bandwidth.nrd(dfB$b2)) *
    bandwidth_mult
  kd <- with(dfB, MASS::kde2d(b1, b2, n = nu, lims = lims, h = h))
  D <- kd$z
  D <- log1p(D / max(D, na.rm = TRUE))
  D <- D / stats::quantile(D, 0.99, na.rm = TRUE)
  D[D > 1] <- 1; D[D < 0] <- 0
  ALPHA <- D^alpha_power
  
  gx <- kd$x
  gy <- kd$y
  G <- expand.grid(b1 = gx, b2 = gy)
  
  Uq <- (sweep(as.matrix(G), 2, mu, "-") %*% t(S_hi)) / s
  theta <- atan2(Uq[, 2], Uq[, 1])
  r <- sqrt(Uq[, 1]^2 + Uq[, 2]^2)
  
  H0 <- 170; L0 <- 60; Cmax <- 110; betaC <- 0.90
  H <- (H0 + theta * 180 / pi) %% 360
  C <- pmin(Cmax * (pmin(r, 1)^betaC), Cmax)
  L <- pmax(0, pmin(100, L0 - 6 * (pmin(r, 1)^1.1)))
  
  G$fill <- grDevices::hcl(H, C, L)
  G$alpha <- as.vector(ALPHA)

  edge_fade <- function(x, lo, hi, frac) {
    width <- max((hi - lo) * frac, .Machine$double.eps)
    z <- pmin(pmax(pmin(x - lo, hi - x) / width, 0), 1)
    sin(z * pi / 2)^2
  }
  G$alpha <- G$alpha *
    edge_fade(G$b1, lims[1], lims[2], edge_fade_frac) *
    edge_fade(G$b2, lims[3], lims[4], edge_fade_frac)
  
  ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = G,
      ggplot2::aes(b1, b2, fill = I(fill), alpha = alpha),
      interpolate = TRUE
    ) +
    ggplot2::scale_alpha(range = c(0, 1), guide = "none") +
    ggplot2::geom_point(
      data = dfB,
      ggplot2::aes(b1, b2),
      shape = 16, size = 0.8, colour = scales::alpha("black", 0.32)
    ) +
    ggplot2::coord_equal(
      xlim = lims[1:2], ylim = lims[3:4],
      expand = FALSE, clip = "on"
    ) +
    ggplot2::labs(x = "b1", y = "b2") +
    theme_pub(12)
}

direction_wheel_plot_labeled <- function(geom,
                                         label_col = "participant_id",
                                         label_size = 0.75,
                                         label_alpha = 0.5,
                                         repel = TRUE) {
  U_sq <- subset(geom$U, insq)
  if (!nrow(U_sq)) return(NULL)
  
  # garantir labels
  if (is.null(rownames(U_sq))) {
    stop("[direction_wheel_plot_labeled] rownames(U) must contain participant_id")
  }
  U_sq$label <- rownames(U_sq)
  
  nu <- 500
  pad <- 0.75
  gx <- seq(-1 - pad, 1 + pad, length.out = nu)
  gy <- seq(-1 - pad, 1 + pad, length.out = nu)
  
  kd <- with(U_sq, MASS::kde2d(u1, u2, n = nu,
                               lims = c(-1 - pad, 1 + pad,
                                        -1 - pad, 1 + pad)))
  D <- kd$z
  D <- log1p(D / max(D, na.rm = TRUE))
  D <- D / quantile(D, 0.99, na.rm = TRUE)
  D[D > 1] <- 1
  D[D < 0] <- 0
  ALPHA <- D^0.70
  
  G <- expand.grid(u1 = gx, u2 = gy)
  theta <- atan2(G$u2, G$u1)
  r <- sqrt(G$u1^2 + G$u2^2)
  
  H0 <- 170; L0 <- 60; Cmax <- 110; betaC <- 0.90
  H <- (H0 + theta * 180 / pi) %% 360
  C <- pmin(Cmax * (pmin(r, 1)^betaC), Cmax)
  L <- pmax(0, pmin(100, L0 - 6 * (pmin(r, 1)^1.1)))
  
  G$fill <- grDevices::hcl(H, C, L)
  G$alpha <- as.vector(ALPHA)
  
  feather1d <- function(x, lo = -1, hi = 1, w = 0.08) {
    tL <- pmin(pmax((x - lo) / w, 0), 1)
    tR <- pmin(pmax((hi - x) / w, 0), 1)
    fL <- (cos(tL * pi / 2))^2
    fR <- (cos(tR * pi / 2))^2
    pmin(fL, fR)
  }
  G$alpha <- G$alpha * feather1d(G$u1) * feather1d(G$u2)
  
  anchor <- data.frame(
    x = c(0.85, -0.85, 0.02, 0.02),
    y = c(0.02, 0.02, 0.85, -0.85),
    lab = c("+b1", "-b1", "+b2", "-b2"),
    col = grDevices::hcl((H0 + c(0, 180, 90, -90)) %% 360, Cmax, L0)
  )
  
  p <- ggplot() +
    geom_raster(
      data = G,
      aes(u1, u2, fill = I(fill), alpha = alpha),
      interpolate = TRUE
    ) +
    scale_alpha(range = c(0, 1), guide = "none") +
    geom_point(
      data = U_sq,
      aes(u1, u2),
      shape = 16, size = 0.7,
      colour = scales::alpha("black", 0.35)
    )
  
  if (repel && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p +
      ggrepel::geom_text_repel(
        data = U_sq,
        aes(u1, u2, label = label),
        size = label_size,
        alpha = label_alpha,
        colour = "black",
        max.overlaps = Inf,
        min.segment.length = 0,
        box.padding = 0.15,
        point.padding = 0.05,
        segment.size = 0.15
      )
  } else {
    p <- p +
      geom_text(
        data = U_sq,
        aes(u1, u2, label = label),
        size = label_size,
        alpha = label_alpha,
        colour = "black",
        vjust = -0.3
      )
  }
  
  p +
    coord_equal(
      xlim = c(-1.5, 1.5),
      ylim = c(-1.5, 1.5),
      expand = FALSE,
      clip = "on"
    ) +
    geom_point(
      data = anchor,
      aes(x, y),
      shape = 15, size = 3,
      colour = anchor$col
    ) +
    geom_text(
      data = anchor,
      aes(x, y, label = lab),
      nudge_x = 0.07,
      size = 3.2
    ) +
    labs(x = "u1 (whitened b1,b2)", y = "u2") +
    theme_pub(12)
}

direction_wheel_plot_base_labeled <- function(Base_A,
                                              cover = 1,
                                              pad_frac = 0.18,
                                              bandwidth_mult = 1.20,
                                              alpha_power = 0.58,
                                              edge_fade_frac = 0.10,
                                              label_size = 0.75,
                                              label_alpha = 0.5,
                                              repel = TRUE) {
  dfB <- data.frame(
    b1 = as.numeric(Base_A[, 1]),
    b2 = as.numeric(Base_A[, 2]),
    label = rownames(Base_A),
    stringsAsFactors = FALSE
  )
  dfB <- dfB[is.finite(dfB$b1) & is.finite(dfB$b2) & nzchar(dfB$label), , drop = FALSE]
  if (!nrow(dfB)) return(NULL)
  
  Xstd <- standardise_to_circle(Base_A, cover = cover)
  mu <- Xstd$mu
  S_hi <- Xstd$S_half_inv
  s <- Xstd$s
  
  nu <- 500
  
  rx <- range(dfB$b1, na.rm = TRUE)
  ry <- range(dfB$b2, na.rm = TRUE)
  wx <- diff(rx); wy <- diff(ry)
  if (!is.finite(wx) || wx <= 0) wx <- 1
  if (!is.finite(wy) || wy <= 0) wy <- 1
  
  lims <- c(
    rx[1] - pad_frac * wx, rx[2] + pad_frac * wx,
    ry[1] - pad_frac * wy, ry[2] + pad_frac * wy
  )
  
  h <- c(MASS::bandwidth.nrd(dfB$b1), MASS::bandwidth.nrd(dfB$b2)) *
    bandwidth_mult
  kd <- with(dfB, MASS::kde2d(b1, b2, n = nu, lims = lims, h = h))
  D <- kd$z
  D <- log1p(D / max(D, na.rm = TRUE))
  D <- D / stats::quantile(D, 0.99, na.rm = TRUE)
  D[D > 1] <- 1; D[D < 0] <- 0
  ALPHA <- D^alpha_power
  
  gx <- kd$x; gy <- kd$y
  G <- expand.grid(b1 = gx, b2 = gy)
  
  Uq <- (sweep(as.matrix(G), 2, mu, "-") %*% t(S_hi)) / s
  theta <- atan2(Uq[, 2], Uq[, 1])
  r <- sqrt(Uq[, 1]^2 + Uq[, 2]^2)
  
  H0 <- 170; L0 <- 60; Cmax <- 110; betaC <- 0.90
  H <- (H0 + theta * 180 / pi) %% 360
  C <- pmin(Cmax * (pmin(r, 1)^betaC), Cmax)
  L <- pmax(0, pmin(100, L0 - 6 * (pmin(r, 1)^1.1)))
  
  G$fill <- grDevices::hcl(H, C, L)
  G$alpha <- as.vector(ALPHA)

  edge_fade <- function(x, lo, hi, frac) {
    width <- max((hi - lo) * frac, .Machine$double.eps)
    z <- pmin(pmax(pmin(x - lo, hi - x) / width, 0), 1)
    sin(z * pi / 2)^2
  }
  G$alpha <- G$alpha *
    edge_fade(G$b1, lims[1], lims[2], edge_fade_frac) *
    edge_fade(G$b2, lims[3], lims[4], edge_fade_frac)
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(
      data = G,
      ggplot2::aes(b1, b2, fill = I(fill), alpha = alpha),
      interpolate = TRUE
    ) +
    ggplot2::scale_alpha(range = c(0, 1), guide = "none") +
    ggplot2::geom_point(
      data = dfB,
      ggplot2::aes(b1, b2),
      shape = 16, size = 0.7,
      colour = scales::alpha("black", 0.35)
    )
  
  if (repel && requireNamespace("ggrepel", quietly = TRUE)) {
    p <- p +
      ggrepel::geom_text_repel(
        data = dfB,
        ggplot2::aes(b1, b2, label = label),
        size = label_size,
        alpha = label_alpha,
        colour = "black",
        max.overlaps = Inf,
        min.segment.length = 0,
        box.padding = 0.15,
        point.padding = 0.05,
        segment.size = 0.15
      )
  } else {
    p <- p +
      ggplot2::geom_text(
        data = dfB,
        ggplot2::aes(b1, b2, label = label),
        size = label_size,
        alpha = label_alpha,
        colour = "black",
        vjust = -0.3
      )
  }
  
  p +
    ggplot2::coord_equal(
      xlim = lims[1:2], ylim = lims[3:4],
      expand = FALSE
      # clip = "off"  # remove this, rely on padding instead
    ) +
    ggplot2::labs(x = "b1", y = "b2") +
    theme_pub(12) +
    ggplot2::theme(plot.margin = ggplot2::margin(8, 12, 8, 8, "pt"))
}

write_density_diagnostics <- function(Base,
                                      Base_A,
                                      geom,
                                      X_selected,
                                      weights = NULL) {
  axis_diag <- function(df, space) {
    out <- lapply(names(df), function(nm) {
      x <- suppressWarnings(as.numeric(df[[nm]]))
      x <- x[is.finite(x)]
      if (!length(x)) return(NULL)
      rx <- range(x)
      width <- diff(rx)
      eps_rel <- if (is.finite(width) && width > 0) 1e-3 * width else 1e-3
      qs <- stats::quantile(
        x,
        probs = c(0.001, 0.01, 0.05, 0.5, 0.95, 0.99, 0.999),
        na.rm = TRUE,
        names = FALSE
      )
      data.frame(
        space = space,
        axis = nm,
        n = length(x),
        n_unique_round_6 = length(unique(round(x, 6))),
        n_unique_round_3 = length(unique(round(x, 3))),
        min = rx[1],
        q001 = qs[1],
        q01 = qs[2],
        q05 = qs[3],
        q50 = qs[4],
        q95 = qs[5],
        q99 = qs[6],
        q999 = qs[7],
        max = rx[2],
        n_at_min = sum(x == rx[1]),
        n_at_max = sum(x == rx[2]),
        n_near_min_0p1pct_range = sum(abs(x - rx[1]) <= eps_rel),
        n_near_max_0p1pct_range = sum(abs(x - rx[2]) <= eps_rel),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, out)
  }
  
  pred_diag <- lapply(names(X_selected), function(nm) {
    v <- X_selected[[nm]]
    vv <- v[!is.na(v)]
    tab <- table(vv)
    n_unique <- length(tab)
    data.frame(
      var = nm,
      weight = if (!is.null(weights) && nm %in% names(weights)) as.numeric(weights[[nm]]) else NA_real_,
      class = paste(class(v), collapse = "|"),
      is_numeric = is.numeric(v),
      is_factor = is.factor(v),
      is_ordered = is.ordered(v),
      n_observed = length(vv),
      n_unique = n_unique,
      min_level_count = if (n_unique) min(as.integer(tab)) else NA_integer_,
      max_level_count = if (n_unique) max(as.integer(tab)) else NA_integer_,
      stringsAsFactors = FALSE
    )
  })
  pred_diag <- do.call(rbind, pred_diag)
  
  coord_diag <- axis_diag(as.data.frame(Base), "base")
  if (!is.null(Base_A)) {
    coord_diag <- dplyr::bind_rows(coord_diag, axis_diag(as.data.frame(Base_A), "base_aligned"))
  }
  if (!is.null(geom$U)) {
    U_df <- as.data.frame(geom$U)
    keep <- intersect(c("u1", "u2"), names(U_df))
    if (length(keep)) {
      u_diag <- axis_diag(U_df[, keep, drop = FALSE], "unit_square")
      if ("insq" %in% names(U_df)) {
        u_diag$n_total_unit_rows <- nrow(U_df)
        u_diag$n_inside_unit_square <- sum(U_df$insq, na.rm = TRUE)
        u_diag$prop_inside_unit_square <- mean(U_df$insq, na.rm = TRUE)
      }
      coord_diag <- dplyr::bind_rows(coord_diag, u_diag)
    }
  }
  
  write_csv(pred_diag, "density_selected_predictor_diagnostics.csv")
  write_csv(coord_diag, "density_coordinate_diagnostics.csv")
  invisible(list(predictors = pred_diag, coordinates = coord_diag))
}

write_density_diagnostics(Base, Base_A, geom, X, w_all)

# Participant density plot
p_dens <- density_plots(geom)
save_plot_gg("FIG_dens_unitsquare_scatter", p_dens, width = 8.0, height = 7.0)

p_dens_b <- density_plots_base(Base_A)
if (!is.null(p_dens_b)) {
  save_plot_gg("FIG_dens_unitsquare_scatter_BASE", p_dens_b, width = 8.0, height = 7.0)
}

if (BASE_DIM >= 3L && ncol(Base) >= 3L) {
  p_dens_b13 <- density_plots_base_pair(Base, 1L, 3L, name = "Base density: b1 vs b3")
  p_dens_b23 <- density_plots_base_pair(Base, 2L, 3L, name = "Base density: b2 vs b3")
  if (!is.null(p_dens_b13)) {
    save_plot_gg("FIG_dens_BASE_b1b3", p_dens_b13, width = 8.0, height = 7.0)
  }
  if (!is.null(p_dens_b23)) {
    save_plot_gg("FIG_dens_BASE_b2b3", p_dens_b23, width = 8.0, height = 7.0)
  }
}

# Direction plot
p_dir <- direction_wheel_plot(geom)
if (!is.null(p_dir)) {
  save_plot_gg("FIG_uv_direction_density_HCL_smooth", p_dir, width = 8.0, height = 7.0)
}

p_dir_b <- direction_wheel_plot_base(Base_A)
if (!is.null(p_dir_b)) {
  save_plot_gg("FIG_uv_direction_density_HCL_smooth_BASE", p_dir_b, width = 8.0, height = 7.0)
}

if (BASE_DIM >= 3L && ncol(Base) >= 3L) {
  Base_A_13 <- Base[, c(1, 3), drop = FALSE]
  colnames(Base_A_13) <- c("b1", "b2")
  p_dir_b13 <- direction_wheel_plot_base(Base_A_13)
  if (!is.null(p_dir_b13)) {
    save_plot_gg("FIG_uv_direction_density_HCL_smooth_BASE_b1b3",
                 p_dir_b13, width = 8.0, height = 7.0)
  }
  
  Base_A_23 <- Base[, c(2, 3), drop = FALSE]
  colnames(Base_A_23) <- c("b1", "b2")
  p_dir_b23 <- direction_wheel_plot_base(Base_A_23)
  if (!is.null(p_dir_b23)) {
    save_plot_gg("FIG_uv_direction_density_HCL_smooth_BASE_b2b3",
                 p_dir_b23, width = 8.0, height = 7.0)
  }
  
  # Fog in whitened base space
  p_fog <- plot_density_fog_3d(
    X3 = Base_w,
    cols = 1:3,
    ng = 60,
    trace = "isosurface",
    h_mult = 3,
    isomin_q = 0.5,
    surface_count = 4,
    opacity = 0.4,
    show_points = TRUE,
    points_size = 3
  )
  p_fog
  
  # Save to disk (HTML widget)
  fog_file <- file.path(OUTPUTS_DIR, "FIG_fog_base_w.html")
  htmlwidgets::saveWidget(plotly::as_widget(p_fog), fog_file, selfcontained = FALSE)
  
}

p_dir_lab <- direction_wheel_plot_labeled(geom)
if (!is.null(p_dir_lab)) {
  save_plot_gg("FIG_uv_direction_density_HCL_smooth_LABELED",
               p_dir_lab,
               width = 12.0, height = 11.0)
}

p_dir_lab_b <- direction_wheel_plot_base_labeled(Base_A)
if (!is.null(p_dir_lab_b)) {
  save_plot_gg("FIG_uv_direction_density_HCL_smooth_LABELED_BASE",
               p_dir_lab_b, width = 12.0, height = 11.0)
}

# ------------------------------------------------------------------------------
# Biplots
# ------------------------------------------------------------------------------

score_item_1d <- function(nm, Z, varmap) {
  idx <- which(varmap == nm)
  if (!length(idx)) return(rep(NA_real_, nrow(Z)))
  v <- if (length(idx) == 1L) {
    as.numeric(Z[, idx])
  } else {
    sc <- try(suppressWarnings(prcomp(Z[, idx, drop = FALSE], rank. = 1)$x[, 1]), silent = TRUE)
    if (inherits(sc, "try-error")) rep(NA_real_, nrow(Z)) else as.numeric(sc)
  }
  as.numeric(scale(v))
}

score_item_1d_raw <- function(nm, Z, varmap) {
  idx <- which(varmap == nm)
  if (!length(idx)) return(rep(NA_real_, nrow(Z)))
  
  if (length(idx) == 1L) return(as.numeric(Z[, idx]))
  
  sc <- try(suppressWarnings(prcomp(Z[, idx, drop = FALSE],
                                    center = TRUE, scale. = FALSE,
                                    rank. = 1)$x[, 1]),
            silent = TRUE)
  if (inherits(sc, "try-error")) return(rep(NA_real_, nrow(Z)))
  as.numeric(sc)
}


pearson_r <- function(a, b) {
  if (length(a) != length(b)) return(NA_real_)
  suppressWarnings(cor(as.numeric(a), as.numeric(b), use = "complete.obs", method = "pearson"))
}

build_biplot_data <- function(Z_A, varmap, levelmap, Base_A, U) {
  B1 <- Base_A[, 1]; B2 <- Base_A[, 2]

  Rtab <- dplyr::bind_rows(lapply(seq_len(ncol(Z_A)), function(j) {
    v <- as.numeric(Z_A[, j])
    if (!any(is.finite(v))) return(NULL)

    nm <- as.character(varmap[j])
    response_level <- as.character(levelmap[j])
    if (is.na(response_level) || !nzchar(response_level)) response_level <- NA_character_
    sd_item <- stats::sd(v, na.rm = TRUE)
    if (!is.finite(sd_item) || sd_item <= 0) sd_item <- 0

    r_b1 <- pearson_r(v, B1); r_b2 <- pearson_r(v, B2)
    r_u1 <- pearson_r(v, U$u1); r_u2 <- pearson_r(v, U$u2)

    data.frame(
      item = nm,
      response_level = response_level,
      encoded_col = colnames(Z_A)[j],
      display_label = if (is.na(response_level)) nm else paste0(nm, " = ", response_level),
      r_b1 = r_b1, r_b2 = r_b2,
      r_u1 = r_u1, r_u2 = r_u2,
      sd_item = sd_item,
      stringsAsFactors = FALSE
    )
  }))
  
  if (!nrow(Rtab)) stop("[biplot] No item correlations computed.")
  
  Rtab |>
    dplyr::mutate(
      across(c(r_b1, r_b2, r_u1, r_u2, sd_item), ~ suppressWarnings(as.numeric(.))),
      r_b1 = ifelse(is.finite(r_b1), r_b1, 0),
      r_b2 = ifelse(is.finite(r_b2), r_b2, 0),
      r_u1 = ifelse(is.finite(r_u1), r_u1, 0),
      r_u2 = ifelse(is.finite(r_u2), r_u2, 0),
      sd_item = ifelse(is.finite(sd_item) & sd_item > 0, sd_item, 0),
      
      # correlation magnitudes (classic biplot)
      mag_r_base = sqrt(r_b1^2 + r_b2^2),
      mag_r_disk = sqrt(r_u1^2 + r_u2^2),
      
      # relevance-weighted vectors: same direction as r, length scaled by sd_item
      a_b1 = sd_item * r_b1,
      a_b2 = sd_item * r_b2,
      a_u1 = sd_item * r_u1,
      a_u2 = sd_item * r_u2,
      
      mag_a_base = sqrt(a_b1^2 + a_b2^2),
      mag_a_disk = sqrt(a_u1^2 + a_u2^2)
    )
}

build_biplot_data_3d <- function(Z_A, varmap, levelmap, Base,
                                 arrow = c("sd", "cov")) {
  arrow <- match.arg(arrow)
  stopifnot(ncol(Base) >= 3L)
  
  B1 <- suppressWarnings(as.numeric(Base[, 1]))
  B2 <- suppressWarnings(as.numeric(Base[, 2]))
  B3 <- suppressWarnings(as.numeric(Base[, 3]))
  
  # axis scales (used only if arrow == "cov")
  sdb1 <- stats::sd(B1, na.rm = TRUE)
  sdb2 <- stats::sd(B2, na.rm = TRUE)
  sdb3 <- stats::sd(B3, na.rm = TRUE)
  
  Rtab <- dplyr::bind_rows(lapply(seq_len(ncol(Z_A)), function(j) {
    v <- as.numeric(Z_A[, j])
    if (!any(is.finite(v))) return(NULL)

    nm <- as.character(varmap[j])
    response_level <- as.character(levelmap[j])
    if (is.na(response_level) || !nzchar(response_level)) response_level <- NA_character_
    r1 <- pearson_r(v, B1)
    r2 <- pearson_r(v, B2)
    r3 <- pearson_r(v, B3)
    
    sdv <- stats::sd(v, na.rm = TRUE)
    if (!is.finite(sdv)) sdv <- 0
    
    if (arrow == "cov") {
      a1 <- r1 * sdv * sdb1
      a2 <- r2 * sdv * sdb2
      a3 <- r3 * sdv * sdb3
    } else {
      # "sd" mode: direction = correlation, length scales with weighted item dispersion only
      a1 <- r1 * sdv
      a2 <- r2 * sdv
      a3 <- r3 * sdv
    }
    
    data.frame(
      item = nm,
      response_level = response_level,
      encoded_col = colnames(Z_A)[j],
      display_label = if (is.na(response_level)) nm else paste0(nm, " = ", response_level),
      r_b1 = r1, r_b2 = r2, r_b3 = r3,
      a_b1 = a1, a_b2 = a2, a_b3 = a3,
      sd_item = sdv,
      stringsAsFactors = FALSE
    )
  }))
  
  if (!nrow(Rtab)) stop("[biplot_3d] No item correlations computed.")
  
  Rtab |>
    dplyr::mutate(
      dplyr::across(c(r_b1, r_b2, r_b3, a_b1, a_b2, a_b3, sd_item),
                    ~ suppressWarnings(as.numeric(.))),
      dplyr::across(c(r_b1, r_b2, r_b3),
                    ~ ifelse(is.finite(.x), .x, 0)),
      dplyr::across(c(a_b1, a_b2, a_b3),
                    ~ ifelse(is.finite(.x), .x, 0)),
      mag3    = sqrt(r_b1^2 + r_b2^2 + r_b3^2),          # correlation magnitude (legacy)
      mag_a3  = sqrt(a_b1^2 + a_b2^2 + a_b3^2)           # relevance-scaled magnitude (use this)
    )
}

plot_biplots <- function(Rtab, Base_A, U, use = c("r", "a"),
                         top_per_octant = 6L,
                         arrow_head_cm = 0.10) {
  use <- match.arg(use)
  
  has_a <- all(c("a_b1","a_b2","a_u1","a_u2","mag_a_base","mag_a_disk") %in% names(Rtab))
  if (use == "a" && !has_a) {
    warning("[biplot] use='a' requested but a_* columns missing; falling back to correlations.")
    use <- "r"
  }
  
  # choose which vectors to plot
  if (use == "a") {
    vx_b1 <- "a_b1"; vx_b2 <- "a_b2"; vmag_base <- "mag_a_base"
    vx_u1 <- "a_u1"; vx_u2 <- "a_u2"; vmag_disk <- "mag_a_disk"
  } else {
    vx_b1 <- "r_b1"; vx_b2 <- "r_b2"; vmag_base <- "mag_r_base"
    vx_u1 <- "r_u1"; vx_u2 <- "r_u2"; vmag_disk <- "mag_r_disk"
  }
  
  B1n <- suppressWarnings(as.numeric(Base_A[, 1]))
  B2n <- suppressWarnings(as.numeric(Base_A[, 2]))
  B1n[!is.finite(B1n)] <- mean(B1n, na.rm = TRUE)
  B2n[!is.finite(B2n)] <- mean(B2n, na.rm = TRUE)
  
  Hidx <- grDevices::chull(B1n, B2n)
  H <- data.frame(b1 = B1n[Hidx], b2 = B2n[Hidx])
  
  cx <- mean(B1n); cy <- mean(B2n)
  Rlim <- 0.80 * min(diff(range(B1n)), diff(range(B2n)))
  
  get_octant <- function(x, y) {
    theta_deg <- atan2(y, x) * 180 / pi
    idx <- floor((theta_deg + 22.5 + 360) %% 360 / 45) + 1
    paste0("Octant_", idx)
  }

  prepare_display_levels <- function(tab, vx_col, vy_col, mag_col) {
    candidates <- tab |>
      dplyr::mutate(
        vx = .data[[vx_col]],
        vy = .data[[vy_col]],
        vmag = .data[[mag_col]],
        response_numeric = suppressWarnings(as.numeric(response_level)),
        short_item = sub("^.*\\.", "", item),
        short_label = ifelse(
          is.na(response_level) | !nzchar(response_level),
          short_item,
          paste0(short_item, " = ", response_level)
        )
      ) |>
      dplyr::filter(is.finite(vx), is.finite(vy), is.finite(vmag), vmag > 0) |>
      dplyr::group_by(item) |>
      dplyr::mutate(n_levels = dplyr::n_distinct(response_level)) |>
      dplyr::ungroup()

    binary_levels <- candidates |>
      dplyr::filter(n_levels == 2L) |>
      dplyr::group_by(item) |>
      dplyr::arrange(dplyr::desc(response_numeric),
                     dplyr::desc(response_level), .by_group = TRUE) |>
      dplyr::slice_head(n = 1L) |>
      dplyr::ungroup()

    multilevels <- candidates |>
      dplyr::filter(n_levels != 2L)

    dplyr::bind_rows(binary_levels, multilevels) |>
      dplyr::mutate(octant = get_octant(vx, vy))
  }

  select_display_levels <- function(tab, vx_col, vy_col, mag_col) {
    prepare_display_levels(tab, vx_col, vy_col, mag_col) |>
      dplyr::group_by(octant) |>
      dplyr::group_modify(function(.x, .y) {
        ranked <- .x |>
          dplyr::arrange(dplyr::desc(vmag), item, response_level) |>
          dplyr::distinct(item, .keep_all = TRUE)

        n_keep <- min(as.integer(top_per_octant), nrow(ranked))
        n_overall <- max(0L, n_keep - 1L)
        selected <- utils::head(ranked, n_overall) |>
          dplyr::mutate(selection_stratum = "top_overall")
        remaining <- ranked |>
          dplyr::filter(!item %in% selected$item)

        if (nrow(selected) && any(selected$n_levels > 2L)) {
          extra <- utils::head(remaining, 1L) |>
            dplyr::mutate(selection_stratum = "next_overall")
        } else {
          extra <- remaining |>
            dplyr::filter(n_levels > 2L) |>
            utils::head(1L) |>
            dplyr::mutate(selection_stratum = "multilevel_anchor")
          if (!nrow(extra)) {
            extra <- utils::head(remaining, 1L) |>
              dplyr::mutate(selection_stratum = "next_overall_no_multilevel")
          }
        }

        dplyr::bind_rows(selected, extra) |>
          utils::head(n_keep)
      }) |>
      dplyr::arrange(dplyr::desc(vmag), .by_group = TRUE) |>
      dplyr::mutate(octant_rank = dplyr::row_number()) |>
      dplyr::ungroup()
  }

  place_labels_on_ring <- function(tab, center_x, center_y, outer_radius) {
    tab |>
      dplyr::mutate(
        octant_number = as.integer(sub("^Octant_", "", octant)),
        octant_center = 45 * (octant_number - 1L),
        vector_angle = atan2(vy, vx) * 180 / pi,
        angle_from_center = ((vector_angle - octant_center + 180) %% 360) - 180
      ) |>
      dplyr::group_by(octant) |>
      dplyr::arrange(angle_from_center, .by_group = TRUE) |>
      dplyr::mutate(
        slot_index = dplyr::row_number(),
        slot_count = dplyr::n(),
        slot_offset = ifelse(
          slot_count == 1L,
          0,
          -18 + 36 * (slot_index - 1L) / (slot_count - 1L)
        ),
        label_angle = (octant_center + slot_offset) * pi / 180,
        label_radius = outer_radius * (1 + 0.045 * ((slot_index - 1L) %% 2L)),
        label_x = center_x + label_radius * cos(label_angle),
        label_y = center_y + label_radius * sin(label_angle),
        label_angle_signed = ((label_angle * 180 / pi + 180) %% 360) - 180,
        label_text_angle = dplyr::case_when(
          label_angle_signed > 90 ~ label_angle_signed - 180,
          label_angle_signed < -90 ~ label_angle_signed + 180,
          TRUE ~ label_angle_signed
        ),
        label_hjust = ifelse(abs(label_angle_signed) <= 90, 0, 1),
        label_vjust = 0.5
      ) |>
      dplyr::ungroup()
  }
  
  # --- BASE ---
  max_mag_base <- max(Rtab[[vmag_base]], na.rm = TRUE)
  if (!is.finite(max_mag_base) || max_mag_base <= 0) max_mag_base <- 1
  Rscale_base <- Rlim / max_mag_base

  A_base <- prepare_display_levels(Rtab, vx_b1, vx_b2, vmag_base) |>
    dplyr::mutate(
      x0 = cx, y0 = cy,
      x1 = cx + Rscale_base * vx,
      y1 = cy + Rscale_base * vy
    )
  S_base <- select_display_levels(Rtab, vx_b1, vx_b2, vmag_base) |>
    dplyr::mutate(
      x0 = cx, y0 = cy,
      x1 = cx + Rscale_base * vx,
      y1 = cy + Rscale_base * vy
    )
  base_outer_radius <- 1.12 * max(
    sqrt((B1n - cx)^2 + (B2n - cy)^2),
    sqrt((S_base$x1 - cx)^2 + (S_base$y1 - cy)^2),
    na.rm = TRUE
    )
  S_base <- place_labels_on_ring(S_base, cx, cy, base_outer_radius)
  A_base_unlabelled <- A_base |>
    dplyr::anti_join(
      S_base |> dplyr::select(item, response_level),
      by = c("item", "response_level")
    )
  
  p_base <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      data = H, ggplot2::aes(b1, b2),
      fill = NA, colour = "black", linewidth = 0.4
    ) +
    ggplot2::geom_segment(
      data = A_base_unlabelled,
      ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1),
      linewidth = 0.18, colour = "grey35", alpha = 0.10,
      arrow = grid::arrow(length = grid::unit(arrow_head_cm * 0.45, "cm"),
                          type = "closed")
    ) +
    ggplot2::geom_segment(
      data = S_base,
      ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1),
      linewidth = 0.65, colour = "firebrick",
      arrow = grid::arrow(length = grid::unit(arrow_head_cm, "cm"),
                          type = "closed")
    ) +
    ggplot2::geom_segment(
      data = S_base,
      ggplot2::aes(x = x1, y = y1, xend = label_x, yend = label_y),
      linewidth = 0.25, colour = "grey45"
    ) +
    ggplot2::geom_text(
      data = S_base,
      ggplot2::aes(x = label_x, y = label_y, label = short_label,
                   hjust = label_hjust, vjust = label_vjust,
                   angle = label_text_angle),
      size = 2.65
    ) +
    ggplot2::coord_equal(clip = "off") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.13)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.13)) +
    ggplot2::labs(x = "b1", y = "b2") +
    theme_pub(12) +
    ggplot2::theme(plot.margin = ggplot2::margin(12, 75, 12, 75))
  
  # --- UNIT DISK ---
  draw_disk_outline <- function() {
    th <- seq(0, 2 * pi, length.out = 361)
    data.frame(x = cos(th), y = sin(th))
  }
  Rdisk <- 0.85
  
  max_mag_disk <- max(Rtab[[vmag_disk]], na.rm = TRUE)
  if (!is.finite(max_mag_disk) || max_mag_disk <= 0) max_mag_disk <- 1
  Rscale_disk <- Rdisk / max_mag_disk

  A_disk <- prepare_display_levels(Rtab, vx_u1, vx_u2, vmag_disk) |>
    dplyr::mutate(
      u0 = 0, v0 = 0,
      u1 = Rscale_disk * vx,
      v1 = Rscale_disk * vy
    )
  S_disk <- select_display_levels(Rtab, vx_u1, vx_u2, vmag_disk) |>
    dplyr::mutate(
      u0 = 0, v0 = 0,
      u1 = Rscale_disk * vx,
      v1 = Rscale_disk * vy
    )
  S_disk <- place_labels_on_ring(S_disk, 0, 0, 1.16)
  A_disk_unlabelled <- A_disk |>
    dplyr::anti_join(
      S_disk |> dplyr::select(item, response_level),
      by = c("item", "response_level")
    )
  
  p_disk <- ggplot2::ggplot() +
    ggplot2::geom_path(data = draw_disk_outline(), ggplot2::aes(x, y)) +
    ggplot2::geom_segment(
      data = A_disk_unlabelled,
      ggplot2::aes(x = u0, y = v0, xend = u1, yend = v1),
      linewidth = 0.18, colour = "grey35", alpha = 0.10,
      arrow = grid::arrow(length = grid::unit(arrow_head_cm * 0.45, "cm"),
                          type = "closed")
    ) +
    ggplot2::geom_segment(
      data = S_disk,
      ggplot2::aes(x = u0, y = v0, xend = u1, yend = v1),
      linewidth = 0.65, colour = "firebrick",
      arrow = grid::arrow(length = grid::unit(arrow_head_cm, "cm"), type = "closed")
    ) +
    ggplot2::geom_segment(
      data = S_disk,
      ggplot2::aes(x = u1, y = v1, xend = label_x, yend = label_y),
      linewidth = 0.25, colour = "grey45"
    ) +
    ggplot2::geom_text(
      data = S_disk,
      ggplot2::aes(x = label_x, y = label_y, label = short_label,
                   hjust = label_hjust, vjust = label_vjust,
                   angle = label_text_angle),
      size = 2.65
    ) +
    ggplot2::coord_equal(xlim = c(-1.55, 1.55), ylim = c(-1.55, 1.55),
                         expand = FALSE, clip = "off") +
    ggplot2::labs(x = "u1 (whitened b1,b2)", y = "u2") +
    theme_pub(12) +
    ggplot2::theme(plot.margin = ggplot2::margin(12, 75, 12, 75))
  
  list(
    p_base = p_base,
    p_disk = p_disk,
    all_base = A_base,
    all_disk = A_disk,
    selected_base = S_base,
    selected_disk = S_disk
  )
}

plot_biplots_3d <- function(R3, Base, top_global = 32L,
                            use = c("a", "r")) {
  suppressPackageStartupMessages(requireNamespace("ggrepel", quietly = TRUE))
  use <- match.arg(use)
  
  B1 <- suppressWarnings(as.numeric(Base[, 1]))
  B2 <- suppressWarnings(as.numeric(Base[, 2]))
  B3 <- suppressWarnings(as.numeric(Base[, 3]))
  
  # choose vectors + magnitude
  has_a <- all(c("a_b1", "a_b2", "a_b3", "mag_a3") %in% names(R3))
  if (use == "a" && !has_a) use <- "r"
  
  if (use == "a") {
    v1 <- R3$a_b1; v2 <- R3$a_b2; v3 <- R3$a_b3
    mag <- R3$mag_a3
  } else {
    v1 <- R3$r_b1; v2 <- R3$r_b2; v3 <- R3$r_b3
    mag <- R3$mag3
  }
  
  # top items globally by chosen magnitude
  o <- order(-mag)
  R3 <- R3[o, , drop = FALSE]
  v1 <- v1[o]; v2 <- v2[o]; v3 <- v3[o]; mag <- mag[o]
  if (nrow(R3) > top_global) {
    R3 <- R3[1:top_global, , drop = FALSE]
    v1 <- v1[1:top_global]; v2 <- v2[1:top_global]; v3 <- v3[1:top_global]; mag <- mag[1:top_global]
  }
  
  # common centres
  cx1 <- mean(B1, na.rm = TRUE)
  cx2 <- mean(B2, na.rm = TRUE)
  cx3 <- mean(B3, na.rm = TRUE)
  
  rx <- diff(range(B1, na.rm = TRUE))
  ry <- diff(range(B2, na.rm = TRUE))
  rz <- diff(range(B3, na.rm = TRUE))
  
  # plane-wise scaling so the longest arrow in that plane uses ~80% of the available span
  len12 <- sqrt(v1^2 + v2^2); m12 <- max(len12, na.rm = TRUE); if (!is.finite(m12) || m12 <= 1e-12) m12 <- 1
  len13 <- sqrt(v1^2 + v3^2); m13 <- max(len13, na.rm = TRUE); if (!is.finite(m13) || m13 <= 1e-12) m13 <- 1
  len23 <- sqrt(v2^2 + v3^2); m23 <- max(len23, na.rm = TRUE); if (!is.finite(m23) || m23 <= 1e-12) m23 <- 1
  
  S12 <- data.frame(
    item = R3$display_label,
    x0 = cx1, y0 = cx2,
    x1 = cx1 + (0.80 * min(rx, ry) / m12) * v1,
    y1 = cx2 + (0.80 * min(rx, ry) / m12) * v2
  )
  
  S13 <- data.frame(
    item = R3$display_label,
    x0 = cx1, y0 = cx3,
    x1 = cx1 + (0.80 * min(rx, rz) / m13) * v1,
    y1 = cx3 + (0.80 * min(rx, rz) / m13) * v3
  )
  
  S23 <- data.frame(
    item = R3$display_label,
    x0 = cx2, y0 = cx3,
    x1 = cx2 + (0.80 * min(ry, rz) / m23) * v2,
    y1 = cx3 + (0.80 * min(ry, rz) / m23) * v3
  )
  
  hull_df <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
    if (!length(x)) return(NULL)
    idx <- grDevices::chull(x, y)
    data.frame(x = x[idx], y = y[idx])
  }
  
  H12 <- hull_df(B1, B2)
  H13 <- hull_df(B1, B3)
  H23 <- hull_df(B2, B3)
  
  p12 <- ggplot2::ggplot() +
    { if (!is.null(H12)) ggplot2::geom_polygon(data = H12, ggplot2::aes(x, y),
                                               fill = NA, colour = "black", linewidth = 0.4) } +
    ggplot2::geom_segment(
      data = S12,
      ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1),
      linewidth = 0.8, colour = "firebrick",
      arrow = grid::arrow(length = grid::unit(0.14, "cm"))
    ) +
    ggrepel::geom_label_repel(
      data = S12,
      ggplot2::aes(x = x1, y = y1, label = item),
      size = 3.1, max.overlaps = Inf,
      label.size = 0, label.padding = grid::unit(0.10, "lines")
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "b1", y = "b2") +
    theme_pub(12)
  
  p13 <- ggplot2::ggplot() +
    { if (!is.null(H13)) ggplot2::geom_polygon(data = H13, ggplot2::aes(x, y),
                                               fill = NA, colour = "black", linewidth = 0.4) } +
    ggplot2::geom_segment(
      data = S13,
      ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1),
      linewidth = 0.8, colour = "firebrick",
      arrow = grid::arrow(length = grid::unit(0.14, "cm"))
    ) +
    ggrepel::geom_label_repel(
      data = S13,
      ggplot2::aes(x = x1, y = y1, label = item),
      size = 3.1, max.overlaps = Inf,
      label.size = 0, label.padding = grid::unit(0.10, "lines")
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "b1", y = "b3") +
    theme_pub(12)
  
  p23 <- ggplot2::ggplot() +
    { if (!is.null(H23)) ggplot2::geom_polygon(data = H23, ggplot2::aes(x, y),
                                               fill = NA, colour = "black", linewidth = 0.4) } +
    ggplot2::geom_segment(
      data = S23,
      ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1),
      linewidth = 0.8, colour = "firebrick",
      arrow = grid::arrow(length = grid::unit(0.14, "cm"))
    ) +
    ggrepel::geom_label_repel(
      data = S23,
      ggplot2::aes(x = x1, y = y1, label = item),
      size = 3.1, max.overlaps = Inf,
      label.size = 0, label.padding = grid::unit(0.10, "lines")
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "b2", y = "b3") +
    theme_pub(12)
  
  list(p12 = p12, p13 = p13, p23 = p23)
}

# --- Save helper for plotly widgets ---
save_plotly_html <- function(stem, p,
                             outdir = OUTPUTS_DIR,
                             width = 1100, height = 900,
                             selfcontained = FALSE) {
  if (!requireNamespace("htmlwidgets", quietly = TRUE)) stop("Need package: htmlwidgets")
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  f <- file.path(outdir, paste0(stem, ".html"))
  htmlwidgets::saveWidget(p, file = f, selfcontained = selfcontained, background = "white")
  invisible(f)
}

# --- 3D octant by sign (coverage analogue to 2D octants) ---
octant3d <- function(x, y, z) {
  paste0("O_",
         ifelse(x >= 0, "p", "n"),
         ifelse(y >= 0, "p", "n"),
         ifelse(z >= 0, "p", "n"))
}

# --- 3D interactive biplot in Base space with true arrowheads (cones) + labels ---
plot_biplot_base_3d_interactive <- function(R3, Base,
                                            top_per_octant = 6L,
                                            show_hull = TRUE,
                                            hull_opacity = 0.06,
                                            points_size = 3.0,
                                            points_opacity = 0.28,
                                            arrow_width = 6,
                                            label_mode = c("hover", "text", "both"),
                                            label_font_size = 12,
                                            label_offset_frac = 0.03,
                                            head_len_frac = 0.06,
                                            head_sizeref = 0.45,
                                            use = c("a", "r")) {
  
  if (!requireNamespace("plotly", quietly = TRUE)) stop("Need package: plotly")
  label_mode <- match.arg(label_mode)
  use <- match.arg(use)
  stopifnot(ncol(Base) >= 3L)
  
  B <- as.data.frame(Base[, 1:3, drop = FALSE])
  names(B) <- c("b1", "b2", "b3")
  for (nm in names(B)) {
    B[[nm]] <- suppressWarnings(as.numeric(B[[nm]]))
    B[[nm]][!is.finite(B[[nm]])] <- mean(B[[nm]], na.rm = TRUE)
  }
  
  R3 <- as.data.frame(R3)
  
  has_a <- all(c("a_b1", "a_b2", "a_b3", "mag_a3") %in% names(R3))
  if (use == "a" && !has_a) use <- "r"
  
  # octants by sign: use correlations (signs match a_* anyway)
  R3$oct <- octant3d(R3$r_b1, R3$r_b2, R3$r_b3)
  
  # rank per-octant by relevance-scaled magnitude if available
  mag <- if (use == "a" && "mag_a3" %in% names(R3)) R3$mag_a3 else R3$mag3
  
  R3s <- R3 |>
    dplyr::mutate(.mag = mag) |>
    dplyr::group_by(oct) |>
    dplyr::arrange(dplyr::desc(.mag)) |>
    dplyr::distinct(item, .keep_all = TRUE) |>
    dplyr::slice_head(n = as.integer(top_per_octant)) |>
    dplyr::ungroup()
  
  # choose vector columns for geometry
  if (use == "a") {
    v1 <- R3s$a_b1; v2 <- R3s$a_b2; v3 <- R3s$a_b3
  } else {
    v1 <- R3s$r_b1; v2 <- R3s$r_b2; v3 <- R3s$r_b3
  }
  
  cx <- mean(B$b1, na.rm = TRUE)
  cy <- mean(B$b2, na.rm = TRUE)
  cz <- mean(B$b3, na.rm = TRUE)
  
  rx <- diff(range(B$b1, na.rm = TRUE))
  ry <- diff(range(B$b2, na.rm = TRUE))
  rz <- diff(range(B$b3, na.rm = TRUE))
  
  # scale so max 3D vector length sits inside data extent
  Lvec <- sqrt(v1^2 + v2^2 + v3^2)
  mL <- max(Lvec, na.rm = TRUE)
  if (!is.finite(mL) || mL <= 1e-12) mL <- 1
  Rscale <- 0.80 * min(rx, ry, rz) / mL
  
  A <- data.frame(
    item = R3s$display_label,
    r1 = R3s$r_b1, r2 = R3s$r_b2, r3 = R3s$r_b3,
    x0 = cx, y0 = cy, z0 = cz,
    x1 = cx + Rscale * v1,
    y1 = cy + Rscale * v2,
    z1 = cz + Rscale * v3
  )
  
  # segment trace uses NA breaks
  xseg <- as.vector(rbind(A$x0, A$x1, NA))
  yseg <- as.vector(rbind(A$y0, A$y1, NA))
  zseg <- as.vector(rbind(A$z0, A$z1, NA))
  
  dx <- A$x1 - A$x0; dy <- A$y1 - A$y0; dz <- A$z1 - A$z0
  L  <- sqrt(dx^2 + dy^2 + dz^2)
  L[L <= 1e-12] <- 1
  ux <- dx / L; uy <- dy / L; uz <- dz / L
  
  head_len <- head_len_frac * min(rx, ry, rz)
  u_cone <- ux * head_len
  v_cone <- uy * head_len
  w_cone <- uz * head_len
  
  off <- label_offset_frac * min(rx, ry, rz)
  xl <- A$x1 + ux * off
  yl <- A$y1 + uy * off
  zl <- A$z1 + uz * off
  
  p <- plotly::plot_ly() |>
    plotly::add_markers(
      data = B,
      x = ~b1, y = ~b2, z = ~b3,
      type = "scatter3d", mode = "markers",
      marker = list(size = points_size, opacity = points_opacity),
      name = "Subjects",
      hoverinfo = "skip"
    )
  
  if (isTRUE(show_hull) && requireNamespace("geometry", quietly = TRUE)) {
    P <- as.matrix(B[, c("b1", "b2", "b3")]); storage.mode(P) <- "double"
    tri <- try(geometry::convhulln(P, options = "Qt"), silent = TRUE)
    if (!inherits(tri, "try-error")) {
      p <- p |>
        plotly::add_trace(
          type = "mesh3d",
          x = P[, 1], y = P[, 2], z = P[, 3],
          i = tri[, 1] - 1, j = tri[, 2] - 1, k = tri[, 3] - 1,
          opacity = hull_opacity,
          name = "Hull",
          showscale = FALSE,
          hoverinfo = "skip"
        )
    }
  }
  
  p <- p |>
    plotly::add_trace(
      x = xseg, y = yseg, z = zseg,
      type = "scatter3d",
      mode = "lines",
      line = list(width = arrow_width, color = "firebrick"),
      name = "Items",
      hoverinfo = "skip"
    ) |>
    plotly::add_trace(
      type = "cone",
      x = A$x1, y = A$y1, z = A$z1,
      u = u_cone, v = v_cone, w = w_cone,
      anchor = "tip",
      sizemode = "absolute",
      sizeref = head_sizeref,
      showscale = FALSE,
      name = "Arrowheads",
      hoverinfo = "skip",
      colorscale = list(c(0, "firebrick"), c(1, "firebrick"))
    ) |>
    plotly::add_markers(
      data = A,
      x = ~x1, y = ~y1, z = ~z1,
      type = "scatter3d", mode = "markers",
      marker = list(size = 3.8, color = "firebrick", opacity = 0.85),
      text = ~sprintf("%s<br>r(b1)=%.3f r(b2)=%.3f r(b3)=%.3f", item, r1, r2, r3),
      hoverinfo = "text",
      name = "Item tips"
    )
  
  if (label_mode %in% c("text", "both")) {
    p <- p |>
      plotly::add_trace(
        type = "scatter3d",
        mode = "text",
        x = xl, y = yl, z = zl,
        text = A$item,
        textfont = list(size = label_font_size, color = "black"),
        showlegend = FALSE,
        hoverinfo = "skip"
      )
  }
  
  p |>
    plotly::layout(
      scene = list(
        xaxis = list(title = "b1"),
        yaxis = list(title = "b2"),
        zaxis = list(title = "b3"),
        aspectmode = "data"
      ),
      legend = list(orientation = "h", x = 0.02, y = 0.98)
    ) |>
    plotly::config(displayModeBar = TRUE)
}

# Biplots
if (!is.null(varmap)) {
  Rtab <- build_biplot_data(Z_A, varmap, levelmap, Base_A, geom$U)
  write_csv(Rtab, "items_vs_base_and_unitdisk_correlations.csv")
  
  BIP  <- plot_biplots(Rtab, Base_A, geom$U, use = "r")      # response-level correlation vectors
  write_csv(BIP$selected_base, "biplot_selected_arrows_BASE.csv")
  write_csv(BIP$selected_disk, "biplot_selected_arrows_UNITDISK.csv")
  save_plot_gg("FIG_biplot_items_BASE", BIP$p_base, width = 9.0, height = 9.0)
  save_plot_gg("FIG_biplot_items_UNITDISK", BIP$p_disk, width = 9.0, height = 9.0)
  
  if (BASE_DIM >= 3L && ncol(Base) >= 3L) {
    R3  <- build_biplot_data_3d(Z_A, varmap, levelmap, Base)
    BIP3 <- plot_biplots_3d(R3, Base, top_global = 32L, use = "a")  # use a_* for lengths
    save_plot_gg("FIG_biplot3D_items_b1b2", BIP3$p12, width = 8.0, height = 7.0)
    save_plot_gg("FIG_biplot3D_items_b1b3", BIP3$p13, width = 8.0, height = 7.0)
    save_plot_gg("FIG_biplot3D_items_b2b3", BIP3$p23, width = 8.0, height = 7.0)
    
    p_bip3d <- plot_biplot_base_3d_interactive(
      R3, Base,
      top_per_octant = 4L,
      label_mode = "both",
      use = "a"
    )
    
    save_plotly_html("FIG_biplot_items_BASE_3D", p_bip3d, outdir = OUTPUTS_DIR)
  }
} else {
  msgf("[biplot] varmap missing; skipping biplots.")
}
