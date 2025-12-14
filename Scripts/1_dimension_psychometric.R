# ==============================================================================
#                               dimension_psychometric.R
# ==============================================================================

# ==============================================================================
# 1. Helper Functions
# ==============================================================================

# --- Data Preparation for Gower Distance ---
prep_X_for_gower <- function(X, rare_prop = 0.01, do_jitter = TRUE) {
  X1 <- as.data.frame(X, check.names = TRUE, stringsAsFactors = FALSE)
  
  # Ensure character cols are factors
  for (nm in names(X1)) {
    v <- X1[[nm]]
    if (is.character(v)) X1[[nm]] <- factor(v)
  }
  
  # Collapses rare factor levels into NA or a residual category
  drop_rare <- function(f, prop) {
    if (!is.factor(f) || is.ordered(f)) {
      return(f)
    }
    tb <- prop.table(table(f))
    keep <- names(tb)[tb >= prop]
    f <- factor(ifelse(f %in% keep, as.character(f), NA), exclude = NULL)
    droplevels(f)
  }
  
  X1 <- as.data.frame(lapply(X1, drop_rare, prop = rare_prop), stringsAsFactors = FALSE)
  
  # Add microscopic jitter to numerics to prevent ties in distance
  if (isTRUE(do_jitter)) {
    for (nm in names(X1)) {
      if (is.numeric(X1[[nm]])) {
        sdv <- stats::sd(X1[[nm]], na.rm = TRUE)
        if (is.finite(sdv) && sdv > 0) {
          X1[[nm]] <- X1[[nm]] + rnorm(length(X1[[nm]]), 0, 1e-6 * sdv)
        }
      }
    }
  }
  
  # Identify column types for daisy()
  ord_cols <- names(X1)[vapply(X1, is.ordered, logical(1))]
  fac_cols <- names(X1)[vapply(X1, function(z) is.factor(z) && !is.ordered(z), logical(1))]
  
  # Internal helper to detect binary
  .is_binary <- function(x) length(unique(na.omit(x))) == 2
  bin_cols <- fac_cols[vapply(X1[fac_cols], .is_binary, logical(1))]
  
  type_list <- list()
  if (length(bin_cols)) type_list$asymm <- bin_cols
  if (length(ord_cols)) type_list$ordratio <- ord_cols
  
  w <- rep(1, ncol(X1))
  names(w) <- names(X1)
  
  list(X = X1, type = type_list, weights = w)
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
.two_nn_from_distvec <- function(D, eps = .Machine$double.eps) {
  # Convert compact dist vector to full matrix (Internal C function, very fast)
  # Memory note: For N=1000, this is only ~8MB.
  M <- as.matrix(D)
  n <- nrow(M)
  
  if (n < 2L) return(list(d1 = numeric(0), d2 = numeric(0)))
  
  # Set diagonal to Inf so we don't find the point itself
  diag(M) <- Inf
  
  # Find 1st NN (min distance is max negative distance)
  negM <- -M
  idx1 <- max.col(negM, ties.method = "first")
  d1 <- M[cbind(seq_len(n), idx1)]
  
  # Find 2nd NN by masking the first
  negM[cbind(seq_len(n), idx1)] <- -Inf
  idx2 <- max.col(negM, ties.method = "first")
  d2 <- M[cbind(seq_len(n), idx2)]
  
  # Safety clamp
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

optimise_gower_weights_constrained <- function(X, init_weights, allow_update,
                                               objective = "TwoNN_all",
                                               w_min = W_MIN,
                                               step_grid = W_STEP_GRID, 
                                               batch_k = W_BATCH_K,
                                               batch_factor = W_BATCH_FACTOR,
                                               max_iter = W_MAX_ITERS,
                                               n_rows_sub = N_ROWS_SUB,
                                               ncores = NULL, 
                                               seed_jitter = SEED_JITTER,
                                               reps_idx = NULL,
                                               core_idx_rep = NULL,
                                               verbose = TRUE,
                                               plot_progress = TRUE) {
  
  # Internal ID calc from numerator/denominator
  calc_id_fast <- function(num, den, n_rows) {
    Dvec <- num / pmax(den, .Machine$double.eps)
    attr(Dvec, "Size")  <- n_rows
    attr(Dvec, "Diag")  <- FALSE
    attr(Dvec, "Upper") <- FALSE
    class(Dvec) <- "dist"
    twonn_id_from_dist(Dvec) 
  }
  
  # Prepare subset
  px  <- prep_X_for_gower(X, rare_prop = RARE_LEVEL_MIN_PROP, do_jitter = TRUE)
  X0  <- px$X
  typ <- px$type
  
  row_pool <- if (!is.null(reps_idx)) reps_idx else seq_len(nrow(X0))
  eff_n_sub <- if (is.null(n_rows_sub)) 1000 else min(n_rows_sub, 1500)
  
  ix_sub_from_reps <- if (length(row_pool) > eff_n_sub) {
    if (isTRUE(FIX_REP_SUBSET)) head(row_pool, eff_n_sub) else sample(row_pool, eff_n_sub)
  } else row_pool
  
  Xs <- X0[ix_sub_from_reps, , drop = FALSE]
  n_sub <- nrow(Xs)
  
  vars <- colnames(Xs)
  w <- init_weights[vars]
  allow_update <- allow_update[vars]
  w[!allow_update] <- pmax(w_min, w[!allow_update])
  
  # Build cache
  cache <- make_NS_cache(Xs, type = typ)
  
  # Initialize State
  num_cur <- Reduce(`+`, Map(`*`, cache$N, as.list(w)))
  den_cur <- Reduce(`+`, Map(`*`, cache$S, as.list(w)))
  
  id0 <- calc_id_fast(num_cur, den_cur, n_sub)
  hist <- data.frame(iter = 0L, ID = id0, changed = NA_character_, note = NA_character_)
  
  if (verbose) cat(sprintf("[optim] Start ID: %.3f | N_sub: %d | Mode: STOCHASTIC\n", id0, n_sub))
  
  max_iter_eff <- if (is.null(max_iter) || !is.finite(max_iter)) 1000L else as.integer(max_iter)
  
  # Optimization Config
  N_SAMPLE_PER_ITER <- 50 
  step_grid <- sort(step_grid, decreasing = FALSE) 
  hot_vars <- integer(0) # Track momentum
  
  id <- id0
  
  for (it in seq_len(max_iter_eff)) {
    
    can_all <- which(allow_update & (w > w_min + 1e-12))
    if (!length(can_all)) break
    
    # 1. Stochastic Sampling: Keep hot vars, fill rest with randoms
    n_rnd <- max(10, N_SAMPLE_PER_ITER - length(hot_vars))
    rnd_vars <- sample(can_all, min(length(can_all), n_rnd))
    can_iter <- unique(c(hot_vars, rnd_vars))
    
    # 2. Evaluate Candidates (Greedy)
    best_res <- list(id = Inf)
    
    for (j in can_iter) {
      w_base <- w[j]
      
      for (factor in step_grid) {
        w_new <- max(w_min, w_base * factor)
        if (w_new >= w_base - 1e-12) next 
        
        delta <- w_new - w_base
        val <- calc_id_fast(num_cur + delta * cache$N[[j]], 
                            den_cur + delta * cache$S[[j]], n_sub)
        
        if (is.finite(val) && val < best_res$id) {
          best_res <- list(id = val, j = j, w_new = w_new)
        }
      }
    }
    
    # 3. Commit Best Move
    changed <- FALSE
    if (is.finite(best_res$id) && best_res$id < id - 1e-6) {
      jbest <- best_res$j
      wbest <- best_res$w_new
      delta <- wbest - w[jbest]
      
      num_cur <- num_cur + delta * cache$N[[jbest]]
      den_cur <- den_cur + delta * cache$S[[jbest]]
      w[jbest] <- wbest
      id <- best_res$id
      changed <- TRUE
      
      # Update "Hot" list (Keep this var, drop oldest if too full)
      hot_vars <- unique(c(jbest, hot_vars))
      if (length(hot_vars) > 10) hot_vars <- head(hot_vars, 10)
      
      hist <- rbind(hist, data.frame(iter = it, ID = id, changed = vars[jbest], note = sprintf("%.3f", wbest)))
      if (verbose) cat(sprintf("   iter %d: %s -> %.3f (ID: %.3f)\n", it, vars[jbest], wbest, id))
      
    } else {
      # No improvement? Flush momentum or quit
      if (length(hot_vars) > 0) {
        hot_vars <- integer(0) 
        if (verbose) cat("   [optim] Momentum lost, flushing hot vars.\n")
      } else {
        if (verbose) cat("[optim] No improvement in random subset.\n")
        break
      }
    }
    
    # 4. Batch Descent (Every 5th iter)
    if (batch_k > 1 && (it %% 5 == 0)) {
      remain <- setdiff(can_all, best_res$j)
      remain <- sample(remain, min(length(remain), N_SAMPLE_PER_ITER))
      
      if (length(remain) > 0) {
        scores <- numeric(length(remain))
        for (i in seq_along(remain)) {
          j <- remain[i]
          delta_b <- (max(w_min, w[j] * batch_factor)) - w[j]
          scores[i] <- calc_id_fast(num_cur + delta_b * cache$N[[j]], 
                                    den_cur + delta_b * cache$S[[j]], n_sub)
        }
        
        ord <- order(scores)
        take_idx <- head(ord, min(batch_k, length(ord)))
        take_vars <- remain[take_idx]
        
        if (length(take_vars)) {
          num_b <- num_cur; den_b <- den_cur
          w_b <- w
          for (j in take_vars) {
            wn <- max(w_min, w[j] * batch_factor)
            d  <- wn - w[j]
            num_b <- num_b + d * cache$N[[j]]
            den_b <- den_b + d * cache$S[[j]]
            w_b[j] <- wn
          }
          id_b <- calc_id_fast(num_b, den_b, n_sub)
          
          if (is.finite(id_b) && id_b < id - 1e-6) {
            num_cur <- num_b; den_cur <- den_b
            w <- w_b
            id <- id_b
            changed <- TRUE
            hist <- rbind(hist, data.frame(iter = it, ID = id, changed = "BATCH", note = paste(length(take_vars), "vars")))
            if (verbose) cat(sprintf("   iter %d: [BATCH] x%.2f on %d vars (ID: %.3f)\n", it, batch_factor, length(take_vars), id))
          }
        }
      }
    }
  }
  
  if (plot_progress && requireNamespace("ggplot2", quietly = TRUE)) {
    gp <- ggplot2::ggplot(hist, ggplot2::aes(iter, ID)) +
      ggplot2::geom_line() + ggplot2::geom_point() +
      ggplot2::labs(title = sprintf("Stochastic Gower Optim (N=%d, Batch=%d)", n_sub, N_SAMPLE_PER_ITER)) +
      ggplot2::theme_minimal()
    print(gp)
  }
  
  list(weights = w, history = hist, final_ID = id, idx_used = ix_sub_from_reps)
}

# --- Post-Optimization Selection ---

knee_triangle <- function(w) {
  y <- sort(as.numeric(w), decreasing = TRUE)
  n <- length(y)
  
  if (n < 3L) {
    return(list(k = n, thr = if (n) y[n] else NA_real_, curve = data.frame(i = seq_len(n), w = y, d = rep(0, n))))
  }
  
  x <- seq_len(n)
  num <- abs((y[n] - y[1]) * x - (n - 1) * y + n * y[1] - y[n])
  den <- sqrt((y[n] - y[1])^2 + (n - 1)^2)
  d <- num / den
  k <- which.max(d)
  
  list(k = k, thr = y[k], curve = data.frame(i = x, w = y, d = d))
}

survivors_from_weights <- function(w, 
                                   w_min = W_MIN,
                                   kmin = NULL,
                                   eps_ceil = 1e-4,
                                   eps_floor = 1e-12,
                                   make_plot = TRUE,
                                   plot_file = "FIG_weight_curve_knee.png") {
  
  if (is.null(names(w))) names(w) <- paste0("V", seq_along(w))
  w[] <- pmax(w_min, as.numeric(w))
  p <- length(w)
  
  cat(sprintf("[weights] p=%d | min=%.4f  q25=%.4f  med=%.4f  q75=%.4f  max=%.4f\n",
              p, min(w), as.numeric(quantile(w, .25)), median(w), as.numeric(quantile(w, .75)), max(w)))
  
  idx_ceil <- which(w >= 1 - eps_ceil)
  idx_tail <- which(w < 1 - eps_ceil & w > w_min + eps_floor)
  thr_tail <- NA_real_
  knee_obj <- NULL
  
  if (length(idx_tail) >= 3L) {
    knee_obj <- knee_triangle(w[idx_tail])
    thr_tail <- knee_obj$thr
  } else if (length(idx_tail) > 0L) {
    thr_tail <- max(w_min + 1e-6, median(w[idx_tail]))
  }
  
  S_ceil <- names(w)[idx_ceil]
  S_tail <- if (is.finite(thr_tail)) names(w)[w >= thr_tail & w < 1 - eps_ceil] else character(0)
  survivors <- union(S_ceil, S_tail)
  
  if (is.null(kmin)) kmin <- max(30L, ceiling(0.10 * p))
  if (length(survivors) < kmin) survivors <- names(sort(w, decreasing = TRUE))[seq_len(kmin)]
  
  cat(sprintf("Ceiling kept: %d | tail knee thr=%s | survivors: %d / %d\n", 
              length(S_ceil), ifelse(is.finite(thr_tail), sprintf("%.3f", thr_tail), "NA"), length(survivors), p))
  
  if (isTRUE(make_plot) && requireNamespace("ggplot2", quietly = TRUE)) {
    ord <- order(w, decreasing = TRUE)
    curve <- data.frame(i = seq_along(ord), w = as.numeric(w[ord]))
    pplt <- ggplot2::ggplot(curve, ggplot2::aes(i, w)) +
      ggplot2::geom_line() +
      ggplot2::labs(x = "rank (sorted ↓)", y = "weight", title = "Weight curve with knee (tail-only)") +
      ggplot2::theme_minimal()
    if (!is.null(knee_obj)) {
      thr_mark <- knee_obj$thr
      knee_row <- curve[which.min(abs(curve$w - thr_mark)), , drop = FALSE]
      pplt <- pplt + ggplot2::geom_point(data = knee_row, size = 2)
    }
    print(pplt)
    save_plot_gg(plot_file, pplt, width = 6, height = 4, dpi = 150)
  }
  
  list(survivors = survivors, thr_tail = thr_tail, w_clamped = w)
}

# --- PCA and Residuals Helpers ---

design_with_map <- function(X) {
  Xg <- as.data.frame(X, check.names = TRUE, stringsAsFactors = FALSE)
  if (!ncol(Xg)) stop("[design_with_map] input has 0 columns")
  
  # Filter constant cols
  keep <- vapply(Xg, function(v) length(unique(na.omit(v))) >= 1L, logical(1))
  if (!any(keep)) stop("[design_with_map] all columns are NA-only")
  Xg <- Xg[, keep, drop = FALSE]
  
  # Type standardization for modeling
  for (nm in names(Xg)) {
    v <- Xg[[nm]]
    if (is.ordered(v)) {
      Xg[[nm]] <- as.numeric(v)
      next
    }
    if (is.numeric(v) || is.integer(v)) next
    if (is.logical(v)) {
      Xg[[nm]] <- factor(v, levels = c(FALSE, TRUE))
      next
    }
    if (is.factor(v)) next
    if (is.matrix(v)) {
      Xg[[nm]] <- as.numeric(v)
      next
    }
    Xg[[nm]] <- factor(as.character(v))
  }
  
  fml <- as.formula(paste("~", paste(colnames(Xg), collapse = " + "), "-1"))
  tm <- terms(fml, data = Xg)
  MM <- model.matrix(tm, data = Xg)
  storage.mode(MM) <- "double"
  
  # Filter zero variance in design matrix
  ok <- apply(MM, 2L, function(col) {
    v <- stats::var(as.numeric(col), na.rm = TRUE)
    is.finite(v) && v > 1e-12
  })
  
  if (!any(ok)) stop("[design_with_map] all encoded columns were ~zero-variance")
  
  assign <- attr(MM, "assign")
  tl <- attr(tm, "term.labels")
  varmap <- tl[assign]
  MM <- MM[, ok, drop = FALSE]
  attr(MM, "varmap") <- varmap[ok]
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
  idx <- which(varmap == nm)
  if (!length(idx)) return(rep(NA_real_, nrow(Z)))
  
  if (length(idx) == 1L) return(as.numeric(E_scaled[, idx]))
  
  # Project E_scaled columns onto the item's original structure
  pc1 <- try(prcomp(Z[, idx, drop = FALSE], rank. = 1), silent = TRUE)
  if (inherits(pc1, "try-error")) return(rep(NA_real_, nrow(Z)))
  
  as.numeric(as.matrix(E_scaled[, idx, drop = FALSE]) %*% pc1$rotation[, 1])
}

r2_base_linear <- function(Base, v) {
  if (nrow(Base) != length(v)) return(NA_real_)
  fit <- try(lm(v ~ Base[, 1] + Base[, 2]), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  summary(fit)$r.squared
}

r2_residual_cv <- function(e, nb, folds = CV_FOLDS, seed = SEED_GLOBAL) {
  set.seed(seed)
  n <- length(e)
  if (n < 10 || var(e, na.rm = TRUE) == 0 || all(is.na(e))) return(NA_real_)
  
  fold_id <- sample(rep(1:folds, length.out = n))
  pred <- rep(NA_real_, n)
  
  for (f in 1:folds) {
    te <- which(fold_id == f)
    if (length(te) == 0) next
    
    # KNN imputation prediction on the test fold
    e_mask <- e
    e_mask[te] <- NA 
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
  res <- oof_R2_two_gams(v, BaseDF, y, K_target = 5, k_gam = 10, seed = 42)
  
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

df <- readr::read_delim(
  PSY_CSV,
  delim = ";",
  locale = readr::locale(decimal_mark = "."),
  progress = FALSE, show_col_types = FALSE
)

id_col <- if ("participant_id" %in% names(df)) "participant_id" else names(df)[1]
ids_all <- as.character(df[[id_col]])

# Clean columns
drop_pattern <- "(?i)diagnosis|SCID|NODIAG"
cols_to_drop <- grep(drop_pattern, names(df), value = TRUE)

X <- dplyr::select(df, -dplyr::all_of(c(id_col, cols_to_drop))) |>
  as.data.frame(stringsAsFactors = FALSE)

# Robust type coercion
is_small_int_scale <- function(v) {
  vn <- suppressWarnings(as.numeric(v))
  if (all(is.na(vn))) return(FALSE)
  u <- sort(unique(na.omit(vn)))
  k <- length(u)
  k >= 3 && k <= 7 && all(abs(u - round(u)) < 1e-8)
}

for (nm in names(X)) {
  v <- X[[nm]]
  if (is.character(v)) {
    vn <- suppressWarnings(as.numeric(v))
    if (!all(is.na(vn))) v <- vn
  }
  if (all(v %in% c(0, 1, NA))) {
    X[[nm]] <- factor(v, levels = c(0, 1), ordered = TRUE)
  } else if (is_small_int_scale(v)) {
    X[[nm]] <- factor(as.integer(round(as.numeric(v))), ordered = TRUE)
  } else if (is.numeric(v)) {
    X[[nm]] <- as.numeric(v)
  } else {
    X[[nm]] <- factor(v)
  }
}

keep <- stats::complete.cases(X) & !is.na(ids_all) & ids_all != ""
X <- X[keep, , drop = FALSE]
ids_all <- ids_all[keep]
rownames(X) <- make.unique(ids_all)

# --- Diagnoses Table Ingest ---
stopifnot(file.exists(DIAG_CSV))
diag_wide_full <- readr::read_delim(
  DIAG_CSV,
  delim = ";",
  col_types = readr::cols(),
  progress = FALSE, show_col_types = FALSE
) |>
  dplyr::mutate(participant_id = as.character(participant_id)) |>
  dplyr::mutate(dplyr::across(
    -participant_id,
    ~ {
      x <- suppressWarnings(as.integer(.))
      x[is.na(x)] <- 0L
      pmin(pmax(x, 0L), 1L)
    }
  ))

# Align diagnoses to X row order
ids_here <- rownames(X)
mm <- match(ids_here, as.character(diag_wide_full$participant_id))
if (anyNA(mm)) stop(sprintf("Diagnoses join failed for %d/%d rows.", sum(is.na(mm)), length(mm)))

DX_wide <- as.data.frame(diag_wide_full[mm, , drop = FALSE])
DX_wide$participant_id <- NULL

if (isTRUE(DX_DENY_NOS)) {
  nos_cols <- names(DX_wide)[tolower(names(DX_wide)) %in% "nodiag"]
  if (length(nos_cols)) DX_wide[nos_cols] <- NULL
}

cat(sprintf(
  "[ingest] X rows=%d, cols=%d | DX rows=%d, cols=%d (post NOS=%s)\n",
  nrow(X), ncol(X), nrow(DX_wide), ncol(DX_wide),
  if (isTRUE(DX_DENY_NOS)) "drop-nodiag" else "kept"
))

# ==============================================================================
# 3. Deduplication and Coreset Selection
# ==============================================================================

PX <- prep_X_for_gower(X, rare_prop = RARE_LEVEL_MIN_PROP, do_jitter = TRUE)
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

gr_all <- if (isTRUE(DO_DEDUP)) complete_groups(Dg, EPS_DEDUP) else seq_len(attr(Dg, "Size"))

# Identify medoids of groups
Dm_g <- as.matrix(Dg)
diag(Dm_g) <- 0
split_idx <- split(seq_len(nrow(Dm_g)), gr_all)
reps <- vapply(split_idx, function(ix) {
  ix[which.min(rowSums(Dm_g[ix, ix, drop = FALSE]))]
}, integer(1))
mult <- as.integer(lengths(split_idx))

Dg_rep <- stats::as.dist(as.matrix(Dg)[reps, reps, drop = FALSE])
core_idx_rep <- twonn_core_by_slope(Dg_rep)

cat(sprintf("[Dedup] eps=%.3f | reps=%d of %d | core_rep=%d\n", 
            EPS_DEDUP, length(reps), nrow(X), length(core_idx_rep)))

if (isTRUE(WRITE_DEDUP_CSV)) {
  mult_df <- data.frame(rep_row = reps, representative_id = 
                          rownames(Dm_g)[reps], multiplicity = mult)
  write_csv(mult_df, sprintf("near_duplicate_groups_complete_eps%g.csv", EPS_DEDUP))
}

# ==============================================================================
# 4. Gower Weight Optimization
# ==============================================================================

# Drop constant columns and optionally trim high-correlation pairs
drop_constant_cols <- function(X) {
  keep <- vapply(X, function(v) length(unique(na.omit(v))) >= 1L, logical(1))
  X[, keep, drop = FALSE]
}
X_pred <- drop_constant_cols(X)

if (DO_CORR_TRIM) {
  corr_trim <- function(X, thr = 0.95) {
    num <- X[, vapply(X, is.numeric, logical(1)), drop = FALSE]
    if (!ncol(num)) return(X)
    
    C <- suppressWarnings(cor(num, use = "pairwise.complete.obs"))
    drop <- character(0)
    for (i in seq_len(ncol(C) - 1)) {
      if (colnames(C)[i] %in% drop) next
      j <- which(abs(C[i, (i + 1):ncol(C)]) >= thr) + i
      drop <- union(drop, colnames(C)[j])
    }
    keep <- setdiff(colnames(X), drop)
    X[, keep, drop = FALSE]
  }
  X_pred <- corr_trim(X_pred, CORR_THRESH)
}
cat(sprintf("Start: X_pred has %d columns after constant-drop%s.\n\n", ncol(X_pred), if (DO_CORR_TRIM) " + corr-trim" else ""))

w_init <- setNames(rep(1, ncol(X)), colnames(X))
allow <- setNames(rep(FALSE, ncol(X)), colnames(X))
allow[colnames(X_pred)] <- TRUE

wopt <- optimise_gower_weights_constrained(X,
                                           init_weights = w_init,
                                           allow_update = allow,
                                           objective = "TwoNN_all",
                                           w_min = W_MIN,
                                           step_grid = W_STEP_GRID,
                                           batch_k = W_BATCH_K,
                                           batch_factor = W_BATCH_FACTOR,
                                           max_iter = W_MAX_ITERS,
                                           n_rows_sub = N_ROWS_SUB,
                                           ncores = NCORES_PAR,
                                           seed_jitter = SEED_JITTER,
                                           reps_idx = reps,
                                           core_idx_rep = core_idx_rep,
                                           verbose = TRUE,
                                           plot_progress = TRUE)

# Select survivor vars via knee on weight tail
sel <- survivors_from_weights(w = wopt$weights, w_min = W_MIN, kmin = max(30L, ceiling(0.10 * length(wopt$weights))), make_plot = TRUE)
w_full <- sel$w_clamped
survivors <- sel$survivors
X <- X[, survivors, drop = FALSE]
w_all <- w_full[survivors]
cat(sprintf("[TwoNN optimiser] Survivors: p=%d\n", ncol(X)))

# Recalculate Distances with final weights
PX <- prep_X_for_gower(X, rare_prop = RARE_LEVEL_MIN_PROP, do_jitter = TRUE)
Xg <- PX$X
type_g <- PX$type
w_use <- w_all[colnames(Xg)]

D_final <- cluster::daisy(Xg, metric = "gower", type = type_g, weights = w_use)
ID_all <- twonn_id_from_dist(D_final)
core_ix <- twonn_core_by_slope(D_final)
DmF <- as.matrix(D_final)
diag(DmF) <- Inf
ID_core <- twonn_id_from_dist(as.dist(DmF[core_ix, core_ix, drop = FALSE]))
ID_LB <- lb_mle_id(DmF[core_ix, core_ix, drop = FALSE], 5, 15)

cat(sprintf("[Constrained %s] TwoNN_all=%.2f | TwoNN_core=%.2f | LB_core=%.2f (n_core=%d, p_active=%d)\n",
            toupper(PREF_TARGET), ID_all, ID_core, ID_LB, length(core_ix), ncol(Xg)))

# ==============================================================================
# 5. PCA, Whitening, and Residuals
# ==============================================================================

# Expand design matrix
Xenc <- design_with_map(X)
varmap <- attr(Xenc, "varmap")
vars <- unique(varmap)

# Distribute Gower weights to encodings; sqrt-scale so PCA respects them
w_enc <- setNames(rep(1, ncol(Xenc)), colnames(Xenc))
alloc <- table(varmap)
for (nm in names(alloc)) {
  idx <- which(varmap == nm)
  wj <- w_all[nm]
  if (!is.finite(wj)) wj <- 1
  w_enc[idx] <- wj / length(idx)
}
Xenc_w <- sweep(Xenc, 2, sqrt(pmax(w_enc, 0)), "*")

m_star <- as.integer(M_STAR_FIXED)
Z <- scale(Xenc_w, center = TRUE, scale = TRUE)

# Robust PCA via PcaHubert (ROBPCA)
k_eff <- max(2L, min(m_star, nrow(Z) - 1L, ncol(Z)))
rp <- rrcov::PcaHubert(
  x        = Z,
  k        = k_eff,        # number of components to keep
  kmax     = k_eff,        # don't compute more than you need
  scale    = FALSE,        # Z is already standardised
  signflip = TRUE          # make loadings comparable to prcomp()
)
Base <- rp@scores[, 1:m_star, drop = FALSE]
colnames(Base) <- paste0("b", seq_len(ncol(Base)))

# --- Orientation Fix ---
# 1. Automate orientation via skewness (align heavy tail to positive)
# m1 <- apply(Base[, 1:2, drop = FALSE], 2, function(x) mean((x - mean(x))^3) / sd(x)^3)

#if (is.finite(m1[1]) && is.finite(m1[2]) && (m1[2] < m1[1])) {
  Base[, 1] <- -Base[, 1]
  rp@loadings[, 1] <- -rp@loadings[, 1]
#}

# 2. Hard override for dim 2 (ext up)
Base[, 2] <- -Base[, 2] 
rp@loadings[, 2] <- -rp@loadings[, 2]

# 3. Create snapshot for plotting
Base_A <- as.data.frame(Base[, 1:2, drop = FALSE])
colnames(Base_A) <- c("b1", "b2")

cat("[Orientation] Base coordinates flipped. Proceeding to Whitening...\n")

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

# --- Residual Extraction (Fibre) ---
make_strat_folds <- function(y, K, seed = SEED_GLOBAL) {
  set.seed(seed)
  y <- as.integer(y)
  n <- length(y)
  folds <- integer(n)
  idx0 <- which(y == 0)
  idx1 <- which(y == 1)
  f0 <- sample(rep(1:K, length.out = length(idx0)))
  f1 <- sample(rep(1:K, length.out = length(idx1)))
  folds[idx0] <- f0
  folds[idx1] <- f1
  folds
}

choose_K <- function(y, K_target = CV_FOLDS, min_per_class = 8) {
  y <- as.integer(y)
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  max(2, min(K_target, floor(n1 / min_per_class), floor(n0 / min_per_class)))
}

# Define Target for CV folds
stopifnot("NODIAG" %in% names(DX_wide))
y_use <- as.integer(DX_wide$NODIAG == 0L)  
stopifnot(length(y_use) == nrow(X))

cat(sprintf("[target] from NODIAG | prev=%.3f | n0=%d (NODIAG=1) | n1=%d (any dx)\n",
            mean(y_use), sum(y_use == 0L), sum(y_use == 1L)))

K_fold <- choose_K(y_use, K_target = CV_FOLDS, min_per_class = 8)
fold_id <- make_strat_folds(y_use, K = K_fold, seed = SEED_GLOBAL)

E <- residualise_foldsafe(Xenc_w, Base, folds = fold_id, k_gam = 6)
E_scaled <- scale(E, center = TRUE, scale = TRUE)

# XR represents the "Residuals" (high-dimensional, post-OOF)
XR <- E_scaled
if (any(!is.finite(XR))) XR[!is.finite(XR)] <- 0

cat(sprintf("[Residuals] XR matrix: %d rows × %d columns (post-OOF, scaled)\n", nrow(XR), ncol(XR)))

# ==============================================================================
# 6. Item Diagnostics (Base vs Fibre Roles)
# ==============================================================================

Z0_std <- scale(Xenc, center = TRUE, scale = TRUE)
vars_diag <- unique(varmap)

# Helper for KNN K selection
choose_k_nb <- function(e, nb_list, folds = CV_FOLDS, seed = SEED_GLOBAL) {
  r2s <- vapply(nb_list, function(nb) r2_residual_cv(e, nb, folds, seed), numeric(1))
  ix <- which.max(r2s)
  list(k = as.integer(sub("^k", "", names(nb_list)[ix])), R2_cv = r2s[ix])
}

RESIDUAL_PERM_B <- N_PERM
MIN_SD_ITEM <- 1e-6

cat(sprintf("[Item Roles] Starting diagnostics on %d items...\n", length(vars_diag)))

roles_rows <- lapply(vars_diag, function(nm) {
  # Score item (Original)
  v <- score_item_base(nm, Z0_std, varmap)
  if (all(is.na(v)) || stats::sd(v, na.rm = TRUE) < MIN_SD_ITEM) return(NULL)
  
  # Score residual (Fibre)
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

# ==============================================================================
# 7. Fibre-Only Decomposition & Clustering
# ==============================================================================

stopifnot(exists("E_scaled"), is.matrix(E_scaled), nrow(E_scaled) >= 3) 
Ef <- scale(E_scaled, center = TRUE, scale = TRUE)
Ef[!is.finite(Ef)] <- 0
n <- nrow(Ef)
pE <- ncol(Ef)
if (pE < 2L || n < 4L) stop("[Fibre-only] insufficient columns/rows in E.")

RESIDUAL_BASE_MAX <- min(6L, pE, n - 1L)

pc_f <- prcomp(Ef, rank. = max(2L, RESIDUAL_BASE_MAX))
Bprime_all <- pc_f$x[, 1:RESIDUAL_BASE_MAX, drop = FALSE]
colnames(Bprime_all) <- paste0("f", seq_len(ncol(Bprime_all)))

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
set.seed(SEED_PRED)
Kf <- max(2L, min(CV_FOLDS, n))
folds_f <- sample(rep(1:Kf, length.out = n))

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

cat(sprintf("[Fibre-only] ID(B'): TwoNN_all=%.2f | TwoNN_core=%.2f | LB_core=%.2f (n_core=%d)\n",
            ID_B_all, ID_B_core, ID_B_LB, length(core_B)))
cat(sprintf("[Fibre-only] ID(F'): TwoNN_all=%.2f | TwoNN_core=%.2f | LB_core=%.2f (n_core=%s)\n",
            ID_Fp_all, ID_Fp_core, ID_Fp_LB, if (exists("core_Fp")) length(core_Fp) else "NA"))

# Clustering in B' (Louvain on kNN graph)
knn_graph <- function(X, k = 15) {
  idx <- RANN::nn2(X, X, k = pmin(k + 1L, nrow(X)))$nn.idx[, -1L, drop = FALSE]
  i <- rep(seq_len(nrow(X)), each = ncol(idx))
  j <- as.vector(idx)
  g <- igraph::graph_from_edgelist(cbind(i, j), directed = FALSE)
  igraph::simplify(g)
}
gF <- knn_graph(Bprime, k = 15)
clF <- igraph::cluster_louvain(gF)$membership

saveRDS(list(
  Bprime = Bprime,
  Fprime = Fprime,
  m_f = m_f,
  ID_Bprime = c(all = ID_B_all, core = ID_B_core, LB = ID_B_LB),
  ID_Fprime = c(all = ID_Fp_all, core = ID_Fp_core, LB = ID_Fp_LB),
  clusters = clF
), file = "residual_only_summary.rds")

cat("[Fibre-only] wrote:", file.path(OUTPUTS_DIR, "residual_only_summary.rds"), "\n")

# ==============================================================================
# 8. Predictive Diagnostics (OOF & Interactions)
# ==============================================================================

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

vars_item_interact <- unique(varmap)
n_probe_count <- if (exists("N_TOP_PER_DX")) N_TOP_PER_DX else 50
vars_probe <- head(vars_item_interact, min(n_probe_count, length(vars_item_interact)))

# Pre-calculate scaled matrix for speed
if (!exists("Z0_std")) Z0_std <- scale(Xenc, center = TRUE, scale = TRUE)

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

oof_R2_two_gams <- function(v, Base, dx, K_target = 5, k_gam = 10, seed = SEED_GLOBAL) {
  n <- length(v)
  if (n != nrow(Base) || n != length(dx)) {
    return(c(R2_add = NA, R2_int = NA, p_like = NA, dR2 = NA))
  }
  
  K <- choose_K_dx(dx, K_target = K_target, min_per_class = 6L)
  if (K < 2) return(c(R2_add = NA, R2_int = NA, p_like = NA, dR2 = NA))
  
  set.seed(seed)
  fold_id <- sample(rep(1:K, length.out = n))
  lev_all <- levels(factor(dx))
  
  y_add <- rep(NA_real_, n)
  y_int <- rep(NA_real_, n)
  
  # Note: Force single-threaded GAM fitting here since we are inside parallel loop
  ctrl <- mgcv::gam.control(maxit = 100, nthreads = 1)
  
  for (k in sort(unique(fold_id))) {
    tr <- which(fold_id != k)
    te <- which(fold_id == k)
    if (length(tr) < 10 || length(te) == 0) next 
    
    dtr <- data.frame(v = v[tr], b1 = Base[tr, 1], b2 = Base[tr, 2], dx = factor(dx[tr], levels = lev_all))
    dte <- data.frame(b1 = Base[te, 1], b2 = Base[te, 2], dx = factor(dx[te], levels = lev_all))
    
    # 1. Additive Model
    f_add <- try(mgcv::gam(v ~ s(b1, b2, k = k_gam, bs = "tp", m = 2),
                           data = dtr, method = "REML", gamma = 1.4, control = ctrl), silent = TRUE)
    
    # 2. Interaction Model
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
    
    if (any(!is.finite(y_add[te]))) y_add[te][!is.finite(y_add[te])] <- mu
    if (any(!is.finite(y_int[te]))) y_int[te][!is.finite(y_int[te])] <- mu
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
      choose_K_dx = choose_K_dx, interact_worker = interact_worker, p = p
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

# ==============================================================================
# 9. Outputs and Session
# ==============================================================================

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

f_means <- sapply(seq_len(ncol(Bprime)), function(j) mean_by_cluster(Bprime[, j], cl_fac))
if (!is.null(dim(f_means))) {
  colnames(f_means) <- paste0("f", seq_len(ncol(Bprime)), "_mean")
} else {
  f_means <- cbind(`f1_mean` = as.numeric(f_means))
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
  
  if (inherits(chi, "try-error")) {
    Eexp_used <- outer(rowSums(tab), colSums(tab), function(r, c) r * c / max(sum(tab), 1))
    Z_enrich <- (tab - Eexp_used) / sqrt(pmax(Eexp_used, 1e-9))
    pmat <- 2 * pnorm(-abs(Z_enrich))
  } else {
    Z_enrich <- chi$stdres
    Eexp_used <- chi$expected
    pmat <- 2 * pnorm(-abs(Z_enrich))
  }
  
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
  write_csv(data.frame(cluster = rownames(Z_enrich), as.data.frame(Z_enrich), check.names = FALSE), "dx_cluster_enrichment_Z.csv")
  write_csv(data.frame(cluster = rownames(pmat), as.data.frame(pmat), check.names = FALSE), "dx_cluster_enrichment_p.csv")
  write_csv(data.frame(cluster = rownames(q_by_dx), as.data.frame(q_by_dx), check.names = FALSE), "dx_cluster_enrichment_q_by_dx.csv")
  write_csv(data.frame(cluster = rownames(q_global_mat), as.data.frame(q_global_mat), check.names = FALSE), "dx_cluster_enrichment_q_global.csv")
  write_csv(data.frame(cluster = rownames(prev_in_cluster), as.data.frame(prev_in_cluster), check.names = FALSE), "dx_cluster_prev_in_cluster.csv")
  
  # Export Long table
  grid <- expand.grid(
    cluster = rownames(tab),
    dx = colnames(tab),
    stringsAsFactors = FALSE
  )
  idx_r <- match(grid$cluster, rownames(tab))
  idx_c <- match(grid$dx, colnames(tab))
  
  grid$n_cluster <- ncl[grid$cluster]
  grid$cases_cluster <- tab[cbind(idx_r, idx_c)]
  grid$prev_cluster <- prev_in_cluster[cbind(idx_r, idx_c)]
  grid$prev_overall <- prev_overall[grid$dx]
  grid$Z <- Z_enrich[cbind(idx_r, idx_c)]
  grid$p <- pmat[cbind(idx_r, idx_c)]
  grid$q_by_dx <- q_by_dx[cbind(idx_r, idx_c)]
  grid$q_global <- q_global_mat[cbind(idx_r, idx_c)]
  grid$expected_cases <- Eexp_used[cbind(idx_r, idx_c)]
  grid$enriched <- (grid$Z > 0) & (grid$q_by_dx < SIG_Q)
  
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
  
  list(mu = mu, S_half = S_half, S_half_inv = S_hi, s = s, fwd = fwd, inv = inv)
}

make_unitdisk_square <- function(nu) {
  u1 <- seq(-1, 1, length.out = nu)
  u2 <- seq(-1, 1, length.out = nu)
  g <- expand.grid(u1 = u1, u2 = u2)
  mask <- with(g, sqrt(u1^2 + u2^2) <= 1 + 1e-9)
  list(u1 = u1, u2 = u2, grid = g, mask_vec = mask)
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

build_geometry <- function(Base_A, GRID_N_B = 600, GRID_N_U = 600) {
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
varmap <- attr(Z, "varmap")
Z_A    <- Z

# Geometry & grids
geom <- build_geometry(Base_A)

# Export Embeddings
emb_base <- data.frame(
  participant_id = ids_base,
  b1 = Base[, 1], b2 = Base[, 2],
  cluster = as.integer(clF),
  target = as.integer(y_use),
  stringsAsFactors = FALSE
)
write_csv(emb_base, "embedding_base_b1b2.csv")

emb_residual <- data.frame(
  participant_id = ids_base,
  Bprime,
  cluster = as.integer(clF),
  stringsAsFactors = FALSE
)
write_csv(emb_residual, "embedding_residual_Bprime.csv")

saveRDS(list(
  participant_id = ids_base, 
  XR = E_scaled, 
  Fprime = Fprime
),
file = "Fprime_matrix.rds"
)

# Export Weights and Maps
w_tbl <- data.frame(
  var = names(w_full),
  weight = as.numeric(w_full),
  selected = names(w_full) %in% survivors,
  stringsAsFactors = FALSE
)
write_csv(w_tbl, "gower_weights_optimised.csv")

enc_map <- data.frame(
  mm_col = colnames(Xenc),
  source_var = as.character(varmap),
  weight_share = as.numeric(w_enc[colnames(Xenc)]),
  stringsAsFactors = FALSE
)
write_csv(enc_map, "encoding_map_and_weight_share.csv")

# Session Info
zz <- file("sessionInfo.txt", open = "wt")
sink(zz)
on.exit(
  {
    sink()
    close(zz)
  },
  add = TRUE
)
print(sessionInfo())
cat("\n\nConfig snapshot:\n")
print(str(cfg, max.level = 1))

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
    { scico::scale_colour_scico(palette = cfg$palette$name, direction = 1) } +
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

# Participant density plot
p_dens <- density_plots(geom)
save_plot_gg("FIG_dens_unitsquare_scatter", p_dens, width = 8.0, height = 7.0)

# Direction plot
p_dir <- direction_wheel_plot(geom)
if (!is.null(p_dir)) {
  save_plot_gg("FIG_uv_direction_density_HCL_smooth", p_dir, width = 8.0, height = 7.0)
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

pearson_r <- function(a, b) {
  if (length(a) != length(b)) return(NA_real_)
  suppressWarnings(cor(as.numeric(a), as.numeric(b), use = "complete.obs", method = "pearson"))
}

build_biplot_data <- function(Z_A, varmap, Base_A, U) {
  B1 <- Base_A[, 1]; B2 <- Base_A[, 2]
  items <- unique(varmap)
  Rtab <- dplyr::bind_rows(lapply(items, function(nm) {
    v <- score_item_1d(nm, Z_A, varmap)
    if (!any(is.finite(v))) return(NULL)
    data.frame(
      item = nm,
      r_b1 = pearson_r(v, B1), r_b2 = pearson_r(v, B2),
      r_u1 = pearson_r(v, U$u1), r_u2 = pearson_r(v, U$u2),
      stringsAsFactors = FALSE
    )
  }))
  if (!nrow(Rtab)) stop("[biplot] No item correlations computed.")
  Rtab |>
    dplyr::mutate(
      across(c(r_b1, r_b2, r_u1, r_u2), ~ suppressWarnings(as.numeric(.))),
      r_b1 = ifelse(is.finite(r_b1), r_b1, 0),
      r_b2 = ifelse(is.finite(r_b2), r_b2, 0),
      r_u1 = ifelse(is.finite(r_u1), r_u1, 0),
      r_u2 = ifelse(is.finite(r_u2), r_u2, 0),
      mag_r_base = sqrt(r_b1^2 + r_b2^2),
      mag_r_disk = sqrt(r_u1^2 + r_u2^2)
    )
}

plot_biplots <- function(Rtab, Base_A, U) {
  suppressPackageStartupMessages(requireNamespace("ggrepel", quietly = TRUE))
  
  B1n <- suppressWarnings(as.numeric(Base_A[, 1]))
  B2n <- suppressWarnings(as.numeric(Base_A[, 2]))
  B1n[!is.finite(B1n)] <- mean(B1n, na.rm = TRUE)
  B2n[!is.finite(B2n)] <- mean(B2n, na.rm = TRUE)
  
  Hidx <- grDevices::chull(B1n, B2n)
  H <- data.frame(b1 = B1n[Hidx], b2 = B2n[Hidx])
  
  cx <- mean(B1n); cy <- mean(B2n)
  Rscale <- 0.80 * min(diff(range(B1n)), diff(range(B2n)))
  n_arrows <- min(30L, nrow(Rtab))
  
  S_base <- Rtab |>
    dplyr::arrange(dplyr::desc(mag_r_base)) |>
    dplyr::slice_head(n = n_arrows) |>
    dplyr::mutate(
      x0 = cx, y0 = cy,
      x1 = cx + Rscale * r_b1,
      y1 = cy + Rscale * r_b2
    )
  
  p_base <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = H, ggplot2::aes(b1, b2),
                          fill = NA, colour = "black", linewidth = 0.4) +
    ggplot2::geom_segment(
      data = S_base,
      ggplot2::aes(x = x0, y = y0, xend = x1, yend = y1),
      linewidth = 0.8, colour = "firebrick",
      arrow = grid::arrow(length = grid::unit(0.14, "cm"))
    ) +
    ggrepel::geom_label_repel(
      data = transform(S_base, xlab = x1, ylab = y1),
      ggplot2::aes(x = xlab, y = ylab, label = item),
      size = 3.1, max.overlaps = Inf, label.size = 0,
      label.padding = grid::unit(0.10, "lines")
    ) +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "b1", y = "b2") +
    theme_pub(12)
  
  draw_disk_outline <- function() {
    th <- seq(0, 2 * pi, length.out = 361)
    data.frame(x = cos(th), y = sin(th))
  }
  
  Rdisk <- 0.85
  S_disk <- Rtab |>
    dplyr::arrange(dplyr::desc(mag_r_disk)) |>
    dplyr::slice_head(n = n_arrows) |>
    dplyr::transmute(item, u0 = 0, v0 = 0,
                     u1 = as.numeric(Rdisk * r_u1),
                     v1 = as.numeric(Rdisk * r_u2))
  
  p_disk <- ggplot2::ggplot() +
    ggplot2::geom_path(data = draw_disk_outline(), ggplot2::aes(x, y)) +
    ggplot2::geom_segment(
      data = S_disk,
      ggplot2::aes(x = u0, y = v0, xend = u1, yend = v1),
      linewidth = 0.8, colour = "firebrick",
      arrow = grid::arrow(length = grid::unit(0.14, "cm"))
    ) +
    ggrepel::geom_text_repel(
      data = transform(S_disk, xlab = u1, ylab = v1),
      ggplot2::aes(x = xlab, y = ylab, label = item),
      size = 3.1, max.overlaps = Inf
    ) +
    ggplot2::coord_equal(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
    ggplot2::labs(x = "u1 (whitened b1,b2)", y = "u2") +
    theme_pub(12)
  
  list(p_base = p_base, p_disk = p_disk)
}

# Biplots
if (!is.null(varmap)) {
  Rtab <- build_biplot_data(Z_A, varmap, Base_A, geom$U)
  write_csv(Rtab, file.path(OUTPUTS_DIR, "items_vs_base_and_unitdisk_correlations.csv"))
  
  BIP <- plot_biplots(Rtab, Base_A, geom$U)
  save_plot_gg("FIG_biplot_items_BASE", BIP$p_base, width = 8.0, height = 7.0)
  save_plot_gg("FIG_biplot_items_UNITDISK", BIP$p_disk, width = 8.0, height = 7.0)
} else {
  msgf("[biplot] varmap missing; skipping biplots.")
}