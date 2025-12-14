# dimension_contplot.R
# Behavioural mapping with GAMs.
#
# Dependencies:
#   Expects 'setup.R' to be sourced previously (provides cfg, palette helpers, IO wrappers, futures).
#   Expects input data frames: Base_A (n x 2), fold_id (optional), and geometry object 'geom'.

options(future.debug = FALSE)
progressr::handlers(global = TRUE)

# ==============================================================================
# 0. Initialisation & Mode Checking
# ==============================================================================

target_name <- basename(behaviour_csv)

is_dx_mode <- identical(target_name, "wide_diagnoses.csv") || identical(OUT_SUBDIR, "wide_diagnoses")
is_cl_mode <- identical(target_name, "cluster_membership_all_participants.csv") || identical(OUT_SUBDIR, "clusters")
is_first_run <- FALSE

if (is_dx_mode) {
    message("Running diagnoses mode.")
    behaviour_csv <- "data/wide_diagnoses.csv"
    OUT_SUBDIR <- "wide_diagnoses"
  } else if (!is_dx_mode && !file.exists(file.path(OUTPUTS_DIR, "wide_diagnoses", "dx_fields.rds"))){
    message("Diagnoses mode has not run yet. Running now.")
    is_dx_mode <- TRUE
    is_first_run <- TRUE
    restore_csv <- behaviour_csv
    restore_subdir <- OUT_SUBDIR
    behaviour_csv <- "data/wide_diagnoses.csv"
    OUT_SUBDIR <- "wide_diagnoses"
  } else if (is_cl_mode) {
    message("Cluster mode detected. Executing cluster membership mapping.")
    behaviour_csv <- "out/cluster_membership_all_participants.csv"
    OUT_SUBDIR <- "clusters"
}

msgf <- function(...) message(sprintf(...))
`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

# ==============================================================================
# 1. Palette Helpers
# ==============================================================================

# Adaptive colour scales supporting scico, paletteer, and manual definitions from global config.
scale_prob_colour <- function(limits = NULL, name = NULL) {
  eng <- tolower(cfg$palette$engine %||% "scico")
  dir <- cfg$palette$direction %||% 1
  pal <- cfg$palette$name %||% "lipari"

  if (eng == "scico" && requireNamespace("scico", quietly = TRUE)) {
    ggplot2::scale_colour_gradientn(
      colours = scico::scico(256, palette = pal, direction = dir),
      limits = limits, oob = scales::squish, name = name,
      na.value = "transparent"
    )
  } else if (eng == "paletteer" && requireNamespace("paletteer", quietly = TRUE)) {
    cols <- as.character(paletteer::paletteer_c(pal, n = 256, direction = dir))
    ggplot2::scale_colour_gradientn(
      colours = cols, limits = limits, oob = scales::squish, name = name,
      na.value = "transparent"
    )
  } else if (eng == "manual" && !is.null(cfg$palette$colours)) {
    ggplot2::scale_colour_gradientn(
      colours = cfg$palette$colours, limits = limits, name = name,
      na.value = "transparent"
    )
  } else {
    ggplot2::scale_colour_gradient(
      low = "#f7fbff", high = "#08306b", limits = limits, name = name,
      na.value = "transparent"
    )
  }
}

scale_prob_fill <- function(limits = NULL, name = NULL) {
  eng <- tolower(cfg$palette$engine %||% "scico")
  dir <- cfg$palette$direction %||% 1
  pal <- cfg$palette$name %||% "lipari"

  if (eng == "scico" && requireNamespace("scico", quietly = TRUE)) {
    ggplot2::scale_fill_gradientn(
      colours = scico::scico(256, palette = pal, direction = dir),
      limits = limits, oob = scales::squish, name = name,
      na.value = "transparent"
    )
  } else if (eng == "paletteer" && requireNamespace("paletteer", quietly = TRUE)) {
    cols <- as.character(paletteer::paletteer_c(pal, n = 256, direction = dir))
    ggplot2::scale_fill_gradientn(
      colours = cols, limits = limits, oob = scales::squish, name = name,
      na.value = "transparent"
    )
  } else if (eng == "manual" && !is.null(cfg$palette$colours)) {
    ggplot2::scale_fill_gradientn(
      colours = cfg$palette$colours, limits = limits, name = name,
      na.value = "transparent"
    )
  } else {
    ggplot2::scale_fill_gradient(
      low = "#f7fbff", high = "#08306b", limits = limits, name = name,
      na.value = "transparent"
    )
  }
}

# ==============================================================================
# 2. Data Preparation
# ==============================================================================

# Robust parsing for mixed-type CSV columns (handles percentage signs, commas, etc.).
char_to_num <- function(v) {
  if (is.numeric(v)) {
    return(v)
  }
  if (is.factor(v) || is.logical(v)) v <- as.character(v)
  if (!is.character(v)) {
    return(v)
  }

  s <- gsub("\\s+", "", v)
  s <- sub("%$", "", s)
  v2 <- suppressWarnings(as.numeric(gsub(",", ".", s)))

  non_blank <- !is.na(s) & nzchar(s)
  prop_num <- if (any(non_blank)) sum(is.finite(v2[non_blank])) / sum(non_blank) else 0
  if (prop_num >= 0.80) v2 else v
}

B <- Base_A
pid <- rownames(B)
base_dt <- data.table::data.table(
  participant_id = as.character(pid),
  b1 = as.numeric(B[, 1]), b2 = as.numeric(B[, 2])
)

DT <- data.table::fread(
  behaviour_csv,
  na.strings = c("", "NA", "N/A", "NaN", "nan", "null", "NULL", ".", "-"),
  strip.white = TRUE
)
DT <- as.data.table(lapply(DT, char_to_num))
stopifnot("participant_id" %in% names(DT))
DT[, participant_id := as.character(participant_id)]

if (is_cl_mode) {
  if (!"cluster" %in% names(DT)) stop("Cluster file must have a 'cluster' column.")
  
  # Expand 'cluster' column (int) into binary columns C1, C2, etc.
  clusters_found <- sort(unique(DT$cluster[DT$cluster != 0]))
  msgf("[cl_mode] Expanding %d clusters into binary targets...", length(clusters_found))
  
  for (cid in clusters_found) {
    col_name <- paste0("C", cid)
    DT[, (col_name) := as.integer(cluster == cid)]
  }
  # Remove original cluster column to avoid processing it as continuous
  DT[, cluster := NULL]
}

DX <- merge(base_dt, DT, by = "participant_id", all = FALSE)

# Load residuals (XR) for stacking analysis, ensuring strict alignment with DX.
XR <- NULL
if (file.exists(file.path(OUTPUTS_DIR, "Fprime_matrix.rds"))) {
  tmp_res <- readRDS("Fprime_matrix.rds")
  if (is.list(tmp_res)) {
    mat_to_use <- if ("XR" %in% names(tmp_res)) tmp_res$XR else tmp_res$Fprime
    ix <- match(DX$participant_id, tmp_res$participant_id)

    if (anyNA(ix)) {
      warning("[setup] Missing participants in residuals file. Filtering DX to intersection.")
      keep_rows <- !is.na(ix)
      DX <- DX[keep_rows, ]
      ix <- match(DX$participant_id, tmp_res$participant_id)

      if (exists("geom") && !is.null(geom$U)) {
        U_full <- geom$U
        if (!is.null(rownames(U_full)) && all(DX$participant_id %in% rownames(U_full))) {
          U_DX <- U_full[DX$participant_id, , drop = FALSE]
        } else {
          U_DX <- NULL
        }
      }
    }

    if (nrow(DX) > 0 && !anyNA(ix)) {
      XR <- mat_to_use[ix, , drop = FALSE]
      message("[setup] Loaded Residuals (XR) for stacking (n=", nrow(DX), ").")
    } else {
      warning("[setup] No participant intersection found. Stacking disabled.")
      XR <- NULL
    }
  } else {
    # Handle raw matrix input
    XR <- tmp_res[match(DX$participant_id, rownames(tmp_res)), , drop = FALSE]
  }
}

# Align whitened geometry (U) to DX.
U_DX <- NULL
if (!is.null(geom$U)) {
  U_all <- geom$U
  if (!is.null(rownames(U_all)) &&
    "participant_id" %in% names(DX) &&
    all(DX$participant_id %in% rownames(U_all))) {
    U_DX <- U_all[DX$participant_id, , drop = FALSE]
  }
}

folds_all <- if (exists("fold_id")) fold_id[match(DX$participant_id, pid)] else rep(1L, nrow(DX))
K <- length(unique(folds_all))

in_whitened_geom <- function(U) {
  if ("insq" %in% names(U)) {
    with(U, insq & is.finite(u1) & is.finite(u2))
  } else if ("rin" %in% names(U)) {
    with(U, rin & is.finite(u1) & is.finite(u2))
  } else {
    with(U, is.finite(u1) & is.finite(u2) & (u1^2 + u2^2 <= 1 + 1e-9))
  }
}

# ==============================================================================
# 3. Interactive Model Selection
# ==============================================================================

# logic: Allows manual override of distribution families via console or CSV.
# Defaults are inferred from data properties (integer, range, contiguity).
get_model_overrides <- function(
  vars,
  DX,
  save_path = file.path(OUTPUTS_DIR, OUT_SUBDIR, "model_overrides.csv")
) {
  
  if (is_cl_mode) {
    mp <- setNames(vector("list", length(vars)), vars)
    for (v in vars) {
      mp[[v]] <- list(choice = "binomial", trials = NA_integer_, trials_var = NA_character_)
    }
    return(mp)
  }
  
  choices <- c("auto", "gauss_z", "gauss_log1p", "poisson", "nb", "binomial", "beta", "ordinal", "nominal", "binomial_trials")
  make_default_entry <- function() list(choice = "auto", trials = NA_integer_, trials_var = NA_character_)

  if (file.exists(save_path)) {
    ov_raw <- tryCatch(data.table::fread(save_path, colClasses = "character"), error = function(e) NULL)
    if (!is.null(ov_raw) && all(c("var", "choice") %in% names(ov_raw))) {
      if (!"trials" %in% names(ov_raw)) ov_raw[, trials := NA_character_]
      if (!"trials_var" %in% names(ov_raw)) ov_raw[, trials_var := NA_character_]

      map <- setNames(vector("list", nrow(ov_raw)), ov_raw$var)
      for (i in seq_len(nrow(ov_raw))) {
        v <- ov_raw$var[i]
        ch <- ov_raw$choice[i]
        tr <- ov_raw$trials[i]
        tv <- ov_raw$trials_var[i]

        e <- make_default_entry()
        e$choice <- ch
        if (identical(ch, "binomial_trials")) {
          if (!is.na(tv) && nzchar(tv) && tv %in% names(DX)) {
            e$trials_var <- tv
          } else if (!is.na(tr) && nzchar(tr)) {
            if (grepl("^[0-9]+$", tr)) {
              e$trials <- as.integer(tr)
            } else if (tr %in% names(DX)) e$trials_var <- tr
          }
        }
        map[[v]] <- e
      }
      for (m in setdiff(vars, names(map))) map[[m]] <- make_default_entry()
      return(map)
    }
  }

  if (!interactive()) {
    mp <- setNames(vector("list", length(vars)), vars)
    for (v in vars) mp[[v]] <- make_default_entry()
    return(mp)
  }

  message("\nModel chooser (console). 0 = accept suggested.\n")

  suggest_choice <- function(v) {
    y <- suppressWarnings(as.numeric(DX[[v]]))
    y[!is.finite(y)] <- NA_real_
    u <- sort(unique(y[is.finite(y)]))
    K <- length(u)
    is_int <- all(is.na(y) | abs(y - round(y)) < 1e-12)
    contig <- is_int && K >= 3 && K <= 12 && all(diff(u) == 1)

    if (grepl("(?:^|_)(rt|latency|lat|time|duration|sec|secs|seconds?|ms|millis)(?:$|_)", v, TRUE)) {
      return("gauss_log1p")
    }
    if (grepl("score", v, TRUE)) {
      return("binomial_trials")
    }
    if (contig) {
      return("ordinal")
    }
    "auto"
  }

  ask_fixed_trials <- function(v) {
    repeat {
      ans <- readline(prompt = paste0("  Fixed number of trials for ", v, " (integer > 1): "))
      val <- suppressWarnings(as.integer(ans))
      if (is.finite(val) && val > 1) {
        return(val)
      }
      cat("    Invalid. Please enter an integer > 1.\n")
    }
  }

  ask_trials_var <- function(v, DX_cols) {
    cand <- DX_cols[
      vapply(DX[, DX_cols, with = FALSE], function(col) {
        if (!is.numeric(col)) {
          return(FALSE)
        }
        all(is.na(col) | (is.finite(col) & col >= 0 & abs(col - round(col)) < 1e-9))
      }, logical(1))
    ]
    if (!length(cand)) {
      cat("  No valid count columns found. Fallback to fixed trials.\n")
      return(NA_character_)
    }
    cat("  Pick denominator column for ", v, " (per-row #trials):\n", sep = "")
    menu_choices <- c("None / use fixed integer", cand)
    idx <- utils::menu(menu_choices, title = "  Denominator variable:")
    if (idx <= 1) NA_character_ else cand[idx - 1L]
  }

  choose_for_var <- function(v) {
    sugg <- suggest_choice(v)
    cat("\n", v, "\n", sep = "")
    cat("  Suggested: ", sugg, "\n", sep = "")
    idx <- utils::menu(choices, title = "Pick model (0 = suggested):", graphics = FALSE)
    pick <- if (idx == 0) sugg else choices[idx]

    out <- make_default_entry()
    out$choice <- pick

    if (identical(pick, "binomial_trials")) {
      use_var <- utils::menu(
        c("Use per-row denominator (e.g., rk_kresponse)", "Use fixed integer"),
        title = "binomial_trials config:"
      )
      if (use_var == 1) {
        tv <- ask_trials_var(v, names(DX))
        if (is.character(tv) && !is.na(tv) && nzchar(tv)) {
          out$trials_var <- tv
        } else {
          out$trials <- ask_fixed_trials(v)
        }
      } else {
        out$trials <- ask_fixed_trials(v)
      }
    }
    out
  }

  out <- setNames(vector("list", length(vars)), vars)
  for (v in vars) out[[v]] <- choose_for_var(v)

  dir.create(dirname(save_path), recursive = TRUE, showWarnings = FALSE)
  data.table::fwrite(
    data.table::data.table(
      var = names(out),
      choice = vapply(out, function(e) e$choice, ""),
      trials = vapply(out, function(e) {
        tr <- e$trials
        if (is.null(tr) || is.na(tr)) NA_character_ else as.character(tr)
      }, character(1)),
      trials_var = vapply(out, function(e) {
        tv <- e$trials_var
        if (is.null(tv) || is.na(tv) || !nzchar(tv)) NA_character_ else tv
      }, character(1))
    ),
    save_path
  )
  message("Saved overrides to: ", normalizePath(save_path, mustWork = FALSE))
  out
}

# ==============================================================================
# 4. Family Resolution & Heuristics
# ==============================================================================

# Helper to handle varying argument names in mgcv::ocat
ocat_family <- function(K) {
  fml <- names(formals(mgcv::ocat))
  if ("K" %in% fml) mgcv::ocat(K = K) else if ("R" %in% fml) mgcv::ocat(R = K) else if ("k" %in% fml) mgcv::ocat(k = K) else mgcv::ocat(K)
}

# Maps override choices to specific GLM/GAM families and link functions.
family_from_choice <- function(choice, y, var_name = "", trials = NA_integer_) {
  y <- suppressWarnings(as.numeric(y))
  y[!is.finite(y)] <- NA_real_
  rtn <- function(y, fam, label, fwd = identity, inv = identity, extra = list()) {
    c(list(y = y, family = fam, label = label, fwd = fwd, inv = inv), extra)
  }

  if (choice == "gauss_z") {
    mu0 <- mean(y, na.rm = TRUE)
    sd0 <- stats::sd(y, na.rm = TRUE)
    if (!is.finite(sd0) || sd0 == 0) sd0 <- 1
    return(rtn(
      (y - mu0) / sd0, gaussian(), "(z/gauss)",
      function(x) (x - mu0) / sd0, function(mu) mu * sd0 + mu0
    ))
  }
  if (choice == "gauss_log1p") {
    return(rtn(
      log1p(pmax(y, 0)), gaussian(), "(log1p/gauss, manual)",
      function(x) log1p(pmax(x, 0)), function(mu) expm1(mu)
    ))
  }
  if (choice == "poisson") {
    return(rtn(y, poisson("log"), "(poisson)"))
  }
  if (choice == "nb") {
    fam <- try(mgcv::nb(), silent = TRUE)
    if (inherits(fam, "try-error")) fam <- poisson("log")
    return(rtn(y, fam, paste0("(", if (inherits(fam, "family")) fam$family else "nb", ")"), extra = list(is_nb = TRUE)))
  }
  if (choice == "binomial") {
    return(rtn(as.integer(y > 0 & is.finite(y)), binomial("logit"), "(binom/logit, manual)"))
  }
  if (choice == "beta") {
    eps <- 1e-6
    yb <- pmin(pmax(y, eps), 1 - eps)
    fam <- try(mgcv::betar("logit"), silent = TRUE)
    if (inherits(fam, "try-error")) fam <- quasibinomial("logit")
    return(rtn(yb, fam, "(beta/quasi, manual)"))
  }
  if (choice == "ordinal") {
    u <- sort(unique(y[is.finite(y)]))
    K <- length(u)
    lev_map <- setNames(seq_len(K), u)
    y_ord <- ifelse(is.finite(y), lev_map[as.character(y)], NA_integer_)
    fam <- ocat_family(K)
    return(rtn(y_ord, fam, sprintf("(ocat K=%d, ordinal)", K),
      fwd = function(x) lev_map[as.character(x)],
      inv = function(k) as.numeric(names(lev_map))[pmax(1, pmin(K, as.integer(round(k))))],
      extra = list(levels = u, is_ordinal = TRUE)
    ))
  }
  if (choice == "nominal") {
    u <- sort(unique(y[is.finite(y)]))
    return(rtn(y, binomial("logit"), "(nominal, one-vs-rest)", extra = list(levels = u, is_nominal = TRUE)))
  }
  if (choice == "binomial_trials") {
    if (!is.finite(trials) || trials <= 1L) stop("binomial_trials selected but 'trials' not provided (>1).")
    y_cnt <- suppressWarnings(as.integer(round(pmin(pmax(y, 0), trials))))
    return(list(
      y = y_cnt, family = binomial("logit"), label = sprintf("(binomial trials=%d)", as.integer(trials)),
      fwd = identity, inv = identity, trials = as.integer(trials), two_col = TRUE
    ))
  }
  NULL
}

# Auto-detection of family based on data distribution (sparse, bounded, integer, etc.).
choose_family <- function(y, var_name = "") {
  var_name <- tolower(var_name)
  y <- suppressWarnings(as.numeric(y))
  y[!is.finite(y)] <- NA_real_
  if (!any(is.finite(y))) {
    return(list(y = y, family = gaussian(), label = "(no data)", fwd = identity, inv = identity))
  }

  in01 <- all(is.na(y) | (y >= 0 & y <= 1))
  is_int <- all(is.na(y) | (abs(y - round(y)) < 1e-12))
  zrate <- mean(y == 0, na.rm = TRUE)
  q <- stats::quantile(y, c(.05, .5, .95), na.rm = TRUE, names = FALSE)
  rtn <- function(y, fam, label, fwd = identity, inv = identity, extra = list()) c(list(y = y, family = fam, label = label, fwd = fwd, inv = inv), extra)

  # Heuristics for RT/Durations
  if (grepl("(?:^|_)(rt|latency|lat|time|duration|sec|secs|seconds?|ms|millis)(?:$|_)", var_name, TRUE)) {
    yp <- pmax(y, 0)
    return(rtn(log1p(yp), gaussian(), "(log1p/gauss, duration)", log1p, function(mu) pmax(expm1(mu), 0)))
  }

  if (grepl("explosion", var_name)) {
    if (is_int && (zrate > 0.5)) {
      return(rtn(log1p(y), gaussian(), "(log1p/gauss, zero-heavy)", log1p, function(mu) pmax(expm1(mu), 0)))
    } else {
      fam <- try(mgcv::nb(), silent = TRUE)
      if (inherits(fam, "try-error")) fam <- poisson("log")
      return(rtn(y, fam, paste0("(", if (inherits(fam, "family")) fam$family else "nb", ")")))
    }
  }

  if (grepl("coefvar", var_name) && in01 && !is_int) {
    eps <- 1e-6
    yb <- pmin(pmax(y, eps), 1 - eps)
    fam <- try(mgcv::betar("logit"), silent = TRUE)
    if (inherits(fam, "try-error")) fam <- quasibinomial("logit")
    return(rtn(yb, fam, "(beta/quasi)"))
  }

  if (grepl("ratio", var_name)) {
    return(rtn(log1p(pmax(y, 0)), gaussian(), "(log1p/gauss, ratio)", function(x) log1p(pmax(x, 0)), function(mu) pmax(expm1(mu), 0)))
  }

  if (all(is.na(y) | y %in% c(0, 1))) {
    return(rtn(y, binomial("logit"), "(binom/logit)"))
  }

  if (in01 && !is_int) {
    eps <- 1e-6
    yb <- pmin(pmax(y, eps), 1 - eps)
    fam <- try(mgcv::betar("logit"), silent = TRUE)
    if (inherits(fam, "try-error")) fam <- quasibinomial("logit")
    return(rtn(yb, fam, "(beta/quasi)"))
  }

  # Heuristics for Counts/Integers
  if (is_int) {
    u <- sort(unique(y[is.finite(y)]))
    K <- length(u)
    contiguous <- K >= 3 && K <= 12 && all(diff(u) == 1) && min(u) %in% c(0, 1)
    if (contiguous) {
      trials <- as.integer(max(u))
      y_cnt <- as.integer(pmin(pmax(round(y), 0), trials))
      return(rtn(y_cnt, binomial("logit"), sprintf("(binomial trials=%d, inferred)", trials),
        identity, identity,
        extra = list(two_col = TRUE, trials = trials)
      ))
    }

    # Large integer counts are treated as continuous (log-normal approximation)
    nuniq <- length(unique(y[is.finite(y)]))
    if (nuniq > 30 && max(y, na.rm = TRUE) > 5) {
      return(rtn(
        log1p(pmax(y, 0)), gaussian(), "(log1p/gauss, int-but-cont)",
        function(x) log1p(pmax(x, 0)), function(mu) pmax(expm1(mu), 0)
      ))
    }

    # Use Negative Binomial if overdispersed, otherwise Poisson
    m <- mean(y, na.rm = TRUE)
    v <- var(y, na.rm = TRUE)
    use_nb <- is.finite(m) && m > 0 && is.finite(v) && v > 1.5 * m
    if (use_nb) {
      fam <- try(mgcv::nb(), silent = TRUE)
      if (!inherits(fam, "try-error")) {
        return(rtn(y, fam, paste0("(", fam$family, ")"), extra = list(is_nb = TRUE)))
      }
    }
    return(rtn(y, poisson("log"), "(poisson)", extra = list(is_nb = FALSE)))
  }

  # Default continuous: log1p if skewed positive, else z-score
  if (all(is.na(y) | (y >= 0))) {
    spread <- if (is.finite(q[2]) && q[2] > 0) (q[3] - q[1]) / q[2] else Inf
    if (is.finite(spread) && spread > 4) {
      return(rtn(log1p(y), gaussian(), "(log1p/gauss)", log1p, function(mu) pmax(expm1(mu), 0)))
    }
  }

  mu0 <- mean(y, na.rm = TRUE)
  sd0 <- stats::sd(y, na.rm = TRUE)
  if (!is.finite(sd0) || sd0 == 0) sd0 <- 1
  rtn((y - mu0) / sd0, gaussian(), "(z/gauss)", function(x) (x - mu0) / sd0, function(mu) mu * sd0 + mu0)
}

# ==============================================================================
# 5. Modelling & Metric Helpers
# ==============================================================================

is_binary_vec <- function(y) {
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  length(y) > 1L && all(y %in% c(0, 1))
}

.y_cmp_from_d <- function(yi, d, is_two_col, is_ord) {
  if (is_two_col) {
    den <- pmax(d$ys + d$yf, 1L)
    pmin(pmax(d$ys / den, 0), 1)
  } else if (is_ord) {
    yi$inv(d$y)
  } else {
    if (is.function(yi$inv)) yi$inv(d$y) else as.numeric(d$y)
  }
}

.permute_d <- function(d, yi, is_two_col) {
  if (is_two_col) {
    n <- nrow(d)
    p <- ifelse(d$ys + d$yf > 0, d$ys / (d$ys + d$yf), 0.5)
    p <- sample(p, n, replace = FALSE)
    tot <- d$ys + d$yf
    ys <- rbinom(n, size = tot, prob = p)
    d$ys <- ys
    d$yf <- tot - ys
    d
  } else {
    d$y <- sample(d$y, length(d$y), replace = FALSE)
    d
  }
}

make_folds_strat_general <- function(y, K, seed = 1L) {
  set.seed(seed)
  n <- length(y)
  fid <- integer(n)
  if (is_binary_vec(y)) {
    i1 <- which(y == 1)
    i0 <- which(y == 0)
    fid[i1] <- sample(rep(seq_len(K), length.out = length(i1)))
    fid[i0] <- sample(rep(seq_len(K), length.out = length(i0)))
  } else {
    rk <- rank(y, ties.method = "first", na.last = "keep")
    q <- cut(rk, breaks = quantile(rk, probs = seq(0, 1, length.out = 6), na.rm = TRUE), include.lowest = TRUE)
    for (g in levels(q)) {
      idx <- which(q == g)
      if (!length(idx)) next
      fid[idx] <- sample(rep(seq_len(K), length.out = length(idx)))
    }
  }
  if (any(fid == 0)) fid[fid == 0] <- sample(seq_len(K), sum(fid == 0), TRUE)
  fid
}

auc_point <- function(y, p) {
  y <- as.integer(y > 0)
  p <- as.numeric(p)
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]
  if (length(y) < 2L || !any(y == 1L) || !any(y == 0L)) {
    return(NA_real_)
  }

  rk <- rank(p, ties.method = "average")
  P <- sum(y == 1)
  N <- sum(y == 0)
  (sum(rk[y == 1]) - P * (P + 1) / 2) / (P * N)
}

auprc_point <- function(y, p) {
  y <- as.integer(y > 0)
  p <- as.numeric(p)
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]
  n <- length(y)
  P <- sum(y == 1L)
  N <- n - P
  if (n < 2L || P == 0L || N == 0L) {
    return(NA_real_)
  }

  o <- order(p, decreasing = TRUE)
  y <- y[o]
  tp <- cumsum(y == 1L)
  fp <- cumsum(y == 0L)
  prec <- tp / pmax(tp + fp, 1L)
  sum(prec[y == 1L]) / P
}

# Precision-Recall-Gain (Flach & Cullignana, 2015)
prg_flach <- function(y, p, eps = 1e-12) {
  y <- as.integer(y > 0)
  p <- as.numeric(p)
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]
  if (length(y) == 0L || length(unique(y)) < 2L) {
    return(list(auprg = NA_real_))
  }

  pi <- mean(y == 1L)
  o <- order(p, decreasing = TRUE)
  y <- y[o]
  tp <- cumsum(y == 1L)
  fp <- cumsum(y == 0L)
  n1 <- sum(y == 1L)
  R <- tp / pmax(n1, 1L)
  Pp <- tp / pmax(tp + fp, 1L)
  R <- c(0, R, 1)
  Pp <- c(pi, Pp, pi)
  PR <- data.frame(R = R, P = Pp)
  PR <- PR[order(PR$R, -PR$P), ]
  PR <- aggregate(P ~ R, PR, max)
  PR <- PR[order(PR$R), ]
  recG <- (PR$R - pi) / ((1 - pi) * pmax(PR$R, eps))
  precG <- (PR$P - pi) / ((1 - pi) * pmax(PR$P, eps))
  recG <- pmin(pmax(recG, 0), 1)
  precG <- pmin(pmax(precG, 0), 1)
  o2 <- order(recG)
  recG <- recG[o2]
  precG <- precG[o2]
  auprg <- if (length(recG) >= 2L) sum(0.5 * (precG[-1] + precG[-length(precG)]) * diff(recG)) else NA_real_
  list(auprg = as.numeric(auprg))
}

cindex_point <- function(y, s, max_pairs = 2e6, seed = 1L) {
  y <- as.numeric(y)
  s <- as.numeric(s)
  ok <- is.finite(y) & is.finite(s)
  y <- y[ok]
  s <- s[ok]
  if (length(unique(y)) < 2L) {
    return(NA_real_)
  }

  n <- length(y)
  set.seed(seed)
  M <- min(max_pairs, as.numeric(n) * max(1, floor(log1p(n))))
  i <- sample.int(n, M, replace = TRUE)
  j <- sample.int(n, M, replace = TRUE)
  keep <- i != j
  i <- i[keep]
  j <- j[keep]
  dy <- y[i] - y[j]
  ds <- s[i] - s[j]
  comp <- dy != 0
  dy <- dy[comp]
  ds <- ds[comp]
  if (!length(dy)) {
    return(NA_real_)
  }

  concord <- sum(dy * ds > 0)
  discord <- sum(dy * ds < 0)
  tie_s <- sum(ds == 0)
  (concord + 0.5 * tie_s) / (concord + discord + tie_s)
}

pr80_20_from_score <- function(y01, score, frac = 0.20, eps = 1e-6) {
  y01 <- as.integer(y01 > 0)
  s <- as.numeric(score)
  ok <- is.finite(y01) & is.finite(s)
  y01 <- y01[ok]
  s <- s[ok]
  n <- length(y01)
  if (n < 10L || length(unique(y01)) < 2L) {
    return(list(pr80_20 = NA_real_, delta_pp = NA_real_))
  }

  rk <- rank(-s, ties.method = "first")
  m <- max(1L, floor(frac * n))
  prev_top <- mean(y01[rk <= m] == 1L)
  prev_bot <- mean(y01[rk > (n - m)] == 1L)

  list(
    pr80_20  = (prev_top + eps) / (prev_bot + eps), # Ratio
    delta_pp = prev_top - prev_bot # Absolute Risk Increase
  )
}

# Calculates clinical utility metrics (Lift, NNS, NPV, Sens, Spec) for high-risk fractions.
clinical_metrics_from_score <- function(y01, score, frac = 0.20) {
  y01 <- as.integer(y01 > 0)
  s   <- as.numeric(score)
  ok  <- is.finite(y01) & is.finite(s)
  y01 <- y01[ok]
  s   <- s[ok]
  n   <- length(y01)
  
  if (n < 20L || sum(y01) == 0 || sum(y01 == 0) == 0) {
    return(list(
      lift = NA, nns_model = NA, nns_base = NA, 
      sens_20 = NA, spec_20 = NA, ppv_20 = NA, npv_20 = NA,
      sens_opt = NA, spec_opt = NA, ppv_opt = NA, npv_opt = NA, thresh_opt = NA
    ))
  }
  
  prev_global <- mean(y01 == 1L)
  nns_base    <- if (prev_global > 0) 1 / prev_global else NA_real_
  
  # --- A. The "Top 20%" Strategy (Screening/Enrichment) ---
  rk_top <- rank(-s, ties.method = "first")
  n_sel  <- max(1L, floor(frac * n))
  is_hr  <- rk_top <= n_sel
  
  TP_20 <- sum(is_hr & y01 == 1L)
  FP_20 <- sum(is_hr & y01 == 0L)
  FN_20 <- sum(!is_hr & y01 == 1L)
  TN_20 <- sum(!is_hr & y01 == 0L)
  
  sens_20 <- TP_20 / (TP_20 + FN_20)
  spec_20 <- TN_20 / (TN_20 + FP_20)
  ppv_20  <- if ((TP_20 + FP_20) > 0) TP_20 / (TP_20 + FP_20) else 0
  npv_20  <- if ((TN_20 + FN_20) > 0) TN_20 / (TN_20 + FN_20) else 0
  
  # --- B. The "Optimal" Strategy (Youden's J) ---
  # Sort by score to compute cumulative TP/FP for every possible threshold
  ord <- order(s, decreasing = TRUE)
  y_sorted <- y01[ord]
  s_sorted <- s[ord]
  
  # Cumulative counts (assuming threshold is just below this point)
  TP_vec <- cumsum(y_sorted == 1L)
  FP_vec <- cumsum(y_sorted == 0L)
  P_tot  <- sum(y01 == 1L)
  N_tot  <- sum(y01 == 0L)
  FN_vec <- P_tot - TP_vec
  TN_vec <- N_tot - FP_vec
  
  # Calculate Sens/Spec for ALL thresholds simultaneously
  sens_vec <- TP_vec / P_tot
  spec_vec <- TN_vec / N_tot
  
  # Youden's J = Sens + Spec - 1
  j_vec <- sens_vec + spec_vec - 1
  best_idx <- which.max(j_vec)
  
  # Extract metrics at the optimal cutpoint
  sens_opt <- sens_vec[best_idx]
  spec_opt <- spec_vec[best_idx]
  # For PPV/NPV at optimal cut
  tp_opt <- TP_vec[best_idx]
  fp_opt <- FP_vec[best_idx]
  fn_opt <- FN_vec[best_idx]
  tn_opt <- TN_vec[best_idx]
  
  ppv_opt <- if ((tp_opt + fp_opt) > 0) tp_opt / (tp_opt + fp_opt) else 0
  npv_opt <- if ((tn_opt + fn_opt) > 0) tn_opt / (tn_opt + fn_opt) else 0
  thresh_opt <- s_sorted[best_idx]
  
  list(
    # Enrichment stats (keep these for Lift/NNS)
    lift      = if (prev_global > 0) ppv_20 / prev_global else NA_real_,
    nns_model = if (ppv_20 > 0) 1 / ppv_20 else NA_real_,
    nns_base  = nns_base,
    
    # 20% Cut metrics
    sens_20 = sens_20, spec_20 = spec_20, ppv_20 = ppv_20, npv_20 = npv_20,
    
    # Optimal Cut metrics
    sens_opt = sens_opt, spec_opt = spec_opt, 
    ppv_opt = ppv_opt, npv_opt = npv_opt, thresh_opt = thresh_opt
  )
}
# --- Bootstrapping Helpers ---

boot_ci_balanced <- function(y, p, FUN, B = 1000L, seed = 42L) {
  set.seed(seed)
  y <- as.integer(y > 0)
  p <- as.numeric(p)
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]
  i0 <- which(y == 0L)
  i1 <- which(y == 1L)
  pt <- suppressWarnings(as.numeric(FUN(y, p)))

  if (!length(i0) || !length(i1) || !is.finite(pt)) {
    return(c(point = pt, lo = NA, hi = NA))
  }

  vals <- replicate(B, {
    ii <- c(sample(i0, length(i0), TRUE), sample(i1, length(i1), TRUE))
    suppressWarnings(as.numeric(FUN(y[ii], p[ii])))
  })
  vals <- vals[is.finite(vals)]
  if (!length(vals)) {
    return(c(point = pt, lo = NA, hi = NA))
  }

  c(point = pt, lo = as.numeric(quantile(vals, 0.025)), hi = as.numeric(quantile(vals, 0.975)))
}

boot_cindex_single <- function(y, s, B = 999L, seed = 1L, max_pairs = 2e6) {
  set.seed(seed)
  y <- as.numeric(y)
  s <- as.numeric(s)
  ok <- is.finite(y) & is.finite(s)
  y <- y[ok]
  s <- s[ok]

  if (length(unique(y)) < 2L) {
    return(c(point = NA, lo = NA, hi = NA, p = NA, win = NA))
  }

  pt <- cindex_point(y, s, max_pairs = max_pairs, seed = seed)
  draws <- replicate(B, {
    ii <- sample.int(length(y), length(y), TRUE)
    cindex_point(y[ii], s[ii], max_pairs = max_pairs, seed = sample.int(1e9, 1))
  })

  vf <- draws[is.finite(draws)]
  if (!length(vf)) {
    return(c(point = pt, lo = NA, hi = NA, p = NA, win = NA))
  }

  lo <- quantile(vf, 0.025)
  hi <- quantile(vf, 0.975)
  win <- mean(vf > 0.5)
  p2 <- min(1, 2 * min(mean(vf >= 0.5), mean(vf <= 0.5)))
  c(point = pt, lo = as.numeric(lo), hi = as.numeric(hi), p = p2, win = win)
}

oof_triplet_ci <- function(y, p, B = 1000L, seed = 42L) {
  list(
    auc = boot_ci_balanced(y, p, auc_point, B = B, seed = seed),
    prc = boot_ci_balanced(y, p, auprc_point, B = B, seed = seed),
    prg = boot_ci_balanced(y, p, function(yy, pp) prg_flach(yy, pp)$auprg, B = B, seed = seed)
  )
}

# --- Out-of-Fold (OOF) & Permutation Logic ---

oof_scores_single <- function(yi, d, form, K, seed = 1L, is_two_col = FALSE, is_ord = FALSE, levels = NULL) {
  set.seed(seed)
  n <- nrow(d)
  if (!n) {
    return(numeric(0))
  }

  y_strat <- if (is_two_col) {
    den <- pmax(d$ys + d$yf, 1L)
    pmin(pmax(d$ys / den, 0), 1)
  } else {
    d$y
  }

  fid <- make_folds_strat_general(y_strat, K = max(2L, min(K, 10L)), seed = seed)
  p_hat <- rep(NA_real_, n)

  for (k in seq_len(max(fid))) {
    tr <- fid != k
    te <- fid == k
    if (!any(tr) || !any(te)) next
    fit <- fit_gam_safe(form, yi$family, d[tr, , drop = FALSE], method = "REML")
    if (!is.null(fit)) {
      pr <- predict(fit, newdata = d[te, , drop = FALSE], type = "response")
      if (is.matrix(pr) && isTRUE(is_ord) && !is.null(levels)) pr <- as.numeric(pr %*% levels)
      p_hat[te] <- as.numeric(pr)
    } else {
      p_hat[te] <- mean(y_strat[tr], na.rm = TRUE)
    }
  }
  if (inherits(yi$family, "family") && grepl("binomial", yi$family$family, TRUE)) p_hat <- pmin(pmax(p_hat, 1e-6), 1 - 1e-6)
  p_hat
}

oof_metric_single <- function(yi, d, form, K, seed, is_two_col, is_ord, levels, max_pairs) {
  p_hat <- oof_scores_single(yi, d, form, K, seed, is_two_col, is_ord, levels)
  y_cmp <- .y_cmp_from_d(yi, d, is_two_col, is_ord)
  if (length(unique(y_cmp)) == 2 && all(y_cmp %in% c(0, 1))) {
    auc_point(as.integer(y_cmp > 0), p_hat)
  } else {
    cindex_point(y_cmp, p_hat, max_pairs = max_pairs, seed = seed + 7L)
  }
}

perm_test_oof <- function(yi, d, form, K, metric_is_bin, B, seed, is_two_col, is_ord, levels, max_pairs) {
  set.seed(seed)
  stat_obs <- oof_metric_single(yi, d, form, K, seed, is_two_col, is_ord, levels, max_pairs)
  stats <- rep(NA_real_, B)
  for (b in seq_len(B)) {
    dd <- .permute_d(d, yi, is_two_col)
    stats[b] <- oof_metric_single(yi, dd, form, K, seed + b, is_two_col, is_ord, levels, max_pairs)
  }
  stats <- stats[is.finite(stats)]
  if (!length(stats) || !is.finite(stat_obs)) {
    return(list(stat_obs = stat_obs, p = NA, null_med = NA, null_q95 = NA))
  }
  p <- (1 + sum(stats >= stat_obs)) / (length(stats) + 1)
  list(stat_obs = stat_obs, p = p, null_med = median(stats, na.rm = TRUE), null_q95 = as.numeric(quantile(stats, 0.95, na.rm = TRUE)))
}

# --- Stacking Helpers ---

fit_glm_or_glmnet <- function(y, X) {
  df <- as.data.frame(X)
  yb <- as.integer(y > 0)
  df$.y <- yb
  n1 <- sum(yb == 1L)
  n0 <- sum(yb == 0L)

  if (requireNamespace("glmnet", quietly = TRUE) && n1 >= 8L && n0 >= 8L) {
    x_mat <- as.matrix(df[setdiff(names(df), ".y")])
    cv <- try(glmnet::cv.glmnet(x_mat, yb, alpha = 0, family = "binomial", parallel = FALSE), silent = TRUE)
    if (!inherits(cv, "try-error")) {
      return(list(type = "glmnet", fit = cv, xnames = colnames(x_mat)))
    }
  }
  list(type = "glm", fit = stats::glm(.y ~ ., data = df, family = stats::binomial()))
}

pred_prob <- function(mod, newX) {
  if (mod$type == "glmnet") {
    x <- as.matrix(as.data.frame(newX)[, mod$xnames, drop = FALSE])
    p <- stats::predict(mod$fit, x, s = "lambda.min", type = "response")
  } else {
    p <- stats::predict(mod$fit, newdata = as.data.frame(newX), type = "response")
  }
  as.numeric(pmin(pmax(p, 1e-6), 1 - 1e-6))
}

oof_prob_stacked_custom <- function(y, Base_df, XR, K, seed) {
  y <- as.integer(y > 0)
  fid <- make_folds_strat_general(y, K, seed)
  pB <- rep(NA, length(y))
  pR <- rep(NA, length(y))
  pBR <- rep(NA, length(y))

  for (k in seq_len(K)) {
    tr <- fid != k
    te <- fid == k
    if (!any(tr) || !any(te)) next
    modB <- fit_glm_or_glmnet(y[tr], Base_df[tr, , drop = FALSE])
    pB[te] <- pred_prob(modB, Base_df[te, , drop = FALSE])
    modR <- fit_glm_or_glmnet(y[tr], XR[tr, , drop = FALSE])
    pR[te] <- pred_prob(modR, XR[te, , drop = FALSE])

    Xm_tr <- data.frame(l1 = qlogis(pmin(pmax(pB[tr], 1e-6), 1 - 1e-6)), l2 = qlogis(pmin(pmax(pR[tr], 1e-6), 1 - 1e-6)))
    m <- try(stats::glm(y[tr] ~ ., data = Xm_tr, family = stats::binomial()), silent = TRUE)
    if (!inherits(m, "try-error")) {
      Xm_te <- data.frame(l1 = qlogis(pmin(pmax(pB[te], 1e-6), 1 - 1e-6)), l2 = qlogis(pmin(pmax(pR[te], 1e-6), 1 - 1e-6)))
      pBR[te] <- as.numeric(stats::predict(m, newdata = Xm_te, type = "response"))
    } else {
      pBR[te] <- (pB[te] + pR[te]) / 2
    }
  }
  list(Base = pB, Resid = pR, Both = pBR)
}

fit_gam_safe <- function(form, family, data, method = "REML", gamma_val = 1.8) {
  fx <- try(mgcv::gam(
    form, family = family, data = data,
    method = method, select = TRUE, gamma = gamma_val,
    control = mgcv::gam.control(maxit = 400)
  ), silent = TRUE)
  
  # Check for convergence failure or infinite coefficients
  if (inherits(fx, "try-error") || any(!is.finite(tryCatch(stats::coef(fx), error=function(e) NA)))) {
    return(NULL)
  }
  fx
}

trim_gam_safe <- function(fit, keep_model = TRUE) {
  fn <- try(getFromNamespace("trim.gam", "mgcv"), silent = TRUE)
  if (is.function(fn)) {
    out <- try(fn(fit), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(out)
    }
  }
  drop <- c("y", "fitted.values", "linear.predictors", "working.weights", "weights", "offset", "residuals", "prior.weights", "hat", "qr", "X", "R")
  if (!keep_model) drop <- c(drop, "model")
  for (nm in intersect(names(fit), drop)) fit[[nm]] <- NULL
  fit
}

dev_expl_from_fit <- function(fit) {
  if (is.null(fit)) {
    return(NA_real_)
  }
  s <- try(summary(fit), silent = TRUE)
  if (inherits(s, "try-error")) {
    return(NA_real_)
  }
  as.numeric(s$dev.expl)
}
edf_total <- function(fit) {
  if (is.null(fit)) {
    return(NA_real_)
  }
  s <- try(summary(fit), silent = TRUE)
  if (inherits(s, "try-error") || is.null(s$s.table)) {
    return(NA_real_)
  }
  sum(suppressWarnings(as.numeric(s$s.table[, "edf"])), na.rm = TRUE)
}
edf_ti <- function(fit) {
  st <- try(suppressWarnings(summary(fit)$s.table), silent = TRUE)
  if (inherits(st, "try-error") || is.null(st)) {
    return(NA_real_)
  }
  ii <- grep("^ti\\(", rownames(st))
  if (!length(ii)) {
    return(0)
  }
  sum(suppressWarnings(as.numeric(st[ii, "edf"])), na.rm = TRUE)
}

akaike_weight <- function(dAIC) if (is.finite(dAIC)) 1 / (1 + exp(-0.5 * dAIC)) else NA_real_
akaike_odds <- function(dAIC) if (is.finite(dAIC)) exp(0.5 * dAIC) else NA_real_

lrt_p_ml_nosel <- function(form_add, form_full, family, data) {
  f1 <- try(mgcv::gam(form_add, family = family, data = data, method = "ML", select = FALSE), silent = TRUE)
  f2 <- try(mgcv::gam(form_full, family = family, data = data, method = "ML", select = FALSE), silent = TRUE)
  if (inherits(f1, "try-error") || inherits(f2, "try-error")) {
    return(NA_real_)
  }

  an <- try(mgcv::anova.gam(f1, f2, test = "Chisq"), silent = TRUE)
  if (inherits(an, "try-error") || nrow(an) < 2) {
    return(NA_real_)
  }
  as.numeric(an[2, ncol(an)])
}

make_folds <- function(y, K, seed = 1) {
  set.seed(seed)
  n <- length(y)
  fid <- integer(n)
  if (length(unique(y)) == 2 && all(y %in% c(0, 1))) {
    i1 <- which(y == 1)
    i0 <- which(y == 0)
    fid[i1] <- sample(rep(seq_len(K), length.out = length(i1)))
    fid[i0] <- sample(rep(seq_len(K), length.out = length(i0)))
  } else {
    rk <- rank(y, ties.method = "first")
    q <- cut(rk, breaks = quantile(rk, seq(0, 1, l = 6)), include.lowest = TRUE)
    for (g in levels(q)) {
      idx <- which(q == g)
      if (length(idx)) fid[idx] <- sample(rep(seq_len(K), length.out = length(idx)))
    }
  }
  if (any(fid == 0)) fid[fid == 0] <- sample(seq_len(K), sum(fid == 0), TRUE)
  fid
}

fit_glm_net <- function(y, X) {
  df <- as.data.frame(X)
  yb <- as.integer(y > 0)
  df$.y <- yb
  n1 <- sum(yb)
  n0 <- sum(yb == 0)

  if (requireNamespace("glmnet", quietly = TRUE) && n1 >= 8 && n0 >= 8) {
    xm <- as.matrix(df[setdiff(names(df), ".y")])
    fid <- make_folds(yb, 5, 42)
    cv <- try(glmnet::cv.glmnet(xm, yb, alpha = 0, family = "binomial", foldid = fid, parallel = FALSE), silent = TRUE)
    if (!inherits(cv, "try-error")) {
      return(list(type = "glmnet", fit = cv, xnames = colnames(xm)))
    }
  }
  list(type = "glm", fit = glm(.y ~ ., data = df, family = binomial))
}

oof_stacked <- function(y, Base, XR, K) {
  y <- as.integer(y > 0)
  fid <- make_folds(y, K, 42)
  pB <- rep(NA, length(y))
  pR <- rep(NA, length(y))
  pBR <- rep(NA, length(y))

  for (k in 1:K) {
    tr <- fid != k
    te <- fid == k
    if (!any(tr) || !any(te)) next
    mB <- fit_glm_net(y[tr], Base[tr, ])
    pB[te] <- pred_prob(mB, Base[te, ])
    mR <- fit_glm_net(y[tr], XR[tr, ])
    pR[te] <- pred_prob(mR, XR[te, ])

    Xst <- data.frame(l1 = qlogis(pB[tr]), l2 = qlogis(pR[tr]), y = y[tr])
    mS <- try(glm(y ~ ., data = Xst, family = binomial), silent = TRUE)
    if (!inherits(mS, "try-error")) {
      pBR[te] <- predict(mS, newdata = data.frame(l1 = qlogis(pB[te]), l2 = qlogis(pR[te])), type = "response")
    } else {
      pBR[te] <- (pB[te] + pR[te]) / 2
    }
  }
  list(Base = pB, Resid = pR, Both = pBR)
}

boot_ci <- function(y, p, FUN, B = 200) {
  pt <- FUN(y, p)
  if (!is.finite(pt)) {
    return(c(NA, NA, NA))
  }
  i0 <- which(y == 0)
  i1 <- which(y == 1)
  v <- replicate(B, {
    ii <- c(sample(i0, length(i0), T), sample(i1, length(i1), T))
    FUN(y[ii], p[ii])
  })
  c(pt, quantile(v, c(0.025, 0.975), na.rm = TRUE))
}

# Calculates calibration intercept, slope, Brier score, and ECE.
calibration_metrics_binary <- function(y, p, min_bin = 25L, max_bins = 8L) {
  y <- as.integer(y > 0)
  p <- as.numeric(p)
  ok <- is.finite(y) & is.finite(p)
  y <- y[ok]
  p <- p[ok]
  if (length(unique(y)) < 2L || length(y) < 30L) {
    return(list(
      intercept = NA_real_, slope = NA_real_, brier = NA_real_, ece = NA_real_, n_bins = 0L, n_total = length(y),
      curve = data.frame(p_hat = numeric(0), obs = numeric(0), n = integer(0)),
      line = data.frame(p_hat = numeric(0), p = numeric(0), lo = numeric(0), hi = numeric(0))
    ))
  }

  p_clip <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  lp <- qlogis(p_clip)
  fit_cal <- try(glm(y ~ lp, family = binomial), silent = TRUE)
  intercept <- slope <- NA_real_
  line <- NULL

  if (!inherits(fit_cal, "try-error")) {
    cf <- coef(fit_cal)
    if (length(cf) >= 2 && all(is.finite(cf[1:2]))) {
      intercept <- unname(cf[1])
      slope <- unname(cf[2])
    }
    gp <- seq(0, 1, l = 200)
    pr <- predict(fit_cal, newdata = data.frame(lp = qlogis(pmin(pmax(gp, 1e-6), 1 - 1e-6))), type = "link", se.fit = TRUE)
    line <- data.frame(p_hat = gp, p = plogis(pr$fit), lo = plogis(pr$fit - 1.96 * pr$se.fit), hi = plogis(pr$fit + 1.96 * pr$se.fit))
  }

  brier <- mean((p_clip - y)^2)
  br <- unique(quantile(p_clip, probs = seq(0, 1, length.out = max(3, min(max_bins, floor(length(p_clip) / min_bin))) + 1), na.rm = TRUE, names = FALSE))
  if (length(br) < 3) br <- c(0, 0.5, 1)
  g <- cut(p_clip, breaks = br, include.lowest = TRUE)
  T <- data.frame(p_hat = tapply(p_clip, g, mean), obs = tapply(y, g, mean), n = tapply(y, g, length))

  list(
    intercept = intercept, slope = slope, brier = brier, ece = weighted.mean(abs(T$p_hat - T$obs), w = T$n),
    n_bins = nrow(T), n_total = length(y), curve = T, line = line
  )
}

# --- helper: build a fresh family object for refits (critical for NB + OCAT) ---
fresh_family <- function(family, data = NULL) {
  if (is.function(family)) return(family())
  
  if (inherits(family, "family")) {
    fam <- family$family %||% ""
    
    # OCAT: rebuild from K (number of categories)
    if (grepl("^Ordered Categorical", fam)) {
      # K derived from data if available, else from fitted y later
      K <- NA_integer_
      if (!is.null(data)) {
        yv <- data$y
        yv <- yv[is.finite(yv)]
        if (length(yv)) K <- as.integer(max(yv))
      }
      if (!is.finite(K) || K < 2L) {
        # fall back to 2 to avoid crash; caller should ensure K is right
        K <- 2L
      }
      return(ocat_family(K))
    }
    
    # NB: rebuild so theta is re-estimated on each refit
    if (grepl("^Negative Binomial", fam)) {
      return(mgcv::nb())
    }
  }
  
  family
}

# --- helper: extract NB theta robustly ---
get_nb_theta <- function(fit) {
  th <- NA_real_
  if (!is.null(fit$family)) {
    if (!is.null(fit$family$getTheta)) {
      th <- suppressWarnings(try(as.numeric(fit$family$getTheta(TRUE)), silent = TRUE))
      if (inherits(th, "try-error")) th <- NA_real_
    }
    if (!is.finite(th) && !is.null(fit$family$theta)) th <- as.numeric(fit$family$theta)
  }
  if (!is.finite(th) && !is.null(fit$theta)) th <- as.numeric(fit$theta)
  th
}

# Parametric simulation for OCAT and NB
simulate_from_null <- function(fit_null, data, B,
                               is_two_col = FALSE, trials = NULL,
                               is_ordinal = FALSE, levels = NULL) {
  fam_name <- ""
  if (!is.null(fit_null$family) && inherits(fit_null$family, "family")) {
    fam_name <- fit_null$family$family %||% ""
  }
  
  n <- nrow(data)
  
  # Binomial trials (two-column): simulate successes directly
  if (is_two_col) {
    if (is.null(trials)) stop("simulate_from_null: trials must be provided when is_two_col=TRUE")
    p <- as.numeric(fitted(fit_null))
    p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
    size <- if (length(trials) == 1L) rep(trials, n) else as.integer(trials)
    return(replicate(B, rbinom(n, size = size, prob = p), simplify = FALSE))
  }
  
  # OCAT: simulate categories from predicted category probabilities
  if (grepl("^Ordered Categorical", fam_name) || isTRUE(is_ordinal)) {
    P <- suppressWarnings(try(predict(fit_null, newdata = data, type = "response"), silent = TRUE))
    if (inherits(P, "try-error") || is.null(P)) return(NULL)
    if (!is.matrix(P)) {
      # If mgcv ever returns something unexpected, abort cleanly
      return(NULL)
    }
    K <- ncol(P)
    # normalise rows
    rs <- rowSums(P)
    rs[!is.finite(rs) | rs <= 0] <- 1
    P <- P / rs
    
    return(replicate(B, {
      ysim <- integer(n)
      for (i in seq_len(n)) {
        pr <- P[i, ]
        pr[!is.finite(pr) | pr < 0] <- 0
        s <- sum(pr)
        if (s <= 0) pr[] <- 1 / K else pr <- pr / s
        ysim[i] <- sample.int(K, size = 1, prob = pr)
      }
      ysim
    }, simplify = FALSE))
  }
  
  # NB: simulate from NB(mu, theta)
  if (grepl("^Negative Binomial", fam_name)) {
    mu <- as.numeric(fitted(fit_null))
    mu[!is.finite(mu) | mu < 0] <- 0
    th <- get_nb_theta(fit_null)
    if (!is.finite(th) || th <= 0) return(NULL)
    
    return(replicate(B, rnbinom(n, mu = mu, size = th), simplify = FALSE))
  }
  
  # Default: try mgcv simulate (gaussian/poisson/binomial etc.)
  S <- suppressWarnings(try(simulate(fit_null, nsim = B), silent = TRUE))
  if (inherits(S, "try-error") || is.null(S)) return(NULL)
  S
}

lrt_boot_gam <- function(fit_null, form_null, form_alt, family, data,
                         B = 199L, method = "ML", select = TRUE,
                         is_two_col = FALSE, trials = NULL,
                         is_ordinal = FALSE, levels = NULL) {
  
  fam2 <- fresh_family(family, data)
  f2 <- try(mgcv::gam(form_alt, family = fam2, data = data, method = method, select = select), silent = TRUE)
  if (inherits(f2, "try-error")) return(NA_real_)
  
  stat_obs <- as.numeric(2 * (logLik(f2) - logLik(fit_null)))
  if (!is.finite(stat_obs)) return(NA_real_)
  
  S <- simulate_from_null(fit_null, data, B,
                          is_two_col = is_two_col, trials = trials,
                          is_ordinal = is_ordinal, levels = levels)
  if (is.null(S) || !length(S)) return(NA_real_)
  
  stat_rep <- rep(NA_real_, length(S))
  
  for (b in seq_along(S)) {
    d_sim <- data
    if (is_two_col) {
      ys <- as.numeric(S[[b]])
      d_sim$ys <- ys
      d_sim$yf <- as.numeric(trials) - ys
    } else {
      d_sim$y <- as.numeric(S[[b]])
    }
    
    fam_b1 <- fresh_family(family, d_sim)
    fam_b2 <- fresh_family(family, d_sim)
    
    f1b <- try(mgcv::gam(form_null, family = fam_b1, data = d_sim, method = method, select = select), silent = TRUE)
    f2b <- try(mgcv::gam(form_alt,  family = fam_b2, data = d_sim, method = method, select = select), silent = TRUE)
    
    if (!inherits(f1b, "try-error") && !inherits(f2b, "try-error")) {
      stat_rep[b] <- as.numeric(2 * (logLik(f2b) - logLik(f1b)))
    }
  }
  
  stat_rep <- stat_rep[is.finite(stat_rep)]
  if (!length(stat_rep)) return(NA_real_)
  
  mean(stat_rep >= stat_obs)
}

# ==============================================================================
# 6. PD / PDP (Principal Direction and Partial Dependence)
# ==============================================================================

principal_direction_from_fit <- function(fit2d, geom, Bobs = NULL) {
  G <- geom$gridB_full
  H <- geom$mask_hull %||% rep(TRUE, nrow(G))
  hb1 <- median(diff(sort(unique(G$b1))), na.rm = TRUE)
  hb2 <- median(diff(sort(unique(G$b2))), na.rm = TRUE)
  if (!is.finite(hb1) || hb1 <= 0) hb1 <- 1e-3
  if (!is.finite(hb2) || hb2 <= 0) hb2 <- 1e-3

  e_p1 <- as.numeric(predict(fit2d, newdata = transform(G, b1 = b1 + hb1), type = "link"))
  e_m1 <- as.numeric(predict(fit2d, newdata = transform(G, b1 = b1 - hb1), type = "link"))
  e_p2 <- as.numeric(predict(fit2d, newdata = transform(G, b2 = b2 + hb2), type = "link"))
  e_m2 <- as.numeric(predict(fit2d, newdata = transform(G, b2 = b2 - hb2), type = "link"))

  g1 <- (e_p1 - e_m1) / (2 * hb1)
  g2 <- (e_p2 - e_m2) / (2 * hb2)
  g1[!is.finite(g1)] <- 0
  g2[!is.finite(g2)] <- 0

  mag <- sqrt(g1^2 + g2^2)
  ok <- is.finite(mag) & H
  if (!any(ok)) {
    return(c(1, 0))
  }

  thr <- stats::quantile(mag[ok], 0.75, na.rm = TRUE)
  S <- which(ok & mag >= thr)
  if (!length(S)) S <- which(ok)

  u <- c(mean(g1[S], na.rm = TRUE), mean(g2[S], na.rm = TRUE))
  if (!all(is.finite(u)) || sum(u^2) <= 0) {
    return(c(1, 0))
  }
  u <- u / sqrt(sum(u^2))

  if (!is.null(Bobs) && all(c("b1", "b2") %in% names(Bobs))) {
    s_raw <- u[1] * Bobs$b1 + u[2] * Bobs$b2
    mu <- tryCatch(as.numeric(predict(fit2d, newdata = Bobs, type = "response")), error = function(e) NULL)
    if (!is.null(mu) && length(mu) == length(s_raw)) {
      rr <- suppressWarnings(cor(s_raw, mu, use = "complete", method = "spearman"))
      if (is.finite(rr) && rr < 0) u <- -u
    }
  }
  u
}

pd_project_z <- function(b1, b2, u, centre = c("median", "mean")) {
  centre <- match.arg(centre)
  s_raw <- u[1] * b1 + u[2] * b2
  c0 <- if (centre == "median") median(s_raw, na.rm = TRUE) else mean(s_raw, na.rm = TRUE)
  s_c <- s_raw - c0
  s_sd <- stats::sd(s_c, na.rm = TRUE)
  if (!is.finite(s_sd) || s_sd == 0) s_sd <- 1
  list(raw = s_raw, z = s_c / s_sd, centre = c0, scale = s_sd)
}

pdp_1d_from_2d <- function(fit2d, Uobs, geom, axis = c("u1", "u2"), n = 200,
                           weights = NULL,
                           is_ordinal = FALSE, levels = NULL, trials = NULL,
                           ci_draws = 300L, ci_obs = 1200L, ci_level = 0.95,
                           seed = 1L) {
  
  axis <- match.arg(axis)
  
  u1_all <- Uobs$u1
  u2_all <- Uobs$u2
  ok <- is.finite(u1_all) & is.finite(u2_all)
  u1_obs <- u1_all[ok]
  u2_obs <- u2_all[ok]
  stopifnot(length(u1_obs) == length(u2_obs))
  
  w <- if (is.null(weights)) rep(1, length(u1_obs)) else as.numeric(weights[ok])
  w[!is.finite(w) | w < 0] <- 0
  if (sum(w) <= 0) w[] <- 1
  w <- w / sum(w)
  
  grid <- seq(-1, 1, length.out = n)
  
  # covariance + coefficients
  Vp <- tryCatch(fit2d$Vp, error = function(e) NULL)
  if (is.null(Vp)) Vp <- tryCatch(vcov(fit2d), error = function(e) NULL)
  beta_hat <- tryCatch(coef(fit2d), error = function(e) NULL)
  
  out_x  <- rep(NA_real_, n)
  out_lo <- rep(NA_real_, n)
  out_hi <- rep(NA_real_, n)
  
  get_ocat_thresholds <- function(fit) {
    fam <- fit$family
    if (is.null(fam) || !inherits(fam, "family")) return(NULL)
    s <- fam$family %||% ""
    # extract numeric cutpoints from "Ordered Categorical(a,b,c,...)"
    nums <- regmatches(s, gregexpr("[-+]?[0-9]*\\.?[0-9]+(?:[eE][-+]?[0-9]+)?", s))[[1]]
    th <- suppressWarnings(as.numeric(nums))
    th <- th[is.finite(th)]
    if (!length(th)) NULL else th
  }
  
  # weighted mean of category probabilities under ordered-logit:
  # P(Y <= j) = plogis(theta_j - eta)
  ordinal_weighted_expected <- function(eta_vec, w_vec, theta, lev) {
    K <- length(theta) + 1L
    if (length(lev) != K) lev <- seq_len(K)
    
    # compute weighted CDFs for each threshold
    wcdf <- numeric(length(theta))
    for (j in seq_along(theta)) {
      wcdf[j] <- sum(w_vec * plogis(theta[j] - eta_vec))
    }
    
    pbar <- numeric(K)
    pbar[1] <- wcdf[1]
    if (K > 2L) {
      for (j in 2:(K - 1L)) pbar[j] <- wcdf[j] - wcdf[j - 1L]
    }
    pbar[K] <- 1 - wcdf[K - 1L]
    
    # clamp numerical noise
    pbar[pbar < 0] <- 0
    s <- sum(pbar)
    if (s > 0) pbar <- pbar / s else pbar <- rep(1 / K, K)
    
    sum(pbar * lev)
  }
  
  ordinal_weighted_expected_draws <- function(Xg, w_vec, theta, lev, beta_hat, Vp, S) {
    K <- length(theta) + 1L
    if (length(lev) != K) lev <- seq_len(K)
    
    # draw coefficients
    if (!requireNamespace("MASS", quietly = TRUE)) {
      stop("MASS package required for simulation-based CIs (MASS::mvrnorm).")
    }
    B <- MASS::mvrnorm(n = S, mu = beta_hat, Sigma = Vp)
    
    # eta draws: (n_obs x S)
    eta_draw <- Xg %*% t(B)
    
    # weighted CDFs per threshold, per draw: each is length S
    wcdf <- matrix(NA_real_, nrow = length(theta), ncol = S)
    for (j in seq_along(theta)) {
      # plogis on matrix is vectorised
      cdf_j <- plogis(theta[j] - eta_draw)
      wcdf[j, ] <- as.numeric(crossprod(w_vec, cdf_j))  # 1 x S
    }
    
    # category probs averaged over obs, per draw
    pbar <- matrix(0, nrow = K, ncol = S)
    pbar[1, ] <- wcdf[1, ]
    if (K > 2L) {
      for (j in 2:(K - 1L)) pbar[j, ] <- wcdf[j, ] - wcdf[j - 1L, ]
    }
    pbar[K, ] <- 1 - wcdf[K - 1L, ]
    
    # clamp numerical noise
    pbar[pbar < 0] <- 0
    col_sums <- colSums(pbar)
    col_sums[col_sums <= 0 | !is.finite(col_sums)] <- 1
    pbar <- sweep(pbar, 2, col_sums, "/")
    
    # expected value per draw
    as.numeric(crossprod(lev, pbar))  # 1 x S
  }
  
  # subsample obs for speed (important: done once)
  set.seed(seed)
  if (length(w) > ci_obs) {
    ii <- sample.int(length(w), size = ci_obs, replace = FALSE, prob = w)
    u1_obs_ci <- u1_obs[ii]
    u2_obs_ci <- u2_obs[ii]
    w_ci <- w[ii]
    w_ci <- w_ci / sum(w_ci)
  } else {
    u1_obs_ci <- u1_obs
    u2_obs_ci <- u2_obs
    w_ci <- w
  }
  
  # ordinal cutpoints if needed
  theta <- NULL
  if (isTRUE(is_ordinal)) {
    theta <- get_ocat_thresholds(fit2d)
    if (is.null(theta)) {
      # if thresholds cannot be recovered, we cannot do ordinal CI correctly
      warning("Ordinal requested but could not extract OCAT thresholds from fit$family$family.")
    }
  }
  
  alpha <- (1 - ci_level) / 2
  qlo <- alpha
  qhi <- 1 - alpha
  
  for (k in seq_along(grid)) {
    
    # build UB for this grid point
    if (axis == "u1") {
      UB <- as.data.frame(geom$Xstd$inv(rep(grid[k], length(u2_obs)), u2_obs))
      UB_ci <- as.data.frame(geom$Xstd$inv(rep(grid[k], length(u2_obs_ci)), u2_obs_ci))
    } else {
      UB <- as.data.frame(geom$Xstd$inv(u1_obs, rep(grid[k], length(u1_obs))))
      UB_ci <- as.data.frame(geom$Xstd$inv(u1_obs_ci, rep(grid[k], length(u1_obs_ci))))
    }
    names(UB) <- c("b1", "b2")
    names(UB_ci) <- c("b1", "b2")
    
    # Delta Method for ocat CIs
    # Calculates the variance of the expected value by propagating the 
    # coefficient covariance matrix (Vp) through the gradient of the 
    # expectation function.
    if (isTRUE(is_ordinal) && !is.null(theta)) {
      
      # 1. Linear Predictor Matrix and Eta
      # Xg is the design matrix for the current grid point marginalized over the secondary axis
      Xg <- predict(fit2d, newdata = UB, type = "lpmatrix")
      eta_hat <- as.numeric(Xg %*% beta_hat)
      
      # 2. Point Estimate: Weighted Expectation
      # E[y] = Sum(Score_k * P(y=k))
      out_x[k] <- ordinal_weighted_expected(eta_hat, w, theta, levels)
      
      # 3. Interval Estimate: Analytical Gradient Calculation
      if (!is.null(Vp) && !is.null(beta_hat)) {
        
        # Calculate the density of the logistic distribution at each cutpoint (theta)
        # for each observation. This represents the derivative of the cumulative 
        # probability function F(theta - eta) with respect to the linear predictor.
        # Dimension: [N_obs x N_cutpoints]
        f_mat <- outer(eta_hat, theta, function(e, t) stats::dlogis(t - e))
        
        # Weight the densities by the marginal integration weights (w)
        w_f_mat <- f_mat * w 
        
        # Compute gradients of the cumulative probabilities w.r.t. beta parameters.
        # By the chain rule: dF/d_beta = dF/d_eta * d_eta/d_beta
        # Since d_eta/d_beta = X, this projection effectively sums gradients 
        # across the marginal distribution.
        # Dimension: [N_coefficients x N_cutpoints]
        Grad_F <- t(Xg) %*% w_f_mat 
        
        # Map cumulative gradients to the gradient of the Expected Value.
        # Using the property: E[y] = Score_K - Sum_{j=1}^{K-1} (Score_{j+1} - Score_j) * F_j
        # We aggregate the gradients of F_j weighted by the step size between scores.
        Grad_E <- numeric(length(beta_hat))
        
        # Ensure levels vector matches theta length + 1 (K categories => K-1 cutpoints)
        K <- length(theta) + 1L
        lev_vec <- if (length(levels) == K) levels else seq_len(K)
        
        for (j in seq_along(theta)) {
          # The derivative contribution of cutpoint j to the Expectation
          score_diff <- lev_vec[j+1] - lev_vec[j]
          Grad_E <- Grad_E + (score_diff * Grad_F[, j])
        }
        
        # 4. Variance Propagation (First-order Taylor Expansion)
        # Var(E[y]) approx Grad_E^T * Vp * Grad_E
        var_E <- as.numeric(crossprod(Grad_E, Vp %*% Grad_E))
        se_E  <- sqrt(max(0, var_E))
        
        # 95% CI
        out_lo[k] <- out_x[k] - 1.96 * se_E
        out_hi[k] <- out_x[k] + 1.96 * se_E
      }
      
    } else { 
      
      # Gaussian/GLM Confidence Intervals
      # Standard error calculation for non-ordinal families using the linear 
      # predictor space and transforming via link inverse.
      fam <- fit2d$family
      linkinv <- fam$linkinv
      
      if (is.null(Vp) || is.null(beta_hat)) {
        pr <- as.numeric(predict(fit2d, newdata = UB, type = "response"))
        out_x[k] <- sum(w * pr)
      } else {
        Xg <- predict(fit2d, newdata = UB, type = "lpmatrix")
        cg <- as.numeric(crossprod(w, Xg))
        eta <- sum(cg * beta_hat)
        
        var_e <- as.numeric(cg %*% Vp %*% cg)
        se_e <- sqrt(max(var_e, 0))
        
        out_x[k]  <- linkinv(eta)
        out_lo[k] <- linkinv(eta - 1.96 * se_e)
        out_hi[k] <- linkinv(eta + 1.96 * se_e)
      }
    }
  }
  
  if (!is.null(trials)) {
    out_x <- out_x * trials
    out_lo <- out_lo * trials
    out_hi <- out_hi * trials
  }
  
  if (axis == "u1") {
    data.frame(u1 = grid, p = out_x, lo = out_lo, hi = out_hi)
  } else {
    data.frame(u2 = grid, p = out_x, lo = out_lo, hi = out_hi)
  }
}

pdp_pd_from_2d <- function(fit2d, Bobs, u, n = 200, weights = NULL) {
  b1 <- as.numeric(Bobs$b1)
  b2 <- as.numeric(Bobs$b2)
  ok <- is.finite(b1) & is.finite(b2)
  b1 <- b1[ok]
  b2 <- b2[ok]
  S <- pd_project_z(b1, b2, u)
  tz <- S$z
  t_grid <- seq(min(tz, na.rm = TRUE), max(tz, na.rm = TRUE), length.out = n)
  w <- if (is.null(weights)) rep(1, length(tz)) else as.numeric(weights[ok])
  w <- w / sum(w)
  fam <- fit2d$family
  linkinv <- fam$linkinv
  Vp <- tryCatch(fit2d$Vp, error = function(e) NULL)
  if (is.null(Vp)) Vp <- tryCatch(vcov(fit2d), error = function(e) NULL)
  beta <- tryCatch(coef(fit2d), error = function(e) NULL)
  r1 <- b1 - u[1] * S$raw
  r2 <- b2 - u[2] * S$raw

  out <- lapply(seq_along(t_grid), function(k) {
    sz <- t_grid[k] * S$scale + S$centre
    UB <- data.frame(b1 = r1 + u[1] * sz, b2 = r2 + u[2] * sz)
    if (is.null(Vp) || is.null(beta)) {
      pr <- as.numeric(predict(fit2d, newdata = UB, type = "response"))
      c(p = sum(w * pr), lo = NA_real_, hi = NA_real_)
    } else {
      Xg <- predict(fit2d, newdata = UB, type = "lpmatrix")
      cg <- as.numeric(crossprod(w, Xg))
      eta <- sum(cg * beta)
      var_e <- as.numeric(cg %*% Vp %*% cg)
      se_e <- sqrt(max(var_e, 0))
      c(p = linkinv(eta), lo = linkinv(eta - 1.96 * se_e), hi = linkinv(eta + 1.96 * se_e))
    }
  })
  D <- as.data.frame(do.call(rbind, out))
  cbind(t = t_grid, D)
}

# ==============================================================================
# 7. Surface Plotting & Subtitles
# ==============================================================================

subtitle_fit <- function(
  n, K, de_add, de_full, edf_used, kmin_used, dAIC, p_LRT, rho_s, used_model,
  r2_label = NULL, r2_value = NA_real_
) {
  star <- if (is.finite(p_LRT) && p_LRT < 0.05) " ★" else ""
  r2_txt <- if (!is.null(r2_label) && is.finite(r2_value)) sprintf(" | %s=%.3f", r2_label, r2_value) else ""
  paste0(
    "n=", n, ", K=", K, "  |  model=", used_model, "\n",
    "Dev.expl(add/full)=", sprintf("%.3f/%.3f", de_add, de_full),
    " | EDF_used=", sprintf("%.1f", edf_used),
    " | k-index(min)=", ifelse(is.finite(kmin_used), sprintf("%.3f", kmin_used), "NA"),
    " | ρ_s=", ifelse(is.finite(rho_s), sprintf("%.3f", rho_s), "NA"),
    r2_txt, "\n",
    "ΔAIC(add→full)=", ifelse(is.finite(dAIC), sprintf("%.1f", dAIC), "NA"),
    " | LRT p=", ifelse(is.finite(p_LRT), formatC(p_LRT, digits = 2, format = "e"), "NA"),
    star
  )
}

nice_decade_breaks <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) {
    return(NULL)
  }
  r <- range(x)
  if (!all(is.finite(r))) {
    return(NULL)
  }
  if (r[1] == r[2]) {
    return(r)
  }
  if (r[1] <= 0) {
    br <- pretty(r, n = 6)
    br <- br[br >= r[1] & br <= r[2]]
    return(if (length(br)) br else r)
  }
  e_min <- floor(log10(r[1]))
  e_max <- ceiling(log10(r[2]))
  bases <- c(1, 2, 5)
  br <- as.vector(outer(bases, 10^(e_min:e_max), `*`))
  br <- br[br >= r[1] & br <= r[2]]
  if (length(br) < 3) {
    br2 <- pretty(r, n = 5)
    br2 <- br2[br2 >= r[1] & br2 <= r[2]]
    if (length(br2) >= 3) br <- br2
  }
  if (!length(br)) NULL else br
}

plots_for_cont <- function(fields, geom, var_name, subtitle_txt, caption_txt, points, lims_override = NULL) {
  y_pts <- points$y
  yf <- y_pts[is.finite(y_pts)]
  is_prob_like <- length(yf) && min(yf) >= -1e-8 && max(yf) <= 1 + 1e-8
  lims_all <- if (is_prob_like) c(0, 1) else if (is.null(lims_override)) range(c(fields$Fbase$p, points$y), finite = TRUE) else lims_override
  base_scale <- scale_prob_fill(limits = lims_all, name = var_name)
  pt_scale <- scale_prob_colour(limits = lims_all, name = paste0(var_name, " (obs)"))
  vf <- fields$Fbase$p
  vf <- vf[is.finite(vf)]

  br <- try(
    {
      if (!length(vf)) {
        NULL
      } else if (is_prob_like) {
        cand <- c(0.10, 0.25, 0.50, 0.75, 0.90)
        levs <- cand[cand >= min(vf) & cand <= max(vf)]
        if (!length(levs)) NULL else levs
      } else {
        nice_decade_breaks(vf)
      }
    },
    silent = TRUE
  )
  if (inherits(br, "try-error")) br <- NULL

  addC <- function(df, xcol, ycol, breaks, is_prob_like) {
    if (is.null(breaks) || !length(breaks)) {
      return(NULL)
    }
    mapping <- if (is_prob_like) {
      ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]], z = p, label = after_stat(sprintf("%d%%", round(..level.. * 100))))
    } else {
      ggplot2::aes(x = .data[[xcol]], y = .data[[ycol]], z = p, label = after_stat(signif(..level.., 2)))
    }
    metR::geom_contour2(
      data = df, mapping = mapping, breaks = breaks, colour = "white", size = 0.28, label_colour = "white",
      skip = 0, label.placer = metR::label_placer_fraction(0.5), na.rm = TRUE, inherit.aes = FALSE
    )
  }

  theme_fix <- ggplot2::theme(
    panel.background = ggplot2::element_rect(fill = "white", colour = NA),
    plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
    plot.title       = ggplot2::element_text(size = 13, face = "bold"),
    plot.subtitle    = ggplot2::element_text(size = 9, lineheight = 1.05),
    plot.caption     = ggplot2::element_text(size = 9)
  )

  pts_disk <- subset(points, is.finite(u1) & is.finite(u2) & (u1^2 + u2^2 <= 1 + 1e-9))
  p_std <- ggplot(fields$Fstd, aes(u1, u2)) +
    geom_raster(aes(fill = p), interpolate = TRUE, na.rm = TRUE) +
    addC(transform(fields$Fstd, dx = var_name), "u1", "u2", br, is_prob_like) +
    geom_point(data = pts_disk, aes(u1, u2, colour = y), inherit.aes = FALSE, shape = 16, size = 1.9, alpha = 0.9) +
    coord_equal(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
    base_scale +
    pt_scale +
    labs(x = "u1 (whitened)", y = "u2", title = var_name, subtitle = subtitle_txt, caption = caption_txt) +
    theme_pub(11) +
    theme_fix

  pts_base <- subset(points, is.finite(b1) & is.finite(b2))
  p_base <- ggplot(fields$Fbase, aes(b1, b2)) +
    geom_raster(aes(fill = p), interpolate = TRUE, na.rm = TRUE) +
    addC(transform(fields$Fbase, dx = var_name), "b1", "b2", br, is_prob_like) +
    geom_point(data = pts_base, aes(b1, b2, colour = y), inherit.aes = FALSE, shape = 16, size = 1.8, alpha = 0.9) +
    coord_equal() +
    base_scale +
    pt_scale +
    labs(x = "b1", y = "b2", title = var_name, subtitle = subtitle_txt, caption = caption_txt) +
    theme_pub(11) +
    theme_fix

  pts_sq <- subset(points, is.finite(u1) & is.finite(u2))
  p_sq <- ggplot(fields$Fsq, aes(u1, u2)) +
    geom_raster(aes(fill = p), interpolate = TRUE, na.rm = TRUE) +
    addC(transform(fields$Fsq, dx = var_name), "u1", "u2", br, is_prob_like) +
    geom_point(data = pts_sq, aes(u1, u2, colour = y), inherit.aes = FALSE, shape = 16, size = 1.9, alpha = 0.9) +
    coord_equal(xlim = c(-1, 1), ylim = c(-1, 1), expand = FALSE) +
    base_scale +
    pt_scale +
    labs(x = "u1 (whitened)", y = "u2", title = var_name, subtitle = subtitle_txt, caption = caption_txt) +
    theme_pub(11) +
    theme_fix

  list(p_std = p_std, p_base = p_base, p_sq = p_sq)
}

# ==============================================================================
# 8. Core Kernel: Fit, Select, Metrics, Fields
# ==============================================================================

analyse_model_kernel <- function(d, yi, v_name, geom, Uobs, XR, is_two_col, K) {

  # --- 1. Fitting Logic ---
  # Define Additive (s + s) and Full Interaction (s + s + ti) formulas.
  k1_base <- min(MAX_K, length(unique(d$b1)) - 1L, floor(nrow(d) / 8))
  k2_base <- min(MAX_K, length(unique(d$b2)) - 1L, floor(nrow(d) / 8))
  k_int <- max(4L, min(8L, floor((k1_base + k2_base) / 3)))
  k1_add <- min(4L, k1_base)
  k2_add <- min(4L, k2_base)

  lhs <- if (is_two_col) "cbind(ys, yf)" else "y"
  f_add <- as.formula(sprintf("%s ~ s(b1, bs='ts', k=%d, m=2) + s(b2, bs='ts', k=%d, m=2)", lhs, k1_add, k2_add))
  f_full <- as.formula(sprintf("%s ~ s(b1, bs='ts', k=%d, m=2) + s(b2, bs='ts', k=%d, m=2) + ti(b1,b2, bs=c('ts','ts'), k=c(%d,%d), m=2)", lhs, k1_base, k2_base, k_int, k_int))

  fit_add <- tryCatch(R.utils::withTimeout(fit_gam_safe(f_add, yi$family, d, method = "REML"), timeout = FIT_TIMEOUT, onTimeout = "error"), error = function(e) NULL)
  if (is.null(fit_add)) {
    return(NULL)
  }
  fit_full <- tryCatch(R.utils::withTimeout(fit_gam_safe(f_full, yi$family, d, method = "REML"), timeout = FIT_TIMEOUT, onTimeout = "error"), error = function(e) NULL)

  # ML fits for Likelihood Ratio Testing
  fit_add_ml <- fit_gam_safe(f_add, yi$family, d, method = "ML")
  fit_full_ml <- fit_gam_safe(f_full, yi$family, d, method = "ML")

  de_add <- dev_expl_from_fit(fit_add)
  de_full <- dev_expl_from_fit(fit_full)
  a_add <- suppressWarnings(AIC(fit_add_ml))
  a_full <- suppressWarnings(AIC(fit_full_ml))
  dAIC <- if (is.finite(a_add) && is.finite(a_full)) a_add - a_full else NA_real_
  gain <- de_full - de_add
  rel_gain <- if (is.finite(de_add) && de_add > 0) gain / de_add else Inf

  # Selection Criteria: Prefer Full model if AIC improves and Deviance Explained gain is substantial.
  use_full <- (!is.null(fit_full) && is.finite(dAIC) && dAIC >= 2 && is.finite(gain) && gain > 0 && (gain >= 0.10 || rel_gain >= 0.50))
  fit_plot <- if (use_full) fit_full else fit_add
  fit_plot_ml <- if (use_full) fit_full_ml else fit_add_ml
  used_model <- if (use_full) sprintf("full s(b1)+s(b2)+ti(b1,b2) [k=(%d,%d;%d)]", k1_base, k2_base, k_int) else sprintf("additive s(b1)+s(b2) [k=(%d,%d)]", k1_add, k2_add)

  # --- 2. Metrics ---
  edf_used <- edf_total(fit_plot)
  kmin_used <- tryCatch(
    {
      kc <- suppressWarnings(mgcv::k.check(fit_plot))
      if (is.list(kc)) min(as.data.frame(kc$k.check)[, grep("kindex|k-index", colnames(kc$k.check))], na.rm = TRUE) else NA_real_
    },
    error = function(e) NA_real_
  )
  rho_s <- {
    y_fit <- tryCatch(if (!is.null(fit_plot$y)) fit_plot$y else model.response(model.frame(fit_plot)), error = function(e) NULL)
    if (is.matrix(y_fit) && ncol(y_fit) == 2) y_fit <- y_fit[, 1] / rowSums(y_fit)
    mu <- tryCatch(as.numeric(fitted(fit_plot)), error = function(e) NULL)
    if (is.null(y_fit) || is.null(mu)) NA_real_ else suppressWarnings(stats::cor(y_fit, mu, method = "spearman", use = "complete.obs"))
  }

  p_ml <- lrt_p_ml_nosel(f_add, f_full, yi$family, d)

  # Bootstrapped Interaction Test
  p_boot <- NA_real_
  if (DO_BOOT) {
    p_boot <- lrt_boot_gam(
      fit_null = fit_add_ml,
      form_null = f_add,
      form_alt  = f_full,
      family    = yi$family,
      data      = d,
      B         = BOOT_B,
      method    = "ML",
      select    = TRUE,
      is_two_col = is_two_col,
      trials     = if (is_two_col) d$M else NULL,
      is_ordinal = isTRUE(yi$is_ordinal),
      levels     = yi$levels %||% NULL
    )
  }
  p_for_fdr <- if (is.finite(p_boot)) p_boot else p_ml

  # R2 Calculation (McFadden for non-Gaussian)
  r2_label <- "McFadden"
  r2_val <- NA_real_
  try(
    {
      if (!is_two_col && inherits(yi$family, "family") && yi$family$family == "gaussian") {
        y_obs <- as.numeric(d$y)
        mu_hat <- as.numeric(fitted(fit_plot))
        if (length(y_obs) > 0) {
          r2_val <- 1 - sum((y_obs - mu_hat)^2) / sum((y_obs - mean(y_obs))^2)
          r2_label <- "R2"
        }
      } else {
        fit_null_ml <- fit_gam_safe(as.formula(paste(lhs, "~1")), yi$family, d, method = "ML")
        if (!is.null(fit_plot_ml) && !is.null(fit_null_ml)) r2_val <- 1 - as.numeric(logLik(fit_plot_ml) / logLik(fit_null_ml))
      }
    },
    silent = TRUE
  )

  # --- 3. OOF Predictions ---
  y_strat <- if (is_two_col) d$ys / pmax(d$ys + d$yf, 1L) else d$y
  is_bin_oof <- isTRUE(yi$is_binary) || (!is_two_col && length(unique(d$y)) == 2)
  fid <- make_folds_strat_general(y_strat, K = 5, seed = 42)
  p_oof <- rep(NA_real_, nrow(d))
  form_oof <- if (use_full) f_full else f_add

  for (k in 1:5) {
    tr <- fid != k
    te <- fid == k
    if (!any(tr) || !any(te)) next
    m <- fit_gam_safe(form_oof, yi$family, d[tr, , drop = FALSE], "REML")
    if (!is.null(m)) {
      pr <- predict(m, newdata = d[te, , drop = FALSE], type = "response")
      if (isTRUE(yi$is_ordinal) && is.matrix(pr) && !is.null(yi$levels)) {
        pr <- as.numeric(pr %*% yi$levels)     # expected value on original scale
      } else {
        pr <- as.numeric(pr)
      }
      p_oof[te] <- pr
    } else {
      p_oof[te] <- mean(y_strat[tr], na.rm = TRUE)
    }
  }

  # Init Vars
  oof_metric <- "C"
  oof_point <- NA
  oof_lo <- NA
  oof_hi <- NA
  win_oof <- NA
  skill_oof <- NA
  p_perm <- NA
  oof_auprc_point <- NA
  oof_auprc_lo <- NA
  oof_auprc_hi <- NA
  oof_auprg_point <- NA
  oof_auprg_lo <- NA
  oof_auprg_hi <- NA
  auc_residual <- NA
  auc_residual_lo <- NA
  auc_residual_hi <- NA
  auc_stacked <- NA
  auc_stacked_lo <- NA
  auc_stacked_hi <- NA
  p_delta_stacked <- NA
  calib_intercept <- NA
  calib_slope <- NA
  calib_brier <- NA
  calib_ece <- NA
  calib_nbins <- NA
  calib_curve <- NULL
  calib_line <- NULL
  calib_ntotal <- NA

  if (is_bin_oof) {
    # --- BINARY 0/1: Uses Bootstrapping + Stacking ---
    oof_metric <- "AUC"
    y_bin <- as.integer(d$y > 0)

    ci3 <- oof_triplet_ci(y_bin, p_oof, B = OOF_BOOT_B, seed = OOF_SEED)
    oof_point <- ci3$auc[1]
    oof_lo <- ci3$auc[2]
    oof_hi <- ci3$auc[3]
    oof_auprc_point <- ci3$prc[1]
    oof_auprc_lo <- ci3$prc[2]
    oof_auprc_hi <- ci3$prc[3]
    oof_auprg_point <- ci3$prg[1]
    oof_auprg_lo <- ci3$prg[2]
    oof_auprg_hi <- ci3$prg[3]

    # Fast Wilcoxon P-value
    wt <- suppressWarnings(stats::wilcox.test(p_oof ~ y_bin))
    p_perm <- wt$p.value

    # Stacking
    if (!is.null(XR)) {
      stack_res <- try(oof_prob_stacked_custom(y_bin, d[, c("b1", "b2")], XR, K = 5, seed = OOF_SEED), silent = TRUE)
      if (!inherits(stack_res, "try-error")) {
        cr <- boot_ci_balanced(y_bin, stack_res$Resid, auc_point, B = OOF_BOOT_B)
        auc_residual <- cr[1]
        auc_residual_lo <- cr[2]
        auc_residual_hi <- cr[3]
        cs <- boot_ci_balanced(y_bin, stack_res$Both, auc_point, B = OOF_BOOT_B)
        auc_stacked <- cs[1]
        auc_stacked_lo <- cs[2]
        auc_stacked_hi <- cs[3]
        i0 <- which(y_bin == 0)
        i1 <- which(y_bin == 1)
        difs <- replicate(OOF_BOOT_B, {
          ii <- c(sample(i0, length(i0), T), sample(i1, length(i1), T))
          auc_point(y_bin[ii], stack_res$Both[ii]) - auc_point(y_bin[ii], stack_res$Base[ii])
        })
        p_delta_stacked <- min(1, 2 * min(mean(difs >= 0), mean(difs <= 0)))
      }
    }
    cal <- calibration_metrics_binary(y_bin, p_oof)
    calib_intercept <- cal$intercept
    calib_slope <- cal$slope
    calib_brier <- cal$brier
    calib_ece <- cal$ece
    calib_nbins <- cal$n_bins
    calib_ntotal <- cal$n_total
    calib_curve <- cal$curve
    calib_line <- cal$line
  } else {
    # --- CONTINUOUS / COUNT / TRIALS ---

    # Fast P-value (Correlation)
    method_cor <- if (nrow(d) > 5000) "spearman" else "kendall"
    ct <- suppressWarnings(cor.test(y_strat, p_oof, method = method_cor))
    p_perm <- ct$p.value

    # C-Index (Point Estimate)
    bst <- boot_cindex_single(y_strat, p_oof, B = OOF_BOOT_B, seed = OOF_SEED)
    oof_point <- bst[1]
    oof_lo <- bst[2]
    oof_hi <- bst[3]
    win_oof <- bst[5]
  }

  # --- 4. Geometry & Plot Data ---
  u_pd <- principal_direction_from_fit(fit_plot, geom, Bobs = d)
  pd_angle <- atan2(u_pd[2], u_pd[1]) * 180 / pi
  pd_z <- pd_project_z(d$b1, d$b2, u_pd)$z
  
  pr_metrics <- list(pr80_20 = NA_real_, delta_pp = NA_real_)
  clin_res <- list(lift = NA_real_, nns_model = NA_real_, nns_base = NA_real_, npv = NA_real_, ppv = NA_real_)
  
  if (is_bin_oof) {
    pr_metrics <- pr80_20_from_score(as.integer(d$y > 0), pd_z)

    # Clinical Metrics (Enrichment, NNS, NPV) using OOF probabilities
    clin_res <- clinical_metrics_from_score(as.integer(d$y > 0), p_oof, frac = 0.20)
  }

  if (!is.null(geom$Xstd) && is.function(geom$Xstd$inv)) {
    UB <- as.data.frame(geom$Xstd$inv(geom$UD$grid$u1, geom$UD$grid$u2))
    names(UB) <- c("b1", "b2")
  } else {
    UB <- geom$gridX_sq[, c("b1", "b2")]
  }
  pU <- predict(fit_plot, newdata = UB, type = "response")
  pB <- predict(fit_plot, newdata = geom$gridB_full, type = "response")
  
  if (isTRUE(yi$is_ordinal) && is.matrix(pU) && !is.null(yi$levels)) {
    pU <- as.numeric(pU %*% yi$levels)
    pB <- as.numeric(pB %*% yi$levels)
  } else {
    pU <- as.numeric(pU)
    pB <- as.numeric(pB)
    if (!is_two_col && !is.null(yi$inv) && is.function(yi$inv)) {
      pU <- yi$inv(pU)
      pB <- yi$inv(pB)
    }
  }
  
  fields <- list(
    Fstd  = data.frame(u1 = geom$UD$grid$u1, u2 = geom$UD$grid$u2, p = pU),
    Fbase = data.frame(b1 = geom$gridB_full$b1, b2 = geom$gridB_full$b2, p = pB),
    Fsq   = data.frame(u1 = geom$UD$grid$u1, u2 = geom$UD$grid$u2, p = pU)
  )
  fields$Fstd$p[!geom$mask_sq] <- NA
  fields$Fbase$p[!geom$mask_hull] <- NA

  y_pts_full <- if (is_two_col) d$ys / d$M else if (!is.null(yi$inv)) yi$inv(d$y) else d$y
  if (!is.null(Uobs)) {
    mask_g <- in_whitened_geom(Uobs)
    pts <- data.frame(u1 = Uobs[mask_g, 1], u2 = Uobs[mask_g, 2], b1 = d$b1[mask_g], b2 = d$b2[mask_g], y = y_pts_full[mask_g])
  } else {
    pts <- data.frame(u1 = NA, u2 = NA, b1 = d$b1, b2 = d$b2, y = y_pts_full)
  }

  pdp_u1 <- pdp_1d_from_2d(
    fit_plot, as.data.frame(Uobs), geom, "u1",
    is_ordinal = isTRUE(yi$is_ordinal),
    levels     = yi$levels,
    ci_draws   = 300L,
    ci_obs     = 1200L,
    seed       = SEED_GLOBAL + 123L
  )
  
  pdp_u2 <- pdp_1d_from_2d(
    fit_plot, as.data.frame(Uobs), geom, "u2",
    is_ordinal = isTRUE(yi$is_ordinal),
    levels     = yi$levels,
    ci_draws   = 300L,
    ci_obs     = 1200L,
    seed       = SEED_GLOBAL + 456L
  )
  pdp_pd <- pdp_pd_from_2d(fit_plot, d, u_pd)
  if (!is_two_col && !is.null(yi$inv) && !isTRUE(yi$is_ordinal)) { 
    for (nm in c("p", "lo", "hi")) {
      if (nm %in% names(pdp_u1)) pdp_u1[[nm]] <- yi$inv(pdp_u1[[nm]])
      if (nm %in% names(pdp_u2)) pdp_u2[[nm]] <- yi$inv(pdp_u2[[nm]])
      if (nm %in% names(pdp_pd)) pdp_pd[[nm]] <- yi$inv(pdp_pd[[nm]])
    }
  }

  max_trials <- if (is_two_col && length(unique(d$ys + d$yf)) == 1) unique(d$ys + d$yf) else NA_integer_

  met_row <- data.table::data.table(
    var = v_name, n = nrow(d), K = K, family = if (inherits(yi$family, "family")) yi$family$family else "other",
    transform = yi$label, model_plotted = used_model,
    dev_expl_add = de_add, dev_expl_full = de_full, dev_expl_gain = gain, edf_used = edf_used, k_index_min_used = kmin_used, spearman_rho_used = rho_s,
    dAIC_add_to_full = dAIC, aic_weight_full = akaike_weight(dAIC), odds_full_to_add = akaike_odds(dAIC),
    p_for_fdr = p_for_fdr, p_ml_nosel = p_ml, p_boot_interaction = p_boot,
    r2_label = r2_label, r2_value = r2_val,
    oof_metric = oof_metric, oof_point = oof_point, oof_lo = oof_lo, oof_hi = oof_hi, p_oof = p_perm, win_oof = win_oof,
    oof_auprc_point = oof_auprc_point, oof_auprc_lo = oof_auprc_lo, oof_auprc_hi = oof_auprc_hi,
    oof_auprg_point = oof_auprg_point, oof_auprg_lo = oof_auprg_lo, oof_auprg_hi = oof_auprg_hi,
    auc_residual = auc_residual, auc_residual_lo = auc_residual_lo, auc_residual_hi = auc_residual_hi,
    auc_stacked = auc_stacked, auc_stacked_lo = auc_stacked_lo, auc_stacked_hi = auc_stacked_hi, p_delta_stacked = p_delta_stacked,
    pd_angle_deg = pd_angle, pd_vec_u1 = u_pd[1], pd_vec_u2 = u_pd[2], pr80_20 = pr_metrics$pr80_20, delta_pp = pr_metrics$delta_pp,
    calib_intercept = calib_intercept, calib_slope = calib_slope, calib_brier = calib_brier, calib_ece = calib_ece, calib_nbins = calib_nbins,

    # New columns for clinical utility
    clin_lift_20 = clin_res$lift, # Enrichment
    clin_nns_20 = clin_res$nns_model, # Workload (High Risk)
    clin_nns_base = clin_res$nns_base, # Workload (Baseline)
    clin_npv_20 = clin_res$npv, # Rule-Out Power
    clin_ppv_20   = clin_res$ppv, # Rule-In
    clin_sens_20 = clin_res$sens_20,
    clin_spec_20 = clin_res$spec_20,
    clin_sens_opt = clin_res$sens_opt,
    clin_spec_opt = clin_res$spec_opt,
    clin_ppv_opt  = clin_res$ppv_opt,
    clin_npv_opt  = clin_res$npv_opt,
    clin_thresh   = clin_res$thresh_opt
  )

  job <- list(
    key = v_name, label = v_name, fields = fields, points = pts, caption = yi$label, geom = geom,
    sub_stats = list(n = nrow(d), K = K, de_add = de_add, de_full = de_full, edf = edf_used, kmin = kmin_used, rho = rho_s, r2_lbl = r2_label, r2_val = r2_val, dAIC = dAIC, used_model = used_model, p_LRT = p_ml),
    oof_stats = list(name = oof_metric, point = oof_point, lo = oof_lo, hi = oof_hi, p = p_perm, auprc = oof_auprc_point, auprg = oof_auprg_point, pr80 = pr_metrics$pr80_20, dpp = pr_metrics$delta_pp, resid = auc_residual, resid_lo = auc_residual_lo, resid_hi = auc_residual_hi, stack = auc_stacked, stack_lo = auc_stacked_lo, stack_hi = auc_stacked_hi, p_delta = p_delta_stacked),
    calib_stats = list(int = calib_intercept, slope = calib_slope, brier = calib_brier, ece = calib_ece, bins = calib_nbins, n_total = calib_ntotal),
    calib_data = list(curve = calib_curve, line = calib_line),
    pd_stats = list(angle = pd_angle, u1 = u_pd[1], u2 = u_pd[2]),
    marginals = list(u1 = pdp_u1, u2 = pdp_u2, pd = pdp_pd),
    pd_points = data.frame(t = pd_z, y = y_pts_full),
    max_trials = max_trials,
    clin_stats = clin_res
  )

  list(metrics = met_row, job = job, fit = trim_gam_safe(fit_plot))
}

# ==============================================================================
# 9. Execution: Run All Variables
# ==============================================================================

num_cols <- names(DX)[vapply(DX, is.numeric, logical(1))]
vars <- if (exists("vars_keep") && !is.null(vars_keep)) vars_keep else setdiff(num_cols, c("b1", "b2"))
OVR <- get_model_overrides(vars, DX)

run_one_var <- function(v) {
  if (!(v %in% names(DX))) {
    return(NULL)
  }

  ov <- OVR[[v]]
  choice <- if (is.list(ov)) ov$choice else (ov %||% "auto")
  trials_fix <- if (is.list(ov)) ov$trials else NA_integer_
  trials_var <- if (is.list(ov)) ov$trials_var else NA_character_

  # Family setup based on overrides or auto-detection
  if (identical(choice, "binomial_trials")) {
    if (is.character(trials_var) && !is.na(trials_var) && nzchar(trials_var) && trials_var %in% names(DX)) {
      M_vec <- suppressWarnings(as.integer(round(DX[[trials_var]])))
      M_vec[!is.finite(M_vec) | M_vec < 1] <- NA_integer_
      y_cnt <- suppressWarnings(as.integer(round(pmin(pmax(DX[[v]] * M_vec, 0), M_vec))))
      yi <- list(
        y = y_cnt, family = binomial("logit"), label = sprintf("(binomial trials from %s)", trials_var),
        fwd = identity, inv = identity, trials = M_vec, two_col = TRUE
      )
    } else if (is.finite(trials_fix) && trials_fix > 1L) {
      M <- as.integer(trials_fix)
      y_raw_num <- suppressWarnings(as.numeric(DX[[v]]))
      ymax <- max(y_raw_num, na.rm = TRUE)
      if (is.finite(ymax) && ymax <= 1.1) {
        y_cnt <- suppressWarnings(as.integer(round(pmin(pmax(y_raw_num * M, 0), M))))
      } else {
        y_cnt <- suppressWarnings(as.integer(round(pmin(pmax(y_raw_num, 0), M))))
      }
      yi <- list(
        y = y_cnt, family = binomial("logit"), label = sprintf("(binomial trials=%d)", M),
        fwd = identity, inv = identity, trials = M, two_col = TRUE
      )
    } else {
      yi <- family_from_choice(choice, DX[[v]], var_name = v, trials = trials_fix)
    }
  } else {
    yi <- family_from_choice(choice, DX[[v]], var_name = v, trials = trials_fix)
    if (is.null(yi)) yi <- choose_family(DX[[v]], var_name = v)
  }

  # Branch: Nominal (One-vs-Rest loop)
  if (isTRUE(yi$is_nominal)) {
    out_list <- list()
    levels <- yi$levels
    for (lv in levels) {
      y_bin <- as.integer(DX[[v]] == lv)
      if (sum(y_bin, na.rm = TRUE) < K) next

      yi_sub <- list(y = y_bin, family = binomial("logit"), label = paste0(v, "==", lv, " (1vR)"), inv = identity)
      v_sub_name <- paste0(v, "==", lv)

      ok <- is.finite(y_bin) & is.finite(DX$b1) & is.finite(DX$b2)
      if (sum(ok) < 20) next
      d_sub <- data.frame(b1 = DX$b1[ok], b2 = DX$b2[ok], y = y_bin[ok])

      if (sum(d_sub$y > 0) < K) next

      u_sub <- if (!is.null(U_DX)) U_DX[ok, , drop = FALSE] else NULL
      
      # Subset XR to match the finite rows of this specific variable
      XR_sub <- if (!is.null(XR)) XR[ok, , drop = FALSE] else NULL
      
      res <- tryCatch(
        analyse_model_kernel(d_sub, yi, v, geom, u_sub, XR_sub, isTRUE(yi$two_col), K),
        error = function(e) {
          cat(sprintf("\n--- FAIL var: %s \nMessage: %s\n", v, e$message))
          NULL
        }
      )
      
      if (!is.null(res)) out_list[[v_sub_name]] <- res
    }
    return(out_list)
  }

  # Branch: Standard (Continuous/Ordinal/Trials)
  ok <- if (isTRUE(yi$two_col)) {
    y <- yi$y
    M <- if (length(yi$trials) > 1) yi$trials else rep(yi$trials, nrow(DX))
    is.finite(y) & is.finite(M) & is.finite(DX$b1) & is.finite(DX$b2) & M > 0
  } else {
    is.finite(yi$y) & is.finite(DX$b1) & is.finite(DX$b2)
  }
  if (sum(ok) < 20) {
    return(NULL)
  }

  d_sub <- data.frame(b1 = DX$b1[ok], b2 = DX$b2[ok])
  if (isTRUE(yi$two_col)) {
    M <- if (length(yi$trials) > 1) yi$trials[ok] else rep(yi$trials, sum(ok))
    d_sub$ys <- yi$y[ok]
    d_sub$yf <- M - yi$y[ok]
    d_sub$M <- M
  } else {
    d_sub$y <- yi$y[ok]
  }

  # For binary/binomial, ensure at least K positive and K negative cases.
  if (isTRUE(yi$is_binary) || (!is.null(yi$family) && yi$family$family == "binomial")) {
    val_check <- if (isTRUE(yi$two_col)) d_sub$ys else d_sub$y
    n_pos <- sum(val_check > 0, na.rm = TRUE)
    n_neg <- if (isTRUE(yi$two_col)) sum(val_check < d_sub$M, na.rm = TRUE) else sum(val_check == 0, na.rm = TRUE)
    if (n_pos < K || n_neg < K) {
      return(NULL)
    }
  }

  u_sub <- if (!is.null(U_DX)) U_DX[ok, , drop = FALSE] else NULL

  # Subset XR to match the finite rows of d_sub
  XR_sub <- if (!is.null(XR)) XR[ok, , drop = FALSE] else NULL
  
  # --- DIAGNOSTIC START ---
  # cat(sprintf("\n[DEBUG] Variable: %s\n", v))
  # cat(sprintf("  Dimensions of DX: %d rows\n", nrow(DX)))
  # cat(sprintf("  Length of 'ok' vector: %d\n", length(ok)))
  # cat(sprintf("  Sum of 'ok' (rows to keep): %d\n", sum(ok, na.rm = TRUE)))
  # cat(sprintf("  Dimensions of d_sub: %d rows\n", nrow(d_sub)))
  # 
  # if (is.null(XR)) {
  #   cat("  Global XR is NULL. (Stacking should be skipped)\n")
  # } else {
  #   cat(sprintf("  Global XR dimensions: %d x %d\n", nrow(XR), ncol(XR)))
  #   
  #   # Check if global XR aligns with global DX
  #   if (nrow(XR) != nrow(DX)) {
  #     cat("  [CRITICAL WARNING] Global XR rows != DX rows. Alignment broken at setup!\n")
  #   }
  # }
  # 
  # # Create XR_sub explicitly here to test dimensions
  # XR_debug <- if (!is.null(XR)) XR[ok, , drop = FALSE] else NULL
  # 
  # if (!is.null(XR_debug)) {
  #   cat(sprintf("  XR_sub dimensions: %d x %d\n", nrow(XR_debug), ncol(XR_debug)))
  #   
  #   if (nrow(XR_debug) != nrow(d_sub)) {
  #     cat("  [FAILURE CAUSE] d_sub and XR_sub row counts do not match!\n")
  #   } else {
  #     cat("  Status: d_sub and XR_sub align perfectly.\n")
  #   }
  # }
  # cat("------------------------------------------------\n")
  # --- DIAGNOSTIC END ---
  
  res <- tryCatch(analyse_model_kernel(d_sub, yi, v, geom, u_sub, XR_sub, isTRUE(yi$two_col), K),
    error = function(e) {
      cat("\n--- FAIL var:", v, "\n")
      cat("Message:", conditionMessage(e), "\n")
      NULL
    }
  )
  if (is.null(res)) {
    return(NULL)
  }
  list(res)
}

message("Pass 1: Model Fitting...")
res_list <- progressr::with_progress({
  p <- progressr::progressor(steps = length(vars))
  
  FUTURE_LAPPLY(
    vars,
    function(v) {
      
      # Explicitly load mgcv for extended families (OCAT, NB)
      library(mgcv)
      
      # Lock down threading to prevent conflicts in parallel workers
      if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1)
      data.table::setDTthreads(1)
      out <- tryCatch(
        {
          run_one_var(v)
        },
        error = function(e) {
          # Log error to console (visible in main session)
          message(sprintf("\n[Worker Error] Variable '%s' failed: %s", v, e$message))
          NULL 
        }
      )
      # Clean up RAM before next iteration
      gc(verbose = FALSE)

      p() 
      out
    },
    future.seed = TRUE
  )
})

# Flatten results and compute FDR
res_list <- unlist(res_list, recursive = FALSE)
res_list <- Filter(Negate(is.null), res_list)

MET <- if (length(res_list)) data.table::rbindlist(lapply(res_list, `[[`, "metrics"), fill = TRUE) else data.table::data.table()
JOBS <- lapply(res_list, `[[`, "job")
names(JOBS) <- sapply(JOBS, `[[`, "key")
DX_FITS <- lapply(res_list, `[[`, "fit")
names(DX_FITS) <- names(JOBS)

if (length(JOBS) > 0 && (is_dx_mode || is_cl_mode)) {
  
  # Determine prefix based on mode
  prefix <- if (is_cl_mode) "cluster" else "dx"
  message(sprintf("Saving %s mode artifacts...", prefix))
  
  # A. Extract metadata for the Fits file (Model objects + metadata)
  DX_PREV <- lapply(JOBS, function(j) mean(j$points$y, na.rm = TRUE))
  DX_USED_MODEL <- lapply(JOBS, function(j) j$sub_stats$used_model)
  
  saveRDS(
    list(fits = DX_FITS, prev = DX_PREV, used_model = DX_USED_MODEL),
    file.path(OUT_SUBDIR, paste0(prefix, "_gam_fits.rds"))
  )
  
  # B. Extract Fields for the Fields file (Surfaces)
  DX_FIELDS <- lapply(JOBS, `[[`, "fields")
  saveRDS(DX_FIELDS, file.path(OUT_SUBDIR, paste0(prefix, "_fields.rds")))
}

if (nrow(MET) > 0 && "p_oof" %in% names(MET)) MET[, q_oof := p.adjust(p_oof, method = "BH")]
write_csv(MET, file.path(OUT_SUBDIR, "behaviour_map_metrics.csv"))

# ==============================================================================
# 10. Plotting Worker & Execution
# ==============================================================================

plot_worker <- function(job, OUT_SUBDIR, MET_dt = NULL, OOF_PERM_B = 200) {
  nm <- job$key
  s <- job$sub_stats
  
  # --- 1. Data Scaling & Limits ---
  scale_factor <- 1.0
  if (!is.null(job$max_trials) && is.finite(job$max_trials) && job$max_trials > 1) {
    scale_factor <- as.numeric(job$max_trials)
  }
  
  scale_vec <- function(v) {
    if (is.numeric(v)) v * scale_factor else v
  }
  
  job$fields$Fstd$p <- scale_vec(job$fields$Fstd$p)
  job$fields$Fbase$p <- scale_vec(job$fields$Fbase$p)
  job$fields$Fsq$p <- scale_vec(job$fields$Fsq$p)
  job$points$y <- scale_vec(job$points$y)
  job$pd_points$y <- scale_vec(job$pd_points$y)
  
  for (m in names(job$marginals)) {
    df <- job$marginals[[m]]
    if (!is.null(df)) {
      df$p <- scale_vec(df$p)
      df$lo <- scale_vec(df$lo)
      df$hi <- scale_vec(df$hi)
      job$marginals[[m]] <- df
    }
  }
  
  lims_override <- NULL
  if (scale_factor > 1) {
    max_obs <- max(job$points$y, na.rm = TRUE)
    top_lim <- max(1.0, ceiling(max_obs))
    lims_override <- c(0, top_lim)
  }
  
  # --- 2. Subtitle Construction ---
  
  # A. Model Stats
  star <- if (is.finite(s$p_LRT) && s$p_LRT < 0.05) " ★" else ""
  r2_txt <- if (!is.null(s$r2_lbl) && is.finite(s$r2_val)) sprintf(" | %s=%.3f", s$r2_lbl, s$r2_val) else ""
  
  base_txt <- paste0(
    "n=", s$n, ", K=", s$K, "  |  model=", s$used_model, "\n",
    "Dev.expl(add/full)=", sprintf("%.3f/%.3f", s$de_add, s$de_full),
    " | EDF_used=", sprintf("%.1f", s$edf),
    " | k-index(min)=", ifelse(is.finite(s$kmin), sprintf("%.3f", s$kmin), "NA"),
    " | ρ_s=", ifelse(is.finite(s$rho), sprintf("%.3f", s$rho), "NA"),
    r2_txt, "\n",
    "ΔAIC(add→full)=", ifelse(is.finite(s$dAIC), sprintf("%.1f", s$dAIC), "NA"),
    " | LRT p=", ifelse(is.finite(s$p_LRT), formatC(s$p_LRT, digits = 2, format = "e"), "NA"),
    star
  )
  
  # B. Stacking Info
  o <- job$oof_stats
  if (!is.null(o$stack) && is.finite(o$stack)) {
    f_str <- if (is.finite(o$resid)) sprintf("Resid=%.3f(%.3f-%.3f)", o$resid, o$resid_lo %||% NA, o$resid_hi %||% NA) else ""
    s_str <- sprintf(" | Stack=%.3f(%.3f-%.3f)", o$stack, o$stack_lo %||% NA, o$stack_hi %||% NA)
    d_str <- if (is.finite(o$p_delta)) sprintf(" | Δp=%.3g", o$p_delta) else ""
    base_txt <- paste0(base_txt, "\n", trimws(paste0(f_str, s_str, d_str)))
  }
  
  # C. OOF P-values & Q-values
  q_val <- NA_real_
  if (!is.null(MET_dt) && "q_oof" %in% names(MET_dt)) {
    row_q <- MET_dt[var == nm, q_oof]
    if (length(row_q)) q_val <- row_q[1]
  }
  
  p_floor <- 1 / (OOF_PERM_B + 1)
  p_txt <- if (!is.finite(o$p)) {
    " · no CV signal"
  } else {
    if (o$p <= p_floor) sprintf(" | perm p<%.3g", p_floor) else sprintf(" | perm p=%.3g", o$p)
  }
  q_txt <- if (is.finite(q_val)) {
    if (q_val <= p_floor) sprintf(" | q<%.3f", p_floor) else sprintf(" | q=%.3f", q_val)
  } else {
    ""
  }
  
  # D. ML Metrics
  fmt_ci <- function(p, l, h) sprintf("%.3f (%.3f–%.3f)", p, l, h)
  oof_line <- paste0("OOF ", o$name, "=", fmt_ci(o$point, o$lo, o$hi), p_txt, q_txt)
  
  if (o$name == "AUC") {
    extras <- character(0)
    if (is.finite(o$auprc)) extras <- c(extras, sprintf("\nAUPRC=%.3f", o$auprc))
    if (is.finite(o$auprg)) extras <- c(extras, sprintf("AUPRG=%.3f", o$auprg))
    
    # [FIX] Clamp PR80 display to avoid "129630"
    if (is.finite(o$pr80)) {
      pr_val <- if(o$pr80 > 100) ">100" else sprintf("%.1f", o$pr80)
      extras <- c(extras, paste0("PR80/20=", pr_val))
    }
    if (is.finite(o$dpp))   extras <- c(extras, sprintf("Δpp=%.3f", o$dpp)) # Optional
    
    if (length(extras)) oof_line <- paste0(oof_line, " | ", paste(extras, collapse = " | "))
  }
  base_txt <- paste0(base_txt, "\n", oof_line, "\n")
  
  # E. Clinical Metrics
  if (!is.null(job$clin_stats) && is.finite(job$clin_stats$lift)) {
    cs <- job$clin_stats
    
    fmt_pct <- function(x) {
      if (!is.finite(x)) return("NA")
      val <- x * 100
      if (val < 0.1 && val > 0) "<0.1%" else if (val > 99.9 && val < 100) ">99.9%" else sprintf("%.0f%%", val)
    }
    
    # Line 1: Optimal Youden Index
    line_opt <- sprintf("Best Balance (p>%.2f): Sens %s | Spec %s | PPV %s | NPV %s", 
                        cs$thresh_opt, fmt_pct(cs$sens_opt), fmt_pct(cs$spec_opt), 
                        fmt_pct(cs$ppv_opt), fmt_pct(cs$npv_opt))
    
    # Line 2: Screening Utility (20/80)
    line_scr <- sprintf("Screening (Top 20%%):  Lift %.1fx  |  NNS %.0f (vs %.0f)", 
                        cs$lift, cs$nns_model, cs$nns_base)
    
    base_txt <- paste0(base_txt, "\n", line_opt, "\n", line_scr)
  }
  
  # F. Calibration
  c <- job$calib_stats
  if (!is.null(c) && is.finite(c$slope)) {
    nb <- if (isTRUE(c$bins > 0)) c$bins else 0
    cal_txt <- sprintf(
      "Calib: intercept=%.3f, slope=%.3f | Brier=%.3f | ECE=%.3f (bins=%d)",
      c$int, c$slope, c$brier, c$ece, nb
    )
    base_txt <- paste0(base_txt, "\n", cal_txt)
  }
  
  # G. PD Angle
  pd <- job$pd_stats
  if (is.finite(pd$angle)) {
    base_txt <- paste0(base_txt, "\n", sprintf("principal direction: θ=%+.1f° (u₁=%.3f, u₂=%.3f)", pd$angle, pd$u1, pd$u2))
  }
  
  # --- 3. Plot Generation ---
  
  P2d <- plots_for_cont(job$fields, job$geom, nm, base_txt, job$caption, job$points, lims_override = lims_override)
  save_plot_gg(file.path(OUT_SUBDIR, paste0("field_", nm, "_disk_full")), P2d$p_std, width = 7.2, height = 5.6)
  save_plot_gg(file.path(OUT_SUBDIR, paste0("field_", nm, "_base_full")), P2d$p_base, width = 7.2, height = 5.6)
  save_plot_gg(file.path(OUT_SUBDIR, paste0("field_", nm, "_square_full")), P2d$p_sq, width = 7.2, height = 5.6)
  
  theme_m <- theme_marginal()
  plot_scat <- function(df_pts, df_line, x_col, title_p) {
    p <- ggplot(df_pts, aes_string(x = x_col, y = "y")) +
      geom_point(shape = 16, size = 0.9, alpha = 0.6, colour = "black", stroke = 0) +
      geom_ribbon(data = df_line, aes_string(x = x_col, ymin = "lo", ymax = "hi"), inherit.aes = FALSE, alpha = 0.2) +
      geom_line(data = df_line, aes_string(x = x_col, y = "p"), inherit.aes = FALSE, linewidth = 1) +
      labs(title = title_p, caption = job$caption) +
      theme_m
    
    if (!is.null(lims_override)) {
      p <- p + coord_cartesian(ylim = lims_override)
    }
    p
  }
  
  p1 <- plot_scat(job$points, job$marginals$u1, "u1", paste("U1:", nm))
  save_plot_gg(file.path(OUT_SUBDIR, paste0("u1_scatter_", nm)), p1, width = 4, height = 4)
  
  p2 <- plot_scat(job$points, job$marginals$u2, "u2", paste("U2:", nm))
  save_plot_gg(file.path(OUT_SUBDIR, paste0("u2_scatter_", nm)), p2, width = 4, height = 4)
  
  ppd <- plot_scat(job$pd_points, job$marginals$pd, "t", paste("PD:", nm)) +
    labs(subtitle = sprintf("Angle: %.1f", pd$angle))
  save_plot_gg(file.path(OUT_SUBDIR, paste0("uPD_scatter_", nm)), ppd, width = 4, height = 4)
  
  if (!is.null(job$calib_data) && !is.null(job$calib_data$curve)) {
    df <- job$calib_data$curve
    line_df <- job$calib_data$line
    if (nrow(df) > 0) {
      df$p_hat <- pmin(pmax(df$p_hat, 0), 1)
      df$obs <- pmin(pmax(df$obs, 0), 1)
      cap_txt <- NULL
      if (isTRUE(c$n_total > 0) && isTRUE(c$bins > 0)) {
        mean_n <- c$n_total / c$bins
        cap_txt <- sprintf("%d bins of ~%.1f patients each (n=%d)", c$bins, mean_n, round(c$n_total))
      }
      pc <- ggplot(df, aes(p_hat, obs)) +
        geom_abline(intercept = 0, slope = 1, linetype = 2, colour = "grey70", linewidth = 0.5)
      if (!is.null(line_df) && nrow(line_df) > 0) {
        pc <- pc + geom_ribbon(data = line_df, aes(x = p_hat, ymin = lo, ymax = hi), inherit.aes = FALSE, alpha = 0.15) +
          geom_line(data = line_df, aes(p_hat, p), inherit.aes = FALSE, linewidth = 0.9)
      }
      pc <- pc + geom_point(shape = 16, size = 1.8, alpha = 0.7) +
        coord_equal(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
        labs(x = "Predicted probability", y = "Observed event rate", title = paste("Calibration —", nm), subtitle = cal_txt, caption = cap_txt) +
        theme_pub(11) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), legend.position = "none")
      save_plot_gg(file.path(OUT_SUBDIR, paste0("calibration_", nm)), pc, width = 5.5, height = 5)
    }
  }
}

# --- Pass 2: Serialized Parallel Plotting ---

if (nrow(MET) > 0) {
  MET_dt <- as.data.table(MET)
  if ("p_oof" %in% names(MET_dt)) {
    MET_dt[, q_oof := p.adjust(p_oof, method = "BH")]
  }
  if ("p_for_fdr" %in% names(MET)) {
    MET[, q_model := p.adjust(p_for_fdr, method = "BH")]
  }
} else {
  MET_dt <- NULL
}

message("Pass 2: Plotting...")

if (length(JOBS) > 0) {
  # Serialize plot jobs to disk to manage memory overhead during parallel plotting
  job_dir <- file.path(OUTPUTS_DIR, OUT_SUBDIR, "plot_jobs")
  dir.create(job_dir, recursive = TRUE, showWarnings = FALSE)

  message(sprintf("Serializing %d plot jobs to disk...", length(JOBS)))

  job_files <- character(length(JOBS))
  names(job_files) <- names(JOBS)
  for (nm in names(JOBS)) {
    fpath <- file.path(job_dir, paste0("job_", nm, ".rds"))
    saveRDS(JOBS[[nm]], fpath)
    job_files[nm] <- fpath
  }

  rm(JOBS, res_list)
  gc()

  progressr::with_progress({
    p <- progressr::progressor(steps = length(job_files))

    invisible(FUTURE_LAPPLY(job_files, function(fpath) {
      if (requireNamespace("RhpcBLASctl", quietly = TRUE)) RhpcBLASctl::blas_set_num_threads(1)
      data.table::setDTthreads(1)

      job_data <- readRDS(fpath)
      tryCatch(plot_worker(job_data, OUT_SUBDIR, MET_dt), error = function(e) {
        message(sprintf("Error plotting %s: %s", job_data$key, e$message))
      })

      p()
      NULL
    }, future.seed = TRUE))
  })
}

if (is_dx_mode) saveRDS(DX_FITS, file.path(OUT_SUBDIR, "dx_gam_fits.rds"))

message("Done.")

if (is_dx_mode || is_cl_mode) {
  lbl <- if (is_cl_mode) "Cluster" else "Diagnostic"
  message(sprintf("%s run complete. Artifacts saved to %s/", lbl, OUT_SUBDIR))
}

if (is_first_run) {
  is_first_run <- FALSE
  behaviour_csv <- restore_csv
  OUT_SUBDIR <- restore_subdir
  message("First run complete. Re-source dimension_contplot to run desired analyses.")
}