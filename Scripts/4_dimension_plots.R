# ==============================================================================
#                               dimension_plots.R
# ==============================================================================

# ------------------------------------------------------------------------------
# Diagnosis artefact readers
# ------------------------------------------------------------------------------

read_dx_fit_bundle <- function(dx_fit_rds) {
  DXFIT_all <- readRDS(dx_fit_rds)
  if (is.null(DXFIT_all)) return(NULL)
  # expected list(fits=..., prev=...)
  DXFIT_all
}

read_dx_grids <- function(dx_fields_rds, keep_dx = NULL) {
  DXGRIDS_raw <- readRDS(dx_fields_rds)
  if (is.null(DXGRIDS_raw)) return(NULL)
  
  if (!is.null(keep_dx)) {
    DXGRIDS_raw <- DXGRIDS_raw[intersect(names(DXGRIDS_raw), keep_dx)]
  }
  
  bind_field <- function(field) {
    pieces <- lapply(names(DXGRIDS_raw), function(dx) {
      f <- DXGRIDS_raw[[dx]][[field]]
      if (is.null(f)) return(NULL)
      if (!"dx" %in% names(f)) f$dx <- dx
      f
    })
    dplyr::bind_rows(pieces)
  }
  
  list(
    raw   = DXGRIDS_raw,
    Fstd  = bind_field("Fstd"),
    Fbase = bind_field("Fbase"),
    Fsq   = bind_field("Fsq"),
    Fsig  = bind_field("Fsig") %||% NULL
  )
}

read_dx_metrics <- function(dx_met_csv, keep_dx = NULL) {
  M <- read_csv(dx_met_csv)
  if (is.null(M)) return(NULL)
  if ("var" %in% names(M) && !"dx" %in% names(M)) M <- dplyr::rename(M, dx = var)
  if (!"dx" %in% names(M)) return(M)
  if (!is.null(keep_dx)) M <- dplyr::filter(M, dx %in% keep_dx)
  M
}

# ------------------------------------------------------------------------------
# Diagnosis secondary outputs
# ------------------------------------------------------------------------------

export_dx_metrics <- function(M, out_dir = OUTPUTS_DIR) {
  if (is.null(M) || !nrow(M)) {
    msgf("[dx-metrics] No metrics to export.")
    return(invisible(NULL))
  }
  
  # Normalise var -> dx if needed
  if ("var" %in% names(M) && !"dx" %in% names(M)) {
    M <- dplyr::rename(M, dx = var)
  }
  
  # --------------------------------------------------------------------------
  # Case 1: already in expected long AUC schema
  # --------------------------------------------------------------------------
  if (all(c("dx", "model", "AUC", "lo", "hi") %in% names(M))) {
    readr::write_csv(M, file.path(out_dir, "AUC_bootstrap_CI.csv"))
    msgf("[dx-metrics] Wrote AUC_bootstrap_CI.csv (dx/model/AUC/lo/hi).")
    return(invisible(TRUE))
  }
  
  # --------------------------------------------------------------------------
  # Case 2: contplot OOF schema
  #   expects: dx, oof_metric, oof_point, oof_lo, oof_hi
  #   optional: model_plotted
  # --------------------------------------------------------------------------
  if (all(c("dx", "oof_metric", "oof_point", "oof_lo", "oof_hi") %in% names(M))) {
    
    M2 <- M
    
    if (!"model" %in% names(M2) && "model_plotted" %in% names(M2)) {
      M2 <- dplyr::rename(M2, model = model_plotted)
    }
    if (!"model" %in% names(M2)) {
      M2$model <- "GAM"
    }
    
    # AUC export
    A <- M2 |>
      dplyr::filter(oof_metric == "AUC") |>
      dplyr::transmute(
        dx, model,
        AUC = oof_point,
        lo  = oof_lo,
        hi  = oof_hi
      )
    
    if (nrow(A)) {
      readr::write_csv(A, file.path(out_dir, "AUC_bootstrap_CI.csv"))
      msgf("[dx-metrics] Wrote AUC_bootstrap_CI.csv from contplot oof_* fields.")
    } else {
      msgf("[dx-metrics] No AUC rows found in contplot metrics.")
    }
    
    # Optional future-proofing if you later add these metrics back to contplot
    PRC <- M2 |>
      dplyr::filter(oof_metric == "AUPRC") |>
      dplyr::transmute(dx, model, AUPRC = oof_point, lo = oof_lo, hi = oof_hi)
    
    if (nrow(PRC)) {
      readr::write_csv(PRC, file.path(out_dir, "AUPRC_bootstrap_CI.csv"))
      msgf("[dx-metrics] Wrote AUPRC_bootstrap_CI.csv from contplot metrics.")
    }
    
    PRG <- M2 |>
      dplyr::filter(oof_metric == "AUPRG") |>
      dplyr::transmute(dx, model, AUPRG = oof_point, lo = oof_lo, hi = oof_hi)
    
    if (nrow(PRG)) {
      readr::write_csv(PRG, file.path(out_dir, "AUPRG_bootstrap_CI.csv"))
      msgf("[dx-metrics] Wrote AUPRG_bootstrap_CI.csv from contplot metrics.")
    }
    
    return(invisible(TRUE))
  }
  
  msgf("[dx-metrics] Metrics schema not recognised; no exports written.")
  invisible(FALSE)
}

lift_map_plot <- function(DxW_A, order = c("original", "prevalence", "cluster"),
                          adjust_p = FALSE, title = "Co-occurrence (lift)") {
  order <- match.arg(order)
  if (is.null(DxW_A) || !ncol(DxW_A)) return(NULL)
  
  M <- as.matrix(DxW_A > 0)
  dx_all <- colnames(M)
  n <- nrow(M)
  p <- colMeans(M)
  
  keep <- which(is.finite(p) & p > 0 & p < 1)
  if (!length(keep)) {
    return(ggplot() +
             labs(title = paste0(title, " — no eligible diagnoses"), x = NULL, y = NULL) +
             theme_minimal(12) +
             theme(panel.grid = element_blank()))
  }
  M <- M[, keep, drop = FALSE]
  dx <- dx_all[keep]
  p <- p[keep]
  
  Pij <- (t(M) %*% M) / n
  L <- Pij / (p %o% p)
  
  ord <- switch(order,
                prevalence = order(p, decreasing = TRUE),
                cluster = {
                  S <- suppressWarnings(cor(M))
                  S[!is.finite(S)] <- 0
                  hclust(as.dist(1 - S), "average")$order
                },
                original = seq_along(dx)
  )
  dx <- dx[ord]
  p <- p[ord]
  L <- L[ord, ord, drop = FALSE]
  
  L_pairs <- L
  diag(L_pairs) <- NA_real_
  L_pairs[lower.tri(L_pairs, TRUE)] <- NA_real_
  
  U <- as.data.frame(as.table(L_pairs), responseName = "lift") |>
    dplyr::rename(r = Var1, c = Var2) |>
    dplyr::filter(!is.na(lift))
  
  # NOTE: We intentionally remove Fisher tests here to avoid any perception
  # of "inferential refitting" beyond simple descriptive matrices.
  U$star <- ""
  
  cap_hi <- stats::quantile(U$lift, 0.99, na.rm = TRUE)
  cap_lo <- min(U$lift, na.rm = TRUE)
  
  U <- U |>
    dplyr::mutate(
      lift_clip = pmin(pmax(lift, cap_lo), cap_hi)
    )
  
  D <- tibble::tibble(r = dx, c = dx, lift = 1, star = "")
  U$r <- factor(U$r, levels = dx)
  U$c <- factor(U$c, levels = dx)
  D$r <- factor(D$r, levels = dx)
  D$c <- factor(D$c, levels = dx)
  
  ggplot() +
    geom_tile(data = D, aes(c, r, fill = lift), colour = NA) +
    geom_tile(data = U, aes(c, r, fill = lift), colour = NA) +
    scale_fill_viridis_c(name = "lift") +
    coord_equal(expand = FALSE) +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1),
      axis.ticks = element_blank(),
      plot.title = element_text(face = "bold", size = 16, hjust = .5)
    )
}

# ------------------------------------------------------------------------------
# Cluster artefact readers
# ------------------------------------------------------------------------------

read_cluster_grids <- function(cluster_fields_rds) {
  readRDS(cluster_fields_rds)
}

read_cluster_archetypes <- function(
    out_dir = OUTPUTS_DIR,
    prefer_wide = TRUE
) {
  if (prefer_wide) {
    fp <- "cluster_archetypes_dx_Kbest.csv"
    Vw <- read_csv(fp)
    if (is.null(Vw)) return(NULL)
    stopifnot("dx" %in% names(Vw))
    rn <- Vw$dx
    V  <- as.matrix(Vw[setdiff(names(Vw), "dx")])
    rownames(V) <- rn
    ord <- order(as.integer(gsub("^C", "", colnames(V))))
    V[, ord, drop = FALSE]
  } else {
    Lfp   <- file.path(op, "cluster_archetypes_dx_by_k.csv")
    Bfp   <- file.path(op, "best_fit_summary.csv")
    L     <- read_csv(Lfp)
    BEST  <- read_csv(Bfp)
    if (is.null(L) || is.null(BEST)) return(NULL)
    Kb    <- as.integer(BEST$K_best[1])
    lb    <- as.numeric(BEST$lambda_best[1])
    tol   <- 1e-6
    L2    <- dplyr::filter(L, K == Kb, abs(lambda - lb) <= tol)
    if (!nrow(L2)) return(NULL)
    W     <- tidyr::pivot_wider(L2, names_from = "cluster", values_from = "loading")
    rn    <- W$dx
    V     <- as.matrix(W[setdiff(names(W), "dx")])
    rownames(V) <- rn
    ord <- order(as.integer(gsub("^C", "", colnames(V))))
    V[, ord, drop = FALSE]
  }
}

wrap_dx <- function(x, width = 28) vapply(x, \(s) paste(strwrap(s, width), collapse = "\n"), "")

plot_cluster_dx_heatmap <- function(V, cap = 2) {
  if (is.null(V)) return(NULL)
  V <- as.matrix(V)
  dx_levels <- rownames(V) %||% seq_len(nrow(V))
  
  col_sums <- colSums(V, na.rm = TRUE)
  P <- sweep(V, 2, ifelse(col_sums > 0, col_sums, 1), "/")
  
  expected_dx <- rowMeans(P, na.rm = TRUE)
  expected_dx[expected_dx <= 0] <- 1e-12
  
  LIFT <- sweep(P, 1, expected_dx, "/")
  LOG2 <- log2(pmax(LIFT, 1e-12))
  VAL  <- pmin(pmax(LOG2, -cap), cap)
  
  H <- as.data.frame(VAL) |>
    tibble::rownames_to_column("dx") |>
    tidyr::pivot_longer(-dx, names_to = "cluster", values_to = "val") |>
    dplyr::mutate(dx = factor(dx, levels = dx_levels))
  
  labs_wrapped <- setNames(wrap_dx(dx_levels), dx_levels)
  
  ggplot(H, aes(cluster, dx, fill = val)) +
    geom_tile() +
    scale_shared_div(name = "log2(lift)", limits = c(-cap, cap)) +
    scale_y_discrete(labels = labs_wrapped) +
    labs(x = "cluster", y = NULL, title = "Diagnosis archetypes per soft cluster") +
    theme_pub(15) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_text(lineheight = 0.9)
    )
}

plots_for_clusters <- function(CL_LIST, geom) {
  if (length(CL_LIST) == 0) return(NULL)
  
  # Helper to plot one cluster field with its specific subtitle
  plot_one <- function(obj, field_type = "Fstd", x_col, y_col, xlim, ylim) {
    df <- obj$fields[[field_type]]
    
    ggplot(df, aes(x = .data[[x_col]], y = .data[[y_col]])) +
      geom_raster(aes(fill = p), interpolate = TRUE) +
      metR::geom_contour2(
        aes(z = p, label = after_stat(sprintf("%.0f%%", ..level..*100))),
        colour = "white", linewidth = 0.2, breaks = c(0.25, 0.50, 0.75, 0.90)
      ) +
      # Use Prob Scale (0-1)
      scale_prob_fill(limits = c(0, 1), name = "Prob") +
      coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
      labs(
        title = paste("Cluster", obj$key), 
        subtitle = obj$sub_base
      ) +
      theme_pub(10) +
      theme(
        axis.title = element_blank(),
        plot.subtitle = element_text(size = 8, lineheight = 1.1, colour = "grey30"),
        legend.position = "none" # Hide legends for grid view to save space
      )
  }
  
  # Generate lists of plots
  pl_std  <- list()
  pl_base <- list()
  pl_sq   <- list()
  
  # Determine limits for Base
  b_rng <- range(c(geom$gridB_full$b1, geom$gridB_full$b2), na.rm=TRUE)
  
  for (nm in names(CL_LIST)) {
    obj <- CL_LIST[[nm]]
    
    pl_std[[nm]]  <- plot_one(obj, "Fstd",  "u1", "u2", c(-1,1), c(-1,1))
    pl_base[[nm]] <- plot_one(obj, "Fbase", "b1", "b2", NULL, NULL)
    pl_sq[[nm]]   <- plot_one(obj, "Fsq",   "u1", "u2", c(-1,1), c(-1,1))
  }
  
  # Combine into grids (using patchwork if available, else gridExtra)
  compose <- function(plist) {
    if (requireNamespace("patchwork", quietly = TRUE)) {
      patchwork::wrap_plots(plist, ncol = 3) # Adjust ncol as needed
    } else {
      do.call(gridExtra::grid.arrange, c(plist, ncol = 3))
    }
  }
  
  list(
    p_std  = compose(pl_std),
    p_base = compose(pl_base),
    p_sq   = compose(pl_sq)
  )
}

# ------------------------------------------------------------------------------
# Model comparison table
# ------------------------------------------------------------------------------

write_dx_auc_comparison_table <- function(
    DxW_A,
    metrics_csv  = "wide_diagnoses/behaviour_map_metrics.csv",
    out_path     = "Diagnosis_Model_AUC_Table"
){
  # 1. Read the High-Quality GAM Metrics from dimension_contplot
  if (!file.exists(paste0("out/", metrics_csv))) {
    warning("High-quality metrics file not found: ", metrics_csv, 
            "\nRun dimension_contplot.R first to generate GAM estimates.")
    return(NULL)
  }
  
  # Read and filter for diagnoses present in the current set
  MET <- read_csv(metrics_csv) %>%
    dplyr::filter(var %in% names(DxW_A))
  
  if (nrow(MET) == 0) {
    warning("No matching diagnoses found in metrics file.")
    return(NULL)
  }
  
  # 2. Add Counts (n_pos / N)
  CNT <- tibble::tibble(
    var   = names(DxW_A),
    n_pos = colSums(DxW_A > 0, na.rm = TRUE),
    N     = nrow(DxW_A)
  )
  
  DF <- MET %>% 
    dplyr::left_join(CNT, by = "var")
  
  # 3. Formatting Helpers
  fmt3 <- function(x) ifelse(is.finite(x), sprintf("%.3f", as.numeric(x)), "NA")
  
  # Helper to format CI: "0.800 [0.750, 0.850]"
  # Note: behaviour_map_metrics.csv has columns like oof_point, oof_lo, oof_hi
  fmt_ci <- function(pt, lo, hi) {
    ifelse(is.finite(pt), 
           sprintf("%.3f [%.3f, %.3f]", pt, lo, hi), 
           "NA")
  }
  
  # 4. Construct the Table using GAM metrics
  # Note: oof_point is the Base GAM AUC. 
  # auc_residual / auc_stacked were added in the recent fix to dimension_contplot.
  
  OUT <- DF %>%
    dplyr::transmute(
      Diagnosis = var,
      `n+/N (%)` = sprintf("%d/%d (%.1f%%)", n_pos, N, 100 * n_pos / N),
      
      # Base Model (The Surface)
      `Base AUC (GAM)` = fmt_ci(oof_point, oof_lo, oof_hi),
      
      # Residual / Stacking Metrics (if available)
      `Fibre AUC (GLM)`    = fmt3(auc_residual),
      `Combined AUC (Stack)` = fmt3(auc_stacked),
      
      # Delta (Combined - Base)
      # Note: We calculated p_delta_stacked in contplot, but didn't calculate CI for delta there yet.
      # showing point estimate delta for now.
      `dAUC (Stack - Base)` = ifelse(is.finite(auc_stacked) & is.finite(oof_point),
                                     sprintf("%+.3f", auc_stacked - oof_point), "NA"),
      
      `p (Stack vs Base)` = fmt3(p_delta_stacked)
    ) %>%
    dplyr::arrange(Diagnosis)
  
  # 5. Write
  write_csv(OUT, paste0(out_path, ".csv"))
  
  # Optional: Write PDF/TIFF using gt if available
  if (requireNamespace("gt", quietly = TRUE)) {
    try({
      tab <- gt::gt(OUT) %>%
        gt::tab_header(title = "Diagnosis Prediction: Base Surface vs. Stacked Residuals") %>%
        gt::fmt_number(columns = -Diagnosis, decimals = 3)
      gt::gtsave(tab, paste0(OUTPUTS_DIR, "/", out_path, ".pdf"))
    }, silent = TRUE)
  }
  
  message("Comparison table written to: ", out_path, ".csv (Source: GAMs from contplot)")
}

# ------------------------------------------------------------------------------
# RUN
# ------------------------------------------------------------------------------

# Read contplot outputs
dx_fit_rds     <- file.path("wide_diagnoses/dx_gam_fits.rds")       # optional for metadata
dx_fields_rds  <- file.path("wide_diagnoses/dx_fields.rds")         # required for field plots
dx_met_csv     <- file.path("wide_diagnoses/behaviour_map_metrics.csv")

DXFIT_all <- read_dx_fit_bundle(dx_fit_rds)
if (is.null(DXFIT_all)) {
  msgf("[dx] dx_gam_fits.rds not found.")
}

DXGRIDS <- read_dx_grids(dx_fields_rds, keep_dx = keep_dx)
if (is.null(DXGRIDS)) {
  msgf("[dx] dx_fields.rds not found; diagnosis surface plots will be skipped.")
}

M_dx <- read_dx_metrics(dx_met_csv, keep_dx = keep_dx)
if (is.null(M_dx)) {
  msgf("[dx] behaviour_map_metrics.csv not found; metrics exports/tables will be skipped.")
}

# Metrics exports from contplot outputs
export_dx_metrics(M_dx, out_dir = OUTPUTS_DIR)

# Co-occurrence lift map
p_co <- lift_map_plot(DxW_A, order = "cluster", adjust_p = FALSE)
if (!is.null(p_co)) {
  save_plot_gg("FIG_cooccurrence_lift_uppertri_diagprev", p_co, width = 8.0, height = 7.0)
}

# Cluster archetypes
V_best <- read_cluster_archetypes()
p_heat <- plot_cluster_dx_heatmap(V_best)
if (!is.null(p_heat)) {
  save_plot_gg("FIG_cluster_by_dx_heatmap_Kbest", p_heat, width = 9.0, height = 7.0)
} else {
  msgf("[clusters] No archetype CSVs found; skipping heatmap.")
}

# Model comparison table
write_dx_auc_comparison_table(DxW_A)

msgf("[dimension_plots] Run complete.")