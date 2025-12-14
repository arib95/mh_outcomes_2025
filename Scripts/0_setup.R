# =============================== Unified Setup ================================
suppressPackageStartupMessages({
  # Base + utils
  library(stats); library(utils); library(grDevices); library(Matrix)
  # Data
  library(readr); library(dplyr); library(tidyr); library(data.table); library(zoo);
  library(tibble);  library(yaml)
  # Maths / ML
  library(cluster); library(RANN); library(FNN); library(RSpectra); library(dbscan)
  library(glmnet); library(MASS); library(expm); library(clue); library(rrcov)
  library(Matrix)
  # Graphs / viz
  library(ggplot2); library(scales); library(ggrepel); library(webshot2);
  library(patchwork)
  # Domain / extras
  library(igraph); library(aricode); library(vegan); library(princurve); library(FactoMineR)
  library(reticulate); library(R.utils); library(pROC); library(PRROC)
  library(robustbase); library(isotone); library(grid); library(metR);
  
  # parallel
  library(withr)
  library(future)
  library(future.apply)
  library(parallelly)
  library(progressr)
  library(RhpcBLASctl)
})

# =============================== Config =======================================
cfg <- yaml::read_yaml("config/config.yml")

# ==================== Build mgcv with OpenMP in Apple Silicon =================
# There is experimental support for building mgcv with OpenMP, inc. for Apple 
# Silicon Macs. A useful Stack Overflow thread can be found at:
# https://stackoverflow.com/questions/78590558/how-can-i-enable-openmp-in-mac-os-for-mgcv-r-package
library(mgcv)

# =============================== Helpers =====================================
# coalesce nulls
`%||%` <- function(a, b) if (!is.null(a) && length(a) > 0) a else b
# relpath check
.is_relpath <- function(p) !grepl("^(/|[A-Za-z]:)", p) && dirname(p) == "."
.is_abs_path <- function(p) grepl("^(/|~|[A-Za-z]:)", p)

# filename safety
safe_file <- function(s) gsub("[^A-Za-z0-9_]+","_", s)
# timing helper
.timeit <- function(expr){
  t0 <- proc.time()[[3]]
  val <- eval.parent(substitute(expr))
  list(value = val, elapsed = proc.time()[[3]] - t0)
}
# binary guards
.is_binary <- function(x) is.numeric(x) && all(x %in% c(0,1,NA))
.check_binary_df <- function(df) stopifnot(all(vapply(df, .is_binary, logical(1))))

# =============================== Outputs & I/O ================================
OUTPUTS_DIR <- cfg$project$outputs_dir %||% "out"
dir.create(OUTPUTS_DIR, recursive = TRUE, showWarnings = FALSE)
PNG_DPI     <- as.integer(cfg$io$images$png_dpi %||% 300)
SAVE_IMAGES <- isTRUE(cfg$images$save %||% TRUE)
SAVE_PLOT_RDS <- isTRUE(cfg$images$save_rds %||% TRUE)

log_msg <- function(..., .cfg) {
  if (isTRUE(.cfg$verbose %||% FALSE)) message(sprintf(...))
}

write_csv <- function(x, file, ...) {
  if (!.is_abs_path(file)) file <- file.path(OUTPUTS_DIR, file)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  readr::write_csv(x, file, ...)
}

saveRDS <- function(object, file, ...) {
  if (!.is_abs_path(file)) file <- file.path(OUTPUTS_DIR, file)
  dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  base::saveRDS(object, file = file, ...)
}

read_csv <- function(file, ...) {
  if (!.is_abs_path(file)) file <- file.path(OUTPUTS_DIR, file)
  utils::read.csv(file, ...)
}

readRDS <- function(file, ...) {
  if (!.is_abs_path(file)) file <- file.path(OUTPUTS_DIR, file)
  base::readRDS(file, ...)
}

# Helper to strip the massive environment from a ggplot object
clean_ggplot_env <- function(p) {
  if (inherits(p, "ggplot")) {
    p$plot_env <- rlang::new_environment()  # Replace heavy environment with an empty one
  }
  p
}

save_plot_gg <- function(name, plot, width, height, dpi = PNG_DPI) {
  # If neither images nor RDS are needed, do nothing
  if (!SAVE_IMAGES && !SAVE_PLOT_RDS) return(invisible())
  
  base_path <- file.path(OUTPUTS_DIR, name)
  png_path  <- paste0(base_path, ".png")
  pdf_path  <- paste0(base_path, ".pdf")
  rds_path  <- paste0(base_path, ".rds")
  
  if (SAVE_IMAGES) {
    if (isTRUE(USE_RAGG) && requireNamespace("ragg", quietly = TRUE)) {
      ggplot2::ggsave(
        png_path, plot = plot,
        width = width, height = height, dpi = dpi,
        device = ragg::agg_png
      )
    } else {
      ggplot2::ggsave(
        png_path, plot = plot,
        width = width, height = height, dpi = dpi
      )
    }
    
    ggplot2::ggsave(
      pdf_path, plot = plot,
      width = width, height = height,
      device = grDevices::cairo_pdf
    )
  }
  
  if (SAVE_PLOT_RDS) {
    # This prevents the RDS from dragging 100MB of global variables (geom, XR) with it.
    plot_clean <- clean_ggplot_env(plot)
    
    rds_file <- paste0(base_path, ".rds")
    saveRDS(plot_clean, file = rds_file)
  }
  
  invisible(plot)
}

save_plot_base <- function(name, width, height, dpi = PNG_DPI, expr) {
  if (!SAVE_IMAGES) return(invisible())
  png_path <- file.path(OUTPUTS_DIR, paste0(name, ".png"))
  pdf_path <- file.path(OUTPUTS_DIR, paste0(name, ".pdf"))
  grDevices::png(png_path, width = width, height = height, units = "in", res = dpi)
  on.exit(grDevices::dev.off(), add = TRUE); eval.parent(substitute(expr)); grDevices::dev.off()
  grDevices::pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
  on.exit(grDevices::dev.off(), add = TRUE); eval.parent(substitute(expr)); grDevices::dev.off()
}

# =============================== RNG & Threads ================================
RNG_KIND <- cfg$rng$rng_kind %||% "L'Ecuyer-CMRG"
do.call(RNGkind, as.list(strsplit(RNG_KIND, ",")[[1]]))
SEED            <- as.integer(cfg$rng$seed %||% 42L);
SEED_GLOBAL     <- as.integer(cfg$rng$seed %||% 42L)
BUNDLE_SEED     <- as.integer(cfg$rng$bundle_seed %||% SEED_GLOBAL)
SEED_PRED       <- as.integer(cfg$rng$seed_pred %||% SEED_GLOBAL)
SEED_JITTER     <- as.integer(cfg$rng$seed_jitter %||% SEED_GLOBAL)
SEED_BOOT       <- as.integer(cfg$rng$seed_boot %||% SEED_GLOBAL)
OOF_SEED      <- as.integer(cfg$contplot$seed %||% SEED_GLOBAL)

set.seed(SEED_GLOBAL)

# Parallel options
BAM_THREADS <- cfg$compute$bam_threads %||% "all_but_one"
if (is.character(BAM_THREADS) && BAM_THREADS == "all_but_one") {
  BAM_THREADS <- max(1L, parallel::detectCores(logical = TRUE) - 1L)
}
options(mc.cores = BAM_THREADS)

# BLAS/OMP env threads
set_env_threads <- isTRUE(cfg$compute$set_env_threads %||% TRUE)
omp_threads  <- as.integer(cfg$compute$omp_threads %||% 1L)
blas_threads <- as.integer(cfg$compute$blas_threads %||% 1L)
if (set_env_threads) {
  Sys.setenv(
    OMP_NUM_THREADS         = as.character(omp_threads),
    OPENBLAS_NUM_THREADS    = as.character(blas_threads),
    MKL_NUM_THREADS         = as.character(blas_threads),
    VECLIB_MAXIMUM_THREADS  = as.character(blas_threads),
    BLAS_NUM_THREADS        = as.character(blas_threads)
  )
}

if (identical(cfg$weights$ncores_par, "all_but_one")) {
  NCORES_PAR <- max(1L, parallel::detectCores(logical=FALSE) - 1L)
} else {
  NCORES_PAR <- as.integer(cfg$weights$ncores_par %||%
                             max(1L, parallel::detectCores(logical=FALSE) - 1L))
}

# ==================================== Plots ===================================
# Plotting
cfg$palette <- list( #ds000030: lipari+vik; ds005237: lapaz+cork
  engine    = "scico",
  name      = "lipari",      # sequential 
  name_div  = "vik",      # diverging
  direction = 1
)

scale_shared_div <- function(name = NULL, limits = NULL){
  if (requireNamespace("scico", quietly = TRUE)) {
    scico::scale_fill_scico(
      palette   = cfg$palette$name_div %||% "vik",
      direction = cfg$palette$direction %||% 1,
      limits    = limits, oob = scales::squish, name = name
    )
  } else {
    ggplot2::scale_fill_gradient2(limits = limits, oob = scales::squish, name = name)
  }
}

scale_prob_fill <- function(limits = NULL, name = NULL) {
  eng <- tolower(cfg$palette$engine %||% "scico")
  dir <- cfg$palette$direction %||% 1
  pal <- cfg$palette$name %||% "lipari"
  if (eng == "scico" && requireNamespace("scico", quietly = TRUE)) {
    scico::scale_fill_scico(
      palette = pal, direction = dir,
      limits = limits, oob = scales::squish, name = name
    )
  } else if (eng == "colorspace" && requireNamespace("colorspace", quietly = TRUE)) {
    colorspace::scale_fill_continuous_sequential(
      palette = pal, rev = (dir == -1),
      limits = limits, name = name
    )
  } else if (eng == "brewer") {
    ggplot2::scale_fill_distiller(
      palette = pal, type = "seq",
      direction = if (dir == 1) 1 else -1,
      limits = limits, oob = scales::squish, name = name
    )
  } else if (eng == "paletteer" && requireNamespace("paletteer", quietly = TRUE)) {
    paletteer::scale_fill_paletteer_c(pal,
                                      direction = dir,
                                      limits = limits, name = name
    )
  } else if (eng == "manual" && !is.null(cfg$palette$colours)) {
    ggplot2::scale_fill_gradientn(
      colours = cfg$palette$colours,
      limits = limits, name = name
    )
  } else {
    ggplot2::scale_fill_gradient(
      low = "#f7fbff", high = "#08306b",
      limits = limits, name = name
    )
  }
}

cluster_colours <- function(levels_vec){
  if (requireNamespace("scico", quietly = TRUE)) {
    vals <- scico::scico(length(levels_vec), palette = cfg$palette$name %||% "vik",
                         direction = cfg$palette$direction %||% 1)
    stats::setNames(vals, levels_vec)
  } else {
    cols <- grDevices::hcl.colors(length(levels_vec), "Dark 3")
    stats::setNames(cols, levels_vec)
  }
}

# ---- Theme & palette ----------------------------------------------------------
theme_pub <- function(base_size = 11, base_family = "sans") {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      panel.grid       = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.background  = ggplot2::element_rect(fill = "white", colour = NA),
      strip.background = ggplot2::element_blank(),
      legend.title     = ggplot2::element_text(size = base_size),
      legend.text      = ggplot2::element_text(size = base_size - 1)
    )
}

theme_marginal <- function(base_size = 11) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      # panel
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major = ggplot2::element_line(
        colour   = "grey85",
        linetype = "dotted",
        linewidth = 0.25
      ),
      panel.grid.minor = ggplot2::element_line(
        colour   = "grey92",
        linetype = "dotted",
        linewidth = 0.25
      ),
      # axes
      axis.line        = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.ticks       = ggplot2::element_line(colour = "black", linewidth = 0.4),
      axis.ticks.length = grid::unit(3, "pt"),
      axis.text        = ggplot2::element_text(colour = "black"),
      axis.title       = ggplot2::element_text(colour = "black"),
      # titles
      plot.title    = ggplot2::element_text(
        face = "bold", hjust = 0, size = ggplot2::rel(1.4)
      ),
      plot.subtitle = ggplot2::element_text(hjust = 0),
      # misc
      legend.position = "none",
      plot.margin = grid::unit(c(5.5, 5.5, 5.5, 5.5), "pt")
    )
}

# ============================= Parallel Computing =============================
options(
  future.rng.onMisuse   = "warn",
  future.globals.maxSize = 1024^3
)

options(progressr.enable = TRUE)
progressr::handlers(global = TRUE)

NWORKERS <- as.integer(NCORES_PAR %||% 1L)

ok <- try({
  if (NWORKERS <= 1L) {
    future::plan(future::sequential)
  } else if (.Platform$OS.type != "windows" && parallelly::supportsMulticore()) {
    future::plan(future::multicore, workers = NWORKERS)
  } else {
    future::plan(future::multisession, workers = NWORKERS)
  }
}, silent = TRUE)

if (inherits(ok, "try-error")) {
  message("[future] fallback to sequential: ", conditionMessage(attr(ok, "condition")))
  future::plan(future::sequential)
}

.drop_future_args <- function(dots){
  if (is.null(names(dots))) return(dots)
  dots[!grepl("^future\\.", names(dots), perl = TRUE)]
}

FUTURE_LAPPLY <- function(X, FUN, ...) {
  dots <- list(...)
  ok <- try(do.call(future.apply::future_lapply, c(list(X, FUN), dots)), silent = TRUE)
  if (!inherits(ok, "try-error")) return(ok)
  message("[parallel] future_lapply failed → base::lapply: ",
          conditionMessage(attr(ok, "condition")))
  do.call(base::lapply, c(list(X, FUN), .drop_future_args(dots)))
}

FUTURE_SAPPLY <- function(X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE) {
  dots <- list(...)
  ok <- try(do.call(future.apply::future_sapply,
                    c(list(X, FUN, simplify = simplify, USE.NAMES = USE.NAMES), dots)),
            silent = TRUE)
  if (!inherits(ok, "try-error")) return(ok)
  message("[parallel] future_sapply failed → base::sapply: ",
          conditionMessage(attr(ok, "condition")))
  do.call(base::sapply,
          c(list(X, FUN, simplify = simplify, USE.NAMES = USE.NAMES), .drop_future_args(dots)))
}

FUTURE_MAPPLY <- function(FUN, ..., MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE) {
  dots <- list(...)
  ok <- try(
    do.call(future.apply::future_mapply,
            c(list(FUN), dots,
              list(MoreArgs = MoreArgs, SIMPLIFY = SIMPLIFY, USE.NAMES = USE.NAMES))),
    silent = TRUE
  )
  if (!inherits(ok, "try-error")) return(ok)
  message("[parallel] future_mapply failed, falling back to mapply: ",
          conditionMessage(attr(ok, "condition")))
  do.call(base::mapply,
          c(list(FUN), .drop_future_args(dots),
            list(MoreArgs = MoreArgs, SIMPLIFY = SIMPLIFY, USE.NAMES = USE.NAMES)))
}

# ========================= Resave plot as RDS =================================
resave_plot <- function(name,
                        width,
                        height,
                        dpi   = PNG_DPI,
                        types = c("png", "pdf"),
                        suffix = NULL) {
  p <- readRDS(file.path(paste0(name, ".rds")))
  
  if (!is.null(suffix) && nzchar(suffix)) {
    name_out <- paste0(name, "_", suffix)
  } else {
    name_out <- paste0(name, "_", width, "x", height)
  }
  
  base_path <- file.path(OUTPUTS_DIR, name_out)
  
  if ("png" %in% types) {
    ggplot2::ggsave(
      paste0(base_path, ".png"),
      plot   = p,
      width  = width,
      height = height,
      dpi    = dpi
    )
  }
  
  if ("pdf" %in% types) {
    ggplot2::ggsave(
      paste0(base_path, ".pdf"),
      plot   = p,
      width  = width,
      height = height,
      device = grDevices::cairo_pdf
    )
  }
  
  invisible(base_path)
}

# resave_plot("plotXYZ", width = 4, height = 6, suffix = "tall/wide/square")

# =============================== Toggles & Params =============================
# Global toggles
DO_PLOTS        <- isTRUE(cfg$toggles$do_plots %||% TRUE)
DO_DIAGNOSTICS  <- isTRUE(cfg$toggles$do_diagnostics %||% TRUE)
DO_SURFACE      <- isTRUE(cfg$toggles$do_surface %||% TRUE)
DO_SWEEP        <- isTRUE(cfg$toggles$do_sweep %||% TRUE)
RUN_PRINSURF    <- isTRUE(cfg$toggles$run_prinsurf %||% TRUE)
USE_SURFACE_CV  <- isTRUE(cfg$toggles$use_surface_cv %||% TRUE)
DO_DELTAQ       <- isTRUE(cfg$toggles$do_corr_trim %||% TRUE)
DO_CORR_TRIM    <- isTRUE(cfg$toggles$do_corr_trim %||% FALSE)

# Embedding (clusters)
KNN_K           <- as.integer(cfg$embedding_cluster$knn_k %||% 9L)
KNN_VARIANT     <- tolower(cfg$embedding_cluster$knn_variant %||% "mutual")   # mutual | union
LOCAL_SCALE     <- isTRUE(cfg$embedding_cluster$local_scale %||% TRUE)
MULT_WEIGHT     <- isTRUE(cfg$embedding_cluster$multiplicity_weight %||% TRUE)

# Embedding (psychometric)
M_STAR_FIXED    <- as.integer(cfg$embedding$m_star_fixed %||% 2L)
K_ID_LO_HI      <- unlist(cfg$embedding$id_k_range %||% list(8L, 20L))
PC_MAX          <- as.integer(cfg$embedding$pc_max %||% 6L)
M_DEFAULT       <- as.integer(cfg$embedding$m_default %||% 4L)

# Community detection
COMMUNITY_ALGO    <- tolower(cfg$graph$community_algo %||% "leiden")
LEIDEN_OBJECTIVE  <- cfg$graph$leiden_objective %||% "modularity"
LEIDEN_GAMMA      <- as.numeric(cfg$graph$leiden_gamma %||% 0.75)
LEIDEN_ITERS      <- as.integer(cfg$graph$leiden_iters %||% -1L)
MIN_CLUSTER_SIZE  <- as.integer(cfg$graph$min_cluster_size %||% 2L)        # on unique profiles
MIN_CLUSTER_WEIGHT<- as.integer(cfg$graph$min_cluster_weight %||% 8L)      # expanded weight

# Surface
K_BASIS_UV       <- as.integer(cfg$surface$k_basis_uv %||% 40L)
SURFACE_CV_FOLDS <- as.integer(cfg$surface$surface_cv_folds %||% 3L)

# Bundle-fibre
BF_NGR_UV       <- as.integer(cfg$bundle_fibre$grid_res_uv %||% 35L)
BF_K_BASE       <- as.integer(cfg$bundle_fibre$k_base %||% 100L)
BF_K_RANK_CAP   <- as.integer(cfg$bundle_fibre$k_rank_cap %||% 4L)
BF_MP_SHRINK    <- as.numeric(cfg$bundle_fibre$mp_shrink %||% 0.90)
BF_MIN_N_NEIGH  <- as.integer(cfg$bundle_fibre$min_n_neigh %||% 50L)
BF_DENS_K       <- as.integer(cfg$bundle_fibre$dens_k %||% 30L)

# Nulls / resampling
N_PERM          <- as.integer(cfg$nulls$perm_B %||% 200L)
DEG_REPS        <- as.integer(cfg$nulls$degree_reps %||% N_PERM)
NULL_SCOPE      <- tolower(cfg$nulls$degree_scope %||% "kept")           # kept|full
N_BOOT          <- as.integer(cfg$cv$bootstrap_B %||% 200L)
SIL_REPS        <- as.integer(cfg$nulls$sil_reps %||% 200L)
COL_SHUFFLE_REPS<- as.integer(cfg$nulls$col_shuffle_reps %||% 200L)
CV_FOLDS        <- as.integer(cfg$cv$folds %||% 5L)
SIG_METHOD      <- tolower(cfg$sig$method %||% "perm")
SIG_FWER        <- as.numeric(cfg$sig$fwer %||% 0.95)
PRED_BOOT       <- as.integer(cfg$sig$pred_boot %||% 200L)

# Intrinsic dimension
CORE_BAND       <- unlist(cfg$dedup$core_band_limits %||% list(0.15, 0.85))
CORE_KNN_K      <- as.integer(cfg$dedup$core_knn_k %||% max(5L, min(20L, KNN_K)))

# Features / guards
INCLUDE_NODIAG  <- isTRUE(cfg$features$include_nodiag %||% TRUE)
DX_MIN_PREV     <- as.numeric(cfg$features$prevalence_min %||% 0.00)
DX_MAX_PREV     <- as.numeric(cfg$features$prevalence_max %||% 0.99)
MIN_CASES_TOTAL <- as.integer(cfg$features$dx_min_cases %||% 10L)
MIN_CASES_IN    <- as.integer(cfg$features$min_in_cases %||% 5L)
MIN_CASES_OUT   <- as.integer(cfg$features$min_out_cases %||% 5L)
SHOW_NS_LABELS  <- isTRUE(cfg$features$show_ns_labels %||% TRUE)

# Targets / features (additional)
PREF_TARGET         <- toupper(cfg$targets$pref_target %||% "ANY_DX")
CORR_THRESH         <- as.numeric(cfg$features$corr_trim_threshold %||% 0.95)
DX_DENY_NOS         <- isTRUE(cfg$features$deny_nos %||% FALSE)
DX_PREV_MIN         <- as.numeric(cfg$features$prevalence_min %||% 0.00)
DX_CASES_MIN        <- as.integer(cfg$features$dx_min_cases %||% 10L)
N_TOP_PER_DX        <- as.integer(cfg$features$top_items_per_dx %||% 80L)
RARE_LEVEL_MIN_PROP <- as.numeric(cfg$features$rare_level_min_prop %||% 0.01)

# Pillars (A,B,C)
ALPHA_FDR       <- as.numeric(cfg$majors$alpha_fdr %||% 0.05)
MIN_PREV_IN_CL  <- as.numeric(cfg$majors$min_prev_in_cluster %||% 0.00)
MIN_OR          <- as.numeric(cfg$majors$min_odds_ratio %||% 2.0)
ALPHA_LOCALIZE  <- as.numeric(cfg$majors$alpha_localize %||% 0.05)
AUC_MIN         <- as.numeric(cfg$majors$auc_min %||% 0.70)
PREV_MIN        <- as.numeric(cfg$majors$prev_min %||% 0.03)
NCASE_MIN       <- as.integer(cfg$majors$ncase_min %||% 5L)
NIN_MIN         <- as.integer(cfg$majors$nin_min %||% 5L)
NOUT_MIN        <- as.integer(cfg$majors$nout_min %||% 5L)
MAJORS_BOOT_B   <- as.integer(cfg$majors$majors_boot_B %||% 200L)
PILLAR_B_REPS   <- as.integer(cfg$majors$pillar_B_reps %||% 200L)

# Clusters
K_GRID          <- as.numeric(cfg$sweep$k_grid       %||% c(5,6,7,8,9,10,11,12,13,14,15))
VARIANT_GRID    <- tolower(   cfg$sweep$variant_grid %||% c("union", "mutual"))
SWEEP_BOOT_B    <- as.integer(cfg$sweep$sweep_boot_B %||% 200L)

# Weights
W_MIN           <- as.numeric(cfg$weights$w_min          %||% 0.00)
W_STEP_GRID     <- as.numeric(cfg$weights$grid           %||% c(0.95, 0.90, 0.75, 0.5, 0.25, 0.10, 0.05, 0.01, 0.00))
W_BATCH_K       <- as.integer(cfg$weights$batch_k        %||% 3L)
W_BATCH_FACTOR  <- as.numeric(cfg$weights$batch_factor   %||% 0.75)
W_MAX_ITERS     <- as.numeric(cfg$weights$max_iters)     %||% NA_integer_
N_ROWS_SUB      <- as.integer(cfg$weights$n_rows_sub)    %||% NULL
FIX_REP_SUBSET  <- isTRUE(    cfg$weights$fix_rep_subset %||% TRUE)

# Dedup / diagnostics
DO_DEDUP        <- isTRUE(    cfg$dedup$enabled               %||% TRUE)
EPS_DEDUP       <- as.numeric(cfg$dedup$eps                 %||% NA_real_)
WRITE_DEDUP_CSV <- isTRUE(    cfg$dedup$write_csv             %||% TRUE)
CORE_BAND       <- unlist(    cfg$dedup$core_band_limits      %||% list(0.20, 0.70))
CORE_KNN_K      <- as.integer(cfg$dedup$core_knn_k            %||% 10L)
SIG_Q           <- as.numeric(cfg$diagnostics$q_alpha         %||% 0.01)
K_FIBRE_CAP     <- as.integer(cfg$diagnostics$fibre_k_cap     %||% 3L)
MIN_N_DX        <- as.integer(cfg$diagnostics$min_n_dx        %||% 25L)

# Core sweep / optimisation
PCL_K_MIN        <- as.integer(cfg$prob_clusters$k_min        %||% 2L)
PCL_K_MAX        <- as.integer(cfg$prob_clusters$k_max        %||% 4L)
PCL_LAMBDA_GRID  <- as.numeric(cfg$prob_clusters$lambda_grid  %||% c(0, 0.05, 0.1, 0.2, 0.4, 0.8))
PCL_KNN_BASE     <- as.integer(cfg$prob_clusters$knn_base     %||% 12L)
PCL_LOCAL_SCALE  <- isTRUE(    cfg$prob_clusters$local_scale  %||% TRUE)   # avoid clash with embedding_cluster$local_scale
PCL_N_START      <- as.integer(cfg$prob_clusters$n_start      %||% 3L)
PCL_N_ITER       <- as.integer(cfg$prob_clusters$n_iter       %||% 700L)
PCL_TOL          <- as.numeric(cfg$prob_clusters$tol          %||% 1e-5)
PCL_EPS          <- as.numeric(cfg$prob_clusters$eps          %||% 1e-10)
PCL_MIN_POS_DEF  <- as.integer(cfg$prob_clusters$min_pos_def  %||% 10L)
PCL_MIN_NEG_DEF  <- as.integer(cfg$prob_clusters$min_neg_def  %||% 10L)

# Usage / regularisation
PCL_MASS_FLOOR_FRAC <- as.numeric(cfg$prob_clusters$mass_floor_frac %||% 0.05)
PCL_BALANCE_EVERY   <- as.integer(cfg$prob_clusters$balance_every   %||% 8L)
PCL_RESEED_MAX      <- as.integer(cfg$prob_clusters$reseed_max      %||% 3L)
PCL_ENTROPY_PUSH    <- as.numeric(cfg$prob_clusters$entropy_push    %||% 1e-3)

# Selection weights / penalties
PCL_W_Q       <- as.numeric(cfg$prob_clusters$w_q        %||% 0.45)
PCL_W_SNN     <- as.numeric(cfg$prob_clusters$w_snn      %||% 0.50)
PCL_W_REC     <- as.numeric(cfg$prob_clusters$w_recon    %||% 0.20)
PCL_W_SIL     <- as.numeric(cfg$prob_clusters$w_sil      %||% 0.15)
PCL_PEN_EMPTY <- as.numeric(cfg$prob_clusters$pen_empty  %||% 0.60)
PCL_PEN_IMB   <- as.numeric(cfg$prob_clusters$pen_imb    %||% 0.15)
PCL_PEN_RED   <- as.numeric(cfg$prob_clusters$pen_red    %||% 0.30)

# Leiden
LEIDEN_MEMBERSHIP_CSV <- cfg$io$inputs$leiden_membership_csv %||% "cluster_membership_all_participants.csv"
USE_LEIDEN_K          <- isTRUE(cfg$prob_clusters$force_leiden_k     %||% FALSE)

# Plots
USE_RAGG    <- isTRUE(cfg$plots$use_ragg %||% FALSE)
MM_UNITS    <- isTRUE(cfg$plots$mm_units %||% FALSE)
GRID_N_B    <- as.integer(cfg$plots$nu_unit %||% 400L)
GRID_N_U    <- as.integer(cfg$plots$n_base_grid %||% 400L)
CONTOUR_PROBS <- as.numeric(cfg$plots$contour_probs %||% c(0.15, 0.25, 0.40, 0.60, 0.80))
CALIB_SHOWPTS <- isTRUE(cfg$plots$calib_showpts %||% TRUE)
CL_DRAW_HULLS <- isTRUE(cfg$plots$cl_draw_hulls %||% TRUE)

# k-NN Stability Plot
STAB_K_RANGE  <- as.integer(cfg$knn_stab$K_range %||% c(10,11,12,13,14,15,16,17,18,19,20))
STAB_SD_GRID  <- as.numeric(cfg$knn_stab$SD_grid %||% c(0, 0.05, 0.10, 0.15, 0.20))
STAB_REPS     <- as.integer(cfg$knn_stab$reps %||% 600L)

# Outcome plots
# behaviour_csv <- "data/cluster_membership_all_participants.csv"
# OUT_SUBDIR    <- "clusters"
behaviour_csv <- "data/wide_diagnoses.csv"
OUT_SUBDIR    <- "wide_diagnoses"
# behaviour_csv <- "data/health.csv"
# OUT_SUBDIR    <- "health"
vars_keep     <- cfg$behaviour$vars_keep                 %||% NULL
write_plots   <- isTRUE(cfg$behaviour$write_plots        %||% DO_PLOTS)

MAX_K         <- as.integer(cfg$behaviour$max_k          %||% 8L)
MIN_POS       <- as.integer(cfg$behaviour$min_pos        %||% 15L)
FIT_TIMEOUT   <- as.numeric(cfg$behaviour$timeout        %||% 60.0)

FDR_METHOD    <- toupper(cfg$sig$fdr_method              %||% "BH")   # BH | BY
DO_BOOT       <- isTRUE(cfg$behaviour$boot$enabled       %||% TRUE)
BOOT_B        <- as.integer(cfg$behaviour$boot$B         %||% 199L)

# OOF significance config (AUC/C with paired bootstrap) ----
OOF_BOOT_B    <- as.integer(cfg$contplot$bootstrap_B          %||% 400L)
OOF_MAX_PAIRS <- as.numeric(cfg$contplot$max_pairs            %||% 2e6)
OOF_PERM_B    <- as.integer(cfg$contplot$perm_B               %||% 200L)

# =============================== Inputs =======================================
PSY_CSV  <- cfg$io$inputs$psychometric_matrix_csv %||% "data/psychometric_matrix.csv"
DIAG_CSV <- cfg$io$inputs$wide_diagnoses_csv     %||% "data/wide_diagnoses.csv"
