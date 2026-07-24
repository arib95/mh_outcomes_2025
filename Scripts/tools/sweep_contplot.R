# sweep_dimension_contplot.R
# Runs dimension_contplot.R across every CSV in the data/ folder that
# setup.R already references via the default `behaviour_csv`.

# --------- config toggles (optional) ----------
SKIP_IF_METRICS_EXIST <- TRUE   # set TRUE to skip if metrics already present
STOP_ON_ERROR         <- TRUE   # set TRUE to halt on first failure

# --------- bootstrap project context ----------
.sweep_file <- tryCatch(normalizePath(sys.frame(1)$ofile, mustWork = TRUE), error = function(e) NA_character_)
if (!is.character(.sweep_file) || length(.sweep_file) != 1L || is.na(.sweep_file)) {
  stop("Could not determine the location of sweep_contplot.R; source it by filename.")
}
PROJECT_ROOT <- dirname(dirname(.sweep_file))
setwd(PROJECT_ROOT)

SETUP_CFG <- list()
source(file.path(PROJECT_ROOT, "0_setup.R"), chdir = FALSE)

# Where to sweep: derive the data directory from the configured behaviour_csv
data_dir <- file.path(PROJECT_ROOT, dirname(behaviour_csv))

# Allow overriding the data_dir from the command line (optional)
args <- commandArgs(TRUE)
if (length(args) >= 1L && nzchar(args[1])) data_dir <- args[1]

# Enumerate CSVs (ignore .xlsx and others)
csv_files <- sort(list.files(data_dir, pattern = "\\.csv$", full.names = TRUE, ignore.case = TRUE)) 
csv_files <- csv_files[basename(csv_files) != "psychometric_matrix.csv"]

if (!length(csv_files)) {
  stop(sprintf("No CSV files found under: %s", normalizePath(data_dir, mustWork = FALSE)))
}

message(sprintf("Sweeping %d CSV files under: %s", length(csv_files), normalizePath(data_dir, mustWork = FALSE)))

# Quick header reader
.header_cols <- function(f) {
  tryCatch(
    names(data.table::fread(f, sep = ";", dec = ",", nrows = 0L, showProgress = FALSE)),
    error = function(e) character(0)
  )
}

# --------- main sweep ----------
for (f in csv_files) {
  out_stub <- tools::file_path_sans_ext(basename(f))
  out_dir  <- file.path(OUTPUTS_DIR, out_stub)
  
  # Minimal sanity: must have participant_id
  cols0 <- .header_cols(f)
  if (!("participant_id" %in% cols0)) {
    message(sprintf("[-] SKIP (no participant_id): %s", basename(f)))
    next
  }
  
  # Optional skip if a prior run exists
  if (isTRUE(SKIP_IF_METRICS_EXIST)) {
    metrics_fp <- file.path(out_dir, "behaviour_map_metrics.csv")
    if (file.exists(metrics_fp)) {
      message(sprintf("[=] SKIP (metrics exist): %s → %s", basename(f), metrics_fp))
      next
    }
  }
  
  # Ensure per-file output folder exists (model_overrides.csv will live here)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Inject the two variables dimension_contplot.R reads
  behaviour_csv <- f
  OUT_SUBDIR    <- out_stub
  
  message(sprintf("[+] RUN  %s  → OUT_SUBDIR=%s", basename(f), OUT_SUBDIR))
  
  # Execute one full run. Keep the global env to reuse setup objects.
  ok <- try({
    source(file.path(PROJECT_ROOT, "3_dimension_contplot.R"), chdir = FALSE, local = FALSE)
  }, silent = !STOP_ON_ERROR)
  
  if (inherits(ok, "try-error")) {
    msg <- conditionMessage(attr(ok, "condition"))
    message(sprintf("[!] ERROR in %s: %s", basename(f), msg))
    if (isTRUE(STOP_ON_ERROR)) stop(ok)
  } else {
    message(sprintf("[✓] DONE  %s", basename(f)))
  }
}

message("Sweep complete.")
