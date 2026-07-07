# ==============================================================================
# Residual signal probe
# ==============================================================================
# Needs to be run after running dimension_contplot as "wide_diagnoses" so
# its XR is kept in environment.

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a

char_to_num_local <- function(v) {
  if (is.numeric(v)) return(v)
  if (is.factor(v) || is.logical(v)) v <- as.character(v)
  if (!is.character(v)) return(v)
  s <- gsub("\\s+", "", v)
  s <- sub("%$", "", s)
  v2 <- suppressWarnings(as.numeric(gsub(",", ".", s)))
  non_blank <- !is.na(s) & nzchar(s)
  prop_num <- if (any(non_blank)) sum(is.finite(v2[non_blank])) / sum(non_blank) else 0
  if (prop_num >= 0.80) v2 else v
}

probe_residual_signal <- function(target_var, dx_data, xr_matrix, geom_u, top_n = 20) {
  
  if (!"participant_id" %in% names(dx_data)) stop("dx_data must contain participant_id.")
  if (is.null(rownames(xr_matrix))) stop("xr_matrix must have participant_id rownames.")
  if (is.null(rownames(geom_u))) stop("geom_u must have participant_id rownames.")

  common_id <- Reduce(intersect, list(
    as.character(dx_data$participant_id),
    rownames(xr_matrix),
    rownames(geom_u)
  ))
  if (!length(common_id)) stop("No overlapping participant_id values across DX, XR, and geom_u.")

  dx_data <- dx_data[match(common_id, as.character(dx_data$participant_id)), , drop = FALSE]
  xr_matrix <- xr_matrix[common_id, , drop = FALSE]
  geom_u <- geom_u[common_id, , drop = FALSE]

  # Diagnosis vector coded as binary 0/1.
  if (!target_var %in% names(dx_data)) stop("Variable not found in DX.")
  y <- as.numeric(dx_data[[target_var]])
  
  # Keep rows available in the diagnosis vector, residual matrix, and map.
  ok <- is.finite(y) & complete.cases(xr_matrix)
  y  <- y[ok]
  XR <- xr_matrix[ok, , drop=FALSE]
  U  <- geom_u[ok, , drop=FALSE]
  
  if (sum(y > 0) < 5) {
    message("Not enough positive cases to probe.")
    return(NULL)
  }
  
  cor_u1 <- cor(y, U[,1], method = "spearman")
  cor_u2 <- cor(y, U[,2], method = "spearman")
  
  message(sprintf("\n--- PROBE: %s ---", target_var))
  message(sprintf("Correlation with map U1: %+.3f", cor_u1))
  message(sprintf("Correlation with map U2: %+.3f", cor_u2))
  
  # Correlate each residual item column with the diagnosis.
  cors <- cor(XR, y, method = "spearman", use = "pairwise.complete.obs")
  
  res <- data.frame(
    item = rownames(cors) %||% colnames(XR),
    cor_residual = as.numeric(cors)
  )
  
  # Sort by absolute correlation strength
  res <- res[order(-abs(res$cor_residual)), ]
  
  message(sprintf("\nTop %d residual-item correlations with %s:", top_n, target_var))
  print(head(res, top_n))
  
  invisible(res)
}

# ==============================================================================
# Run the probe
# ==============================================================================

B <- Base_A
pid <- rownames(B)
base_dt <- data.table::data.table(
  participant_id = as.character(pid),
  b1 = as.numeric(B[, 1]), b2 = as.numeric(B[, 2])
)

DT <- data.table::fread(
  "data/wide_diagnoses.csv",
  sep = ";",
  dec = ",",
  na.strings = c("", "NA", "N/A", "NaN", "nan", "null", "NULL", ".", "-"),
  strip.white = TRUE
)
DT <- as.data.table(lapply(DT, char_to_num_local))
stopifnot("participant_id" %in% names(DT))
DT[, participant_id := as.character(participant_id)]

DX <- merge(base_dt, DT, by = "participant_id", all = FALSE, sort = FALSE)

# 1. Probe schizophrenia.
probe_residual_signal("SCID.DIAG.Schizophrenia", DX, XR, geom$U)

# 2. Compare with alcohol abuse.
# probe_residual_signal("SCID.DIAG.AlcAbuse", DX, XR, geom$U)
# 
# # 3. Compare with no diagnosis (NODIAG).
# probe_residual_signal("NODIAG", DX, XR, geom$U)
