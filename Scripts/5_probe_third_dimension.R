# ==============================================================================
# Residual signal probe
# ==============================================================================
# Needs to be run after running dimension_contplot as "wide_diagnoses" so
# its XR is kept in environment.

probe_residual_signal <- function(target_var, dx_data, xr_matrix, geom_u, top_n = 20) {
  
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
    item = rownames(cors),
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
DT <- as.data.table(lapply(DT, char_to_num))
stopifnot("participant_id" %in% names(DT))
DT[, participant_id := as.character(participant_id)]

DX <- merge(base_dt, DT, by = "participant_id", all = FALSE)

# 1. Probe schizophrenia.
probe_residual_signal("SCID.DIAG.Schizophrenia", DX, XR, geom$U)

# 2. Compare with alcohol abuse.
# probe_residual_signal("SCID.DIAG.AlcAbuse", DX, XR, geom$U)
# 
# # 3. Compare with no diagnosis (NODIAG).
# probe_residual_signal("NODIAG", DX, XR, geom$U)
