# ==============================================================================
# PROBE: What lives in the "Third Dimension"?
# ==============================================================================
# Needs to be run after running dimension_contplot as "wide_diagnoses" so
# its XR is kept in environment.

probe_residual_signal <- function(target_var, dx_data, xr_matrix, geom_u, top_n = 20) {
  
  # 1. Get the diagnosis vector (binary 0/1)
  if (!target_var %in% names(dx_data)) stop("Variable not found in DX.")
  y <- as.numeric(dx_data[[target_var]])
  
  # Filter to complete cases (intersection of diagnosis and residuals)
  # Assuming XR rows match DX rows (if not, we'd need ID matching)
  ok <- is.finite(y) & complete.cases(xr_matrix)
  y  <- y[ok]
  XR <- xr_matrix[ok, , drop=FALSE]
  U  <- geom_u[ok, , drop=FALSE]
  
  if (sum(y > 0) < 5) {
    message("Not enough positive cases to probe.")
    return(NULL)
  }
  
  # 2. Check correlation with the 2D Map (Base)
  # This confirms "how much does the map know?"
  cor_u1 <- cor(y, U[,1], method = "spearman")
  cor_u2 <- cor(y, U[,2], method = "spearman")
  
  message(sprintf("\n--- PROBE: %s ---", target_var))
  message(sprintf("Correlation with Map (U1): %+.3f", cor_u1))
  message(sprintf("Correlation with Map (U2): %+.3f", cor_u2))
  
  # 3. Check correlation with Residuals (Fibre)
  # This finds "what did we throw away that predicts this?"
  # We use cor() on the matrix for speed (vectorized)
  cors <- cor(XR, y, method = "spearman", use = "pairwise.complete.obs")
  
  # 4. Format results
  res <- data.frame(
    item = rownames(cors),
    cor_residual = as.numeric(cors)
  )
  
  # Sort by absolute correlation strength
  res <- res[order(-abs(res$cor_residual)), ]
  
  # 5. Display
  message(sprintf("\nTop %d items in the Residuals predicting %s:", top_n, target_var))
  print(head(res, top_n))
  
  invisible(res)
}

# ==============================================================================
# RUN THE PROBE
# ==============================================================================

B <- Base_A
pid <- rownames(B)
base_dt <- data.table::data.table(
  participant_id = as.character(pid),
  b1 = as.numeric(B[, 1]), b2 = as.numeric(B[, 2])
)

DT <- data.table::fread(
  "data/wide_diagnoses.csv",
  na.strings = c("", "NA", "N/A", "NaN", "nan", "null", "NULL", ".", "-"),
  strip.white = TRUE
)
DT <- as.data.table(lapply(DT, char_to_num))
stopifnot("participant_id" %in% names(DT))
DT[, participant_id := as.character(participant_id)]

DX <- merge(base_dt, DT, by = "participant_id", all = FALSE)

# 1. Probe Schizophrenia (The hidden signal)
probe_residual_signal("SCID.DIAG.Schizophrenia", DX, XR, geom$U)

# 2. Compare with Alcohol Abuse (Another high-residual variable)
# probe_residual_signal("SCID.DIAG.AlcAbuse", DX, XR, geom$U)
# 
# # 3. Compare with "Healthy" (NODIAG) - should be driven by the Map, not residuals
# probe_residual_signal("NODIAG", DX, XR, geom$U)