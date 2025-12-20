# scripts/surprisal_vs_scz_logit.R
# Fits: SCID schizophrenia membership ~ surprisal
# Inputs:
#   - data/acquiescence.csv (participant_id, surprise)
#   - data/wide_diagnoses.csv (participant_id, SCID.DIAG.Schizophrenia)
# Outputs:
#   - out/surprisal_scz_logit_curve.csv (prediction grid with 95% CI)
#   - out/surprisal_scz_logit.png (optional plot)

suppressPackageStartupMessages({
  library(data.table)
})

# ==================================================
# Helpers
# ==================================================

# Coerce character-like numeric fields (incl. commas/percent) to numeric when mostly parsable
char_to_num <- function(v) {
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

# ==================================================
# Load and merge inputs
# ==================================================

aq <- fread(
  "data/acquiescence.csv",
  na.strings = c("", "NA", "N/A", "NaN", "nan", "null", "NULL", ".", "-"),
  strip.white = TRUE
)
aq <- as.data.table(lapply(aq, char_to_num))
stopifnot(all(c("participant_id", "surprise") %in% names(aq)))
aq[, participant_id := as.character(participant_id)]

dx <- fread(
  "data/wide_diagnoses.csv",
  na.strings = c("", "NA", "N/A", "NaN", "nan", "null", "NULL", ".", "-"),
  strip.white = TRUE
)
dx <- as.data.table(lapply(dx, char_to_num))
stopifnot("participant_id" %in% names(dx), "SCID.DIAG.BipolarI" %in% names(dx))
dx[, participant_id := as.character(participant_id)]

D <- merge(
  aq[, .(participant_id, surprise = as.numeric(surprise))],
  dx[, .(participant_id, SCID.DIAG.BipolarI)],
  by = "participant_id",
  all = FALSE
)

# ==================================================
# Outcome construction and QC
# ==================================================

raw_dx <- D[["SCID.DIAG.BipolarI"]]

D[, y := {
  if (is.numeric(raw_dx)) {
    as.integer(raw_dx > 0)
  } else {
    as.integer(tolower(trimws(as.character(raw_dx))) %in% c("1", "yes", "true", "y"))
  }
}]

D <- D[is.finite(surprise) & !is.na(y)]

if (nrow(D) < 10L || length(unique(D$y)) < 2L) {
  stop("Insufficient data or no case/control variation after merge.")
}

# ==================================================
# Model fit: logit(y ~ surprise)
# ==================================================

fit <- glm(y ~ surprise, data = D, family = binomial())

# ==================================================
# Prediction grid with 95% CI
# ==================================================

xgrid <- data.frame(
  surprise = seq(
    min(D$surprise, na.rm = TRUE),
    max(D$surprise, na.rm = TRUE),
    length.out = 300
  )
)

pr <- predict(fit, newdata = xgrid, type = "link", se.fit = TRUE)
eta <- pr$fit
se  <- pr$se.fit
inv <- fit$family$linkinv

curve_dt <- data.table(
  surprise = xgrid$surprise,
  p  = inv(eta),
  lo = inv(eta - 1.96 * se),
  hi = inv(eta + 1.96 * se)
)

# ==================================================
# Export artefacts
# ==================================================

dir.create("out", showWarnings = FALSE, recursive = TRUE)
fwrite(curve_dt, "out/surprisal_scz_logit_curve.csv")

# Optional plot
if (requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  
  p <- ggplot() +
    geom_point(
      data = D,
      aes(x = surprise, y = y),
      alpha = 0.25,
      shape = 16,
      position = position_jitter(height = 0.04, width = 0)
    ) +
    geom_ribbon(
      data = curve_dt,
      aes(x = surprise, ymin = lo, ymax = hi),
      alpha = 0.18
    ) +
    geom_line(
      data = curve_dt,
      aes(x = surprise, y = p),
      linewidth = 1
    ) +
    scale_y_continuous("Schizophrenia membership (probability)", limits = c(0, 1)) +
    scale_x_continuous("Surprisal") +
    labs(title = "Logit curve: surprisal → schizophrenia membership") +
    theme_minimal(base_size = 12)
  
  ggsave("out/surprisal_scz_logit.png", p, width = 7, height = 4.6, dpi = 200)
}

# ==================================================
# Console summary
# ==================================================

s  <- summary(fit)$coef["surprise", ]
OR <- exp(unname(s["Estimate"]))
CI <- tryCatch(exp(confint(fit, "surprise")), error = function(e) c(NA, NA))

cat(sprintf(
  "Logit(y ~ surprisal): OR per +1 unit = %.3f (95%% CI %.3f–%.3f), p=%.3g; n=%d\n",
  OR, CI[1], CI[2], unname(s["Pr(>|z|)"]), nrow(D)
))