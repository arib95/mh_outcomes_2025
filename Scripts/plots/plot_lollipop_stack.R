OUTER_K <- 5L
INNER_K <- 4L
MIN_CLASS <- 10L
BOOT_B <- 1000L
SEED <- 42L
WEIGHT_GRID <- seq(0, 1, by = 0.05)

clamp_probability <- function(x) pmin(pmax(as.numeric(x), 1e-6), 1 - 1e-6)

stratified_folds <- function(y, k, seed) {
  k <- min(k, sum(y == 0L), sum(y == 1L))
  set.seed(seed)
  fold <- integer(length(y))
  for (value in 0:1) {
    index <- sample(which(y == value))
    fold[index] <- rep(seq_len(k), length.out = length(index))
  }
  fold
}

auc <- function(y, score) {
  n1 <- sum(y == 1L)
  n0 <- sum(y == 0L)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

fit_map <- function(y, coordinates) {
  data <- data.frame(y = y, b1 = coordinates[, 1], b2 = coordinates[, 2])
  mgcv::gam(
    y ~ s(b1, bs = "ts", k = 4, m = 2) + s(b2, bs = "ts", k = 4, m = 2),
    data = data,
    family = stats::binomial(),
    method = "REML",
    select = TRUE,
    gamma = 1.8
  )
}

predict_map <- function(model, coordinates) {
  clamp_probability(stats::predict(
    model,
    newdata = data.frame(b1 = coordinates[, 1], b2 = coordinates[, 2]),
    type = "response"
  ))
}

fit_residual <- function(y, residuals, fold) {
  keep <- apply(residuals, 2, function(x) all(is.finite(x)) && stats::sd(x) > 0)
  model <- suppressWarnings(glmnet::cv.glmnet(
    residuals[, keep, drop = FALSE],
    y,
    family = "binomial",
    alpha = 0,
    foldid = fold,
    type.measure = "deviance",
    standardize = TRUE,
    keep = TRUE,
    parallel = FALSE,
    nlambda = 60
  ))
  lambda_index <- which.min(abs(model$lambda - model$lambda.min))
  list(
    model = model,
    keep = keep,
    inner = clamp_probability(model$fit.preval[, lambda_index])
  )
}

predict_residual <- function(model, residuals) {
  clamp_probability(stats::predict(
    model$model,
    newx = residuals[, model$keep, drop = FALSE],
    s = "lambda.min",
    type = "response"
  ))
}

select_weight <- function(y, map_probability, residual_probability) {
  scores <- vapply(WEIGHT_GRID, function(weight) {
    auc(y, (1 - weight) * map_probability + weight * residual_probability)
  }, numeric(1))
  WEIGHT_GRID[which.max(scores)]
}

cross_fitted_predictions <- function(y, coordinates, residuals, seed) {
  outer_fold <- stratified_folds(y, OUTER_K, seed)
  p_map <- p_residual <- p_combined <- rep(NA_real_, length(y))
  selected_weight <- rep(NA_real_, OUTER_K)

  for (outer in seq_len(OUTER_K)) {
    train <- outer_fold != outer
    test <- !train
    inner_fold <- stratified_folds(y[train], INNER_K, seed + 100L * outer)

    residual_model <- fit_residual(
      y[train],
      residuals[train, , drop = FALSE],
      inner_fold
    )
    inner_residual <- residual_model$inner
    inner_map <- rep(NA_real_, sum(train))

    train_coordinates <- coordinates[train, , drop = FALSE]
    for (inner in seq_len(INNER_K)) {
      inner_train <- inner_fold != inner
      inner_test <- !inner_train
      map_model <- fit_map(y[train][inner_train], train_coordinates[inner_train, , drop = FALSE])
      inner_map[inner_test] <- predict_map(map_model, train_coordinates[inner_test, , drop = FALSE])
    }

    weight <- select_weight(y[train], inner_map, inner_residual)
    selected_weight[outer] <- weight

    map_model <- fit_map(y[train], train_coordinates)
    p_map[test] <- predict_map(map_model, coordinates[test, , drop = FALSE])
    p_residual[test] <- predict_residual(residual_model, residuals[test, , drop = FALSE])
    p_combined[test] <- (1 - weight) * p_map[test] + weight * p_residual[test]
  }

  list(
    map = p_map,
    residual = p_residual,
    combined = p_combined,
    fold = outer_fold,
    weight = selected_weight
  )
}

bootstrap_summary <- function(y, predictions, seed) {
  point <- vapply(predictions, function(x) auc(y, x), numeric(1))
  cases <- which(y == 1L)
  controls <- which(y == 0L)
  set.seed(seed)
  draws <- replicate(BOOT_B, {
    index <- c(
      sample(cases, length(cases), replace = TRUE),
      sample(controls, length(controls), replace = TRUE)
    )
    vapply(predictions, function(x) auc(y[index], x[index]), numeric(1))
  })
  delta <- draws["combined", ] - draws["map", ]
  lower_tail <- (1 + sum(delta <= 0)) / (BOOT_B + 1)
  upper_tail <- (1 + sum(delta >= 0)) / (BOOT_B + 1)
  list(
    point = point,
    lo = apply(draws, 1, stats::quantile, probs = 0.025),
    hi = apply(draws, 1, stats::quantile, probs = 0.975),
    delta = point["combined"] - point["map"],
    delta_lo = stats::quantile(delta, 0.025),
    delta_hi = stats::quantile(delta, 0.975),
    p = min(1, 2 * min(lower_tail, upper_tail))
  )
}

state <- readRDS(file.path(OUTPUTS_DIR, "contplot_state.rds"))
residual_state <- readRDS(file.path(OUTPUTS_DIR, "Fprime_matrix.rds"))
diagnoses <- data.table::fread("data/wide_diagnoses.csv", sep = ";", dec = ",")

participant_id <- as.character(residual_state$participant_id)
diagnosis_index <- match(participant_id, as.character(diagnoses$participant_id))
if (anyNA(diagnosis_index)) stop("Participant IDs do not align with the diagnosis table.")
diagnoses <- as.data.frame(diagnoses)[diagnosis_index, , drop = FALSE]
coordinates <- as.matrix(state$Base_A[, c("b1", "b2"), drop = FALSE])
if (!is.null(rownames(coordinates))) {
  coordinate_index <- match(participant_id, rownames(coordinates))
  if (anyNA(coordinate_index)) stop("Participant IDs do not align with the coordinates.")
  coordinates <- coordinates[coordinate_index, , drop = FALSE]
}
residuals <- as.matrix(residual_state$XR)
if (nrow(diagnoses) != nrow(coordinates) || nrow(diagnoses) != nrow(residuals)) {
  stop("Saved artefacts have different participant counts.")
}

outcomes <- names(diagnoses)[grepl("SCID|NODIAG", names(diagnoses))]
outcomes <- outcomes[vapply(outcomes, function(outcome) {
  y <- diagnoses[[outcome]]
  all(stats::na.omit(unique(y)) %in% c(0, 1)) &&
    sum(y == 1, na.rm = TRUE) >= MIN_CLASS &&
    sum(y == 0, na.rm = TRUE) >= MIN_CLASS
}, logical(1))]

metric_rows <- vector("list", length(outcomes))
prediction_rows <- vector("list", length(outcomes))

for (index in seq_along(outcomes)) {
  outcome <- outcomes[index]
  y <- as.integer(diagnoses[[outcome]])
  keep <- is.finite(y) & is.finite(coordinates[, 1]) & is.finite(coordinates[, 2])
  y <- y[keep]
  message(sprintf("[%d/%d] %s", index, length(outcomes), outcome))

  predictions <- cross_fitted_predictions(
    y,
    coordinates[keep, , drop = FALSE],
    residuals[keep, , drop = FALSE],
    SEED + 1000L * index
  )
  summary <- bootstrap_summary(
    y,
    predictions[c("map", "residual", "combined")],
    SEED + 2000L * index
  )

  metric_rows[[index]] <- data.table::data.table(
    var = outcome,
    n = length(y),
    n_cases = sum(y == 1L),
    auc_map = unname(summary$point["map"]),
    auc_map_lo = unname(summary$lo["map"]),
    auc_map_hi = unname(summary$hi["map"]),
    auc_residual = unname(summary$point["residual"]),
    auc_residual_lo = unname(summary$lo["residual"]),
    auc_residual_hi = unname(summary$hi["residual"]),
    auc_combined = unname(summary$point["combined"]),
    auc_combined_lo = unname(summary$lo["combined"]),
    auc_combined_hi = unname(summary$hi["combined"]),
    delta_auc = unname(summary$delta),
    delta_auc_lo = unname(summary$delta_lo),
    delta_auc_hi = unname(summary$delta_hi),
    p_delta = summary$p,
    mean_residual_weight = mean(predictions$weight),
    weight_min = min(predictions$weight),
    weight_max = max(predictions$weight),
    weight_range = diff(range(predictions$weight)),
    weight_unstable = any(predictions$weight <= 0.25) && any(predictions$weight >= 0.75),
    residual_weights = paste(sprintf("%.2f", predictions$weight), collapse = ",")
  )
  prediction_rows[[index]] <- data.table::data.table(
    participant_id = participant_id[keep],
    outcome = outcome,
    observed = y,
    fold = predictions$fold,
    p_map = predictions$map,
    p_residual = predictions$residual,
    p_combined = predictions$combined
  )
}

dt <- data.table::rbindlist(metric_rows)
data.table::fwrite(
  dt,
  file.path(OUTPUTS_DIR, "wide_diagnoses", "residual_contribution_metrics.csv"),
  sep = ";",
  dec = ","
)
saveRDS(
  data.table::rbindlist(prediction_rows),
  file.path(OUTPUTS_DIR, "wide_diagnoses", "residual_contribution_predictions.rds")
)

dt[, label := gsub("^SCID[._]DIAG[._]?", "", var)]
dt[, label := gsub("_", " ", label)]
dt[weight_unstable == TRUE, label := paste0(label, "†")]
dt <- dt[order(delta_auc)]
dt$label <- factor(dt$label, levels = dt$label)
dt[, significance := data.table::fcase(
  p_delta < 0.001, "***",
  p_delta < 0.01, "**",
  p_delta < 0.05, "*",
  default = ""
)]
dt[, colour := data.table::fifelse(p_delta < 0.05 & delta_auc > 0, "Improved", "Other")]

format_ci <- function(point, lo, hi) sprintf("%.3f\n(%.3f–%.3f)", point, lo, hi)
dt[, cases := sprintf("n=%d", n_cases)]
dt[, map_text := format_ci(auc_map, auc_map_lo, auc_map_hi)]
dt[, combined_text := format_ci(auc_combined, auc_combined_lo, auc_combined_hi)]
dt[, p_text := sprintf("%.3f", p_delta)]
dt[p_delta < 0.001, p_text := "<.001"]

p_table <- ggplot(dt, aes(y = label)) +
  geom_text(aes(x = 0, label = label), hjust = 0, fontface = "bold", size = 3.2) +
  geom_text(aes(x = 2.7, label = cases), hjust = 0, size = 3, color = "grey40") +
  geom_text(aes(x = 3.7, label = map_text), hjust = 0, size = 2.8, lineheight = 0.9) +
  geom_text(aes(x = 5.1, label = combined_text), hjust = 0, size = 2.8, lineheight = 0.9, fontface = "bold") +
  geom_text(aes(x = 6.5, label = p_text), hjust = 0, size = 2.8, fontface = "italic", color = "grey30") +
  scale_x_continuous(
    limits = c(0, 6.9),
    expand = c(0, 0),
    position = "top",
    breaks = c(0, 2.7, 3.7, 5.1, 6.5),
    labels = c("Diagnosis", "Cases", "Map AUC\n(95% CI)", "Map+residual\nAUC (95% CI)", "p")
  ) +
  theme_void() +
  theme(
    axis.text.x.top = element_text(face = "bold", size = 8, hjust = 0, vjust = 0),
    plot.margin = margin(r = 0, l = 5)
  )

p_forest <- ggplot(dt, aes(x = delta_auc, y = label, color = colour)) +
  geom_vline(xintercept = 0, color = "grey80") +
  geom_segment(aes(x = 0, xend = delta_auc, yend = label), linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(aes(label = significance), vjust = -0.5, size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("Improved" = "#2E86C1", "Other" = "grey70")) +
  scale_x_continuous(
    labels = function(x) sprintf("%+.2f", x),
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(x = "AU-ROC gain from residual terms", color = NULL) +
  theme_classic() +
  theme(
    axis.line.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.grid.major.x = element_line(color = "grey95"),
    plot.margin = margin(l = 0)
  )

combined_plot <- p_table + p_forest +
  plot_layout(widths = c(1.7, 1)) +
  plot_annotation(
    caption = "† Residual weight varied from ≤0.25 to ≥0.75 across outer folds.",
    theme = theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.caption = element_text(size = 8, color = "grey35", hjust = 0)
    )
  )

save_plot_gg(
  "Fig_Resid_Contribution",
  combined_plot,
  width = 9.5,
  height = 0.25 * nrow(dt) + 2,
  save_rds = TRUE
)

print(combined_plot)
