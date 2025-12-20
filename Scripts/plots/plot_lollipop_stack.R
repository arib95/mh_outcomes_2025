# ==================================================
# 1. Load inputs
# ==================================================

csv_path <- "wide_diagnoses/behaviour_map_metrics.csv"
full_csv_path <- file.path(getwd(), OUTPUTS_DIR, csv_path)

if (!file.exists(full_csv_path)) {
  stop(paste("Could not find:", full_csv_path))
}

MET <- data.table::fread(full_csv_path)

dt <- MET[var %like% "SCID|NODIAG"]

# ==================================================
# 2. Add case counts from raw diagnosis file
# ==================================================

diag_path <- "data/wide_diagnoses.csv"

if (file.exists(diag_path)) {
  raw_df <- data.table::fread(diag_path)
  
  get_cases <- function(v) {
    if (v == "NODIAG") {
      if ("NODIAG" %in% names(raw_df)) return(sum(raw_df$NODIAG == 1, na.rm = TRUE))
      return(NA)
    }
    if (v %in% names(raw_df)) return(sum(raw_df[[v]] == 1, na.rm = TRUE))
    NA
  }
  
  dt[, n_cases := sapply(var, get_cases)]
} else {
  warning("Raw diagnosis file not found; n_cases will be NA.")
  dt[, n_cases := NA]
}

# ==================================================
# 3. Derivations and display fields
# ==================================================

# Added value of residuals
dt[, delta_auc := auc_stacked - oof_point]

# Display label
dt[, label := gsub("SCID.DIAG.", "", var)]

# Order (largest negative to largest positive; consistent with current code)
dt <- dt[order(delta_auc)]
dt$label <- factor(dt$label, levels = dt$label)

# Significance annotations
dt[, sig_cat := fcase(
  p_delta_stacked < 0.001, "***",
  p_delta_stacked < 0.01,  "**",
  p_delta_stacked < 0.05,  "*",
  default = ""
)]
dt[, color_grp := fifelse(p_delta_stacked < 0.05 & delta_auc > 0, "Sig", "Ns")]

# Text columns for the left “table” panel
fmt_ci <- function(pt, lo, hi) sprintf("%.3f\n(%.3f–%.3f)", pt, lo, hi)

dt[, txt_cases := sprintf("n=%d", n_cases)]
dt[, txt_base  := fmt_ci(oof_point, oof_lo, oof_hi)]

if ("auc_stacked_lo" %in% names(dt)) {
  dt[, txt_stack := fmt_ci(auc_stacked, auc_stacked_lo, auc_stacked_hi)]
} else {
  dt[, txt_stack := sprintf("%.3f", auc_stacked)]
}

dt[, txt_pval := sprintf("%.3f", p_delta_stacked)]
dt[p_delta_stacked < 0.001, txt_pval := "<.001"]

# ==================================================
# 4. Left panel: table-style annotation plot
# ==================================================

p_table <- ggplot(dt, aes(y = label)) +
  geom_text(aes(x = 0,   label = label),     hjust = 0, fontface = "bold", size = 3.2) +
  geom_text(aes(x = 2.2, label = txt_cases), hjust = 0, size = 3,   color = "grey40") +
  geom_text(aes(x = 3.2, label = txt_base),  hjust = 0, size = 2.8, lineheight = 0.9) +
  geom_text(aes(x = 4.5, label = txt_stack), hjust = 0, size = 2.8, lineheight = 0.9, fontface = "bold") +
  geom_text(aes(x = 5.8, label = txt_pval),  hjust = 0, size = 2.8, fontface = "italic", color = "grey30") +
  scale_x_continuous(
    limits = c(0, 6.2),
    expand = c(0, 0),
    position = "top",
    breaks = c(0, 2.2, 3.2, 4.5, 5.8),
    labels = c("Diagnosis", "Cases", "Base AUC\n(95% CI)", "Full AUC\n(Base+Resid)", "p")
  ) +
  theme_void() +
  theme(
    axis.text.x.top = element_text(face = "bold", size = 8, hjust = 0, vjust = 0),
    plot.margin = margin(r = 0, l = 5)
  )

# ==================================================
# 5. Right panel: delta “lollipop” plot
# ==================================================

p_forest <- ggplot(dt, aes(x = delta_auc, y = label, color = color_grp)) +
  geom_vline(xintercept = 0, linetype = "solid", color = "grey80") +
  geom_segment(aes(x = 0, xend = delta_auc, y = label, yend = label), linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(aes(label = sig_cat), vjust = -0.5, size = 3.5, show.legend = FALSE) +
  scale_color_manual(values = c("Sig" = "#2E86C1", "Ns" = "grey70")) +
  scale_x_continuous(
    labels = function(x) sprintf("+%.2f", x),
    expand = expansion(mult = c(0, 0.1))
  ) +
  labs(x = "Δ AU-ROC (Added Value of Residuals)", color = NULL) +
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

# ==================================================
# 6. Combine and save
# ==================================================

combined_plot <- p_table + p_forest +
  plot_layout(widths = c(1.4, 1)) +
  plot_annotation(theme = theme(plot.title = element_text(face = "bold", size = 14)))

save_plot_gg(
  "Fig_Resid_Contribution",
  combined_plot,
  width = 8,
  height = 0.25 * nrow(dt) + 2,
  save_rds = TRUE
)

print(combined_plot)