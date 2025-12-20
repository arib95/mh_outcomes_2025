library(tidyverse)
library(gtsummary)
library(flextable)
library(readr)

# ==================================================
# 0. Base U-space coordinates and quadrant assignment
# ==================================================

base_df <- tibble(
  participant_id = trimws(as.character(ids_base)),
  u1 = Base[, 1],
  u2 = Base[, 2]
) %>%
  mutate(
    Quadrant = case_when(
      u1 > 0  & u2 > 0  ~ "Q1 (Top-Right)",
      u1 <= 0 & u2 > 0  ~ "Q2 (Top-Left)",
      u1 <= 0 & u2 <= 0 ~ "Q3 (Bottom-Left)",
      u1 > 0  & u2 <= 0 ~ "Q4 (Bottom-Right)"
    ),
    Quadrant = factor(
      Quadrant,
      levels = c(
        "Q1 (Top-Right)",
        "Q2 (Top-Left)",
        "Q3 (Bottom-Left)",
        "Q4 (Bottom-Right)"
      )
    )
  )

# ======================================
# 1. Demographics and diagnostic status
# ======================================

demo_raw <- read_delim(
  "data/demo.csv",
  delim = ";",
  show_col_types = FALSE
) %>%
  mutate(participant_id = trimws(as.character(participant_id)))

DX_wide <- read_delim(
  "data/wide_diagnoses.csv",
  delim = ";",
  show_col_types = FALSE
) %>%
  mutate(participant_id = trimws(as.character(participant_id)))

clinical_status <- DX_wide %>%
  select(participant_id, NODIAG)

# ======================================
# 2. Master LA5c analysis dataframe
# ======================================

# Extract U coordinates from geometry object
u_df <- geom$U %>%
  as.data.frame() %>%
  rownames_to_column("participant_id")

la5c_df <- u_df %>%
  mutate(
    Quadrant = case_when(
      u1 > 0  & u2 > 0  ~ "Q1 (Top-Right)",
      u1 <= 0 & u2 > 0  ~ "Q2 (Top-Left)",
      u1 <= 0 & u2 <= 0 ~ "Q3 (Bottom-Left)",
      u1 > 0  & u2 <= 0 ~ "Q4 (Bottom-Right)"
    ),
    Quadrant = factor(
      Quadrant,
      levels = c(
        "Q1 (Top-Right)",
        "Q2 (Top-Left)",
        "Q3 (Bottom-Left)",
        "Q4 (Bottom-Right)"
      )
    )
  ) %>%
  inner_join(demo_raw, by = "participant_id") %>%
  inner_join(clinical_status, by = "participant_id") %>%
  mutate(
    Status = factor(
      if_else(NODIAG == 1, "No diagnosis", "Any lifetime diagnosis"),
      levels = c("Any lifetime diagnosis", "No diagnosis")
    ),
    gender  = factor(gender, levels = c(1, 2), labels = c("Male", "Female")),
    smoking = factor(cigs, levels = c(0, 1), labels = c("Non-smoker", "Smoker")),
    race = case_when(
      race_main == 1 ~ "AI/AN",
      race_main == 2 ~ "Asian",
      race_main == 3 ~ "NH/PI",
      race_main == 4 ~ "Black",
      race_main == 5 ~ "White",
      TRUE ~ "Other"
    ),
    race  = factor(race, levels = c("AI/AN", "Asian", "Black", "Other", "White")),
    group = factor(group)
  )

# Basic sanity checks
nrow(la5c_df)
summary(la5c_df$u1)
summary(la5c_df$u2)

# ======================================
# 3. Table 1: Demographics by quadrant
# ======================================

table_quads <- la5c_df %>%
  select(
    Quadrant, age, gender, race,
    school_yrs, smoking, group, Status
  ) %>%
  tbl_summary(
    by = Quadrant,
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      age        ~ "Age (Years)",
      gender     ~ "Sex",
      race       ~ "Race",
      school_yrs ~ "Education (Years)",
      smoking    ~ "Smoking Status",
      group      ~ "Recruitment Group",
      Status     ~ "Clinical Status (SCID)"
    ),
    missing = "no"
  ) %>%
  add_overall() %>%
  add_p(
    test = list(
      all_continuous()  ~ "kruskal.test",
      all_categorical() ~ "fisher.test"
    ),
    test.args = list(
      race ~ list(workspace = 2e8)
    ),
    pvalue_fun = ~ style_pvalue(.x, digits = 3)
  ) %>%
  bold_labels() %>%
  modify_header(stat_0 ~ "**Overall (N = {N})**") %>%
  modify_footnote(everything() ~ NA_character_)

table_quads %>%
  as_flex_table() %>%
  save_as_docx(path = "out/Table1_Demographics_by_Quadrant.docx")

# ======================================
# 4. U-space quadrant plot
# ======================================

plot_df <- la5c_df %>%
  mutate(
    PlotStatus = factor(
      if_else(
        NODIAG == 1,
        "No SCID Diagnosis",
        "Any Lifetime SCID Diagnosis"
      ),
      levels = c(
        "Any Lifetime SCID Diagnosis",
        "No SCID Diagnosis"
      )
    )
  )

# Plot origin (centroid)
cm_u1 <- 0
cm_u2 <- 0

p_quadrants <- ggplot(
  plot_df,
  aes(x = u1, y = u2, colour = PlotStatus)
) +
  annotate(
    "path",
    x = cm_u1 + cos(seq(0, 2 * pi, length.out = 200)),
    y = cm_u2 + sin(seq(0, 2 * pi, length.out = 200)),
    colour = "grey90",
    linewidth = 1
  ) +
  geom_vline(
    xintercept = cm_u1,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.6
  ) +
  geom_hline(
    yintercept = cm_u2,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.6
  ) +
  geom_point(alpha = 0.65, size = 2.5, stroke = 0) +
  annotate(
    "label",
    x = cm_u1 + c(0.75, -0.75, -0.75, 0.75),
    y = cm_u2 + c(0.75, 0.75, -0.75, -0.75),
    label = c("Q1", "Q2", "Q3", "Q4"),
    fontface = "bold",
    colour = "black",
    alpha = 0.8
  ) +
  scale_color_manual(
    values = c(
      "Any Lifetime SCID Diagnosis" = "#E5D8BD",
      "No SCID Diagnosis" = "black"
    )
  ) +
  coord_equal() +
  labs(x = "u1", y = "u2", colour = "Diagnostic Group") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

save_plot_gg(
  "FIG_participant_quadrant_map",
  p_quadrants,
  width = 8,
  height = 7.5,
  dpi = 300
)