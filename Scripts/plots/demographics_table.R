library(tidyverse)
library(readr)
library(flextable)
library(tibble)
library(officer)

cohort_name <- sub("_2$", "", basename(normalizePath(getwd())))
if (!cohort_name %in% c("LA5c", "TCP")) {
  stop("Run this script from the LA5c or TCP project root.")
}

demo_config <- list(
  cohort = cohort_name,
  demo_path = "data/demo.csv",
  dx_path = "data/wide_diagnoses.csv",
  output_path = "out/Table1_Cohort_Characteristics.docx",

  id_col = "participant_id",
  age_col = "age",
  sex_col = "gender",
  nodiag_col = "NODIAG",
  group_col = "group"
)

clean_id <- function(x) trimws(as.character(x))

normalise_binary <- function(x) {
  if (is.logical(x)) return(if_else(is.na(x), NA_integer_, as.integer(x)))
  if (is.numeric(x)) return(if_else(is.na(x), NA_integer_, as.integer(x != 0)))
  
  x <- str_to_lower(clean_id(x))
  case_when(
    x %in% c("1", "yes", "y", "true", "present", "case", "current") ~ 1L,
    x %in% c("0", "no", "n", "false", "absent", "control", "none") ~ 0L,
    TRUE ~ NA_integer_
  )
}

standardise_sex <- function(x) {
  x <- str_to_lower(clean_id(x))
  case_when(
    x %in% c("1", "m", "male") ~ "Male",
    x %in% c("2", "f", "female") ~ "Female",
    x %in% c("o", "other") ~ "Other",
    TRUE ~ NA_character_
  )
}

derive_group_from_id <- function(ids) {
  factor(
    if_else(str_detect(clean_id(ids), "^sub-1"), "General population", "Clinical group"),
    levels = c("General population", "Clinical group")
  )
}

standardise_group <- function(data, config) {
  if (!config$group_col %in% names(data)) {
    return(derive_group_from_id(data[[config$id_col]]))
  }

  x <- str_to_lower(clean_id(data[[config$group_col]]))
  factor(
    case_when(
      x %in% c("genpop", "general population") ~ "General population",
      x %in% c("patient", "clinical group") ~ "Clinical group",
      TRUE ~ NA_character_
    ),
    levels = c("General population", "Clinical group")
  )
}

fmt_mean_sd <- function(x, digits = 1) {
  x <- x[!is.na(x)]
  if (!length(x)) return("")
  sprintf(paste0("%.", digits, "f (%.", digits, "f)"), mean(x), sd(x))
}

fmt_npct <- function(x, digits = 0) {
  x <- x[!is.na(x)]
  d <- length(x)
  if (!d) return("")
  n <- sum(x)
  pct <- 100 * n / d
  paste0(n, " (", formatC(pct, format = "f", digits = digits), "%)")
}

read_delim_clean <- function(path, id_col) {
  read_csv2(path, show_col_types = FALSE) %>%
    mutate("{id_col}" := clean_id(.data[[id_col]]))
}

derive_clinical_status <- function(data, demo_config) {
  if (!demo_config$nodiag_col %in% names(data)) {
    stop("Missing no-diagnosis column: ", demo_config$nodiag_col)
  }

  nodiag <- normalise_binary(data[[demo_config$nodiag_col]])
  factor(
    case_when(
      nodiag == 1L ~ "No diagnosis",
      nodiag == 0L ~ "Any diagnosis",
      TRUE ~ NA_character_
    ),
    levels = c("No diagnosis", "Any diagnosis")
  )
}

demo_raw <- read_delim_clean(demo_config$demo_path, demo_config$id_col)
dx_raw   <- read_delim_clean(demo_config$dx_path, demo_config$id_col)

u_mat <- as.data.frame(geom$U)

u_ids_from_rownames <- rownames(u_mat)
u_ids_from_rownames <- if (!is.null(u_ids_from_rownames)) clean_id(u_ids_from_rownames) else character(0)

demo_ids <- clean_id(demo_raw[[demo_config$id_col]])
overlap_rownames <- sum(u_ids_from_rownames %in% demo_ids)

if (overlap_rownames > 0) {
  u_df <- u_mat %>%
    rownames_to_column(demo_config$id_col) %>%
    mutate("{demo_config$id_col}" := clean_id(.data[[demo_config$id_col]]))
} else if (exists("ids_base", inherits = TRUE) && length(ids_base) == nrow(u_mat)) {
  u_df <- u_mat %>%
    mutate("{demo_config$id_col}" := clean_id(ids_base)) %>%
    relocate(all_of(demo_config$id_col))
} else {
  stop(
    "Could not match geometry IDs to demo.csv IDs. ",
    "Overlap from rownames(geom$U): ", overlap_rownames, ". ",
    "First few demo IDs: ", paste(head(demo_ids, 5), collapse = ", "), "."
  )
}

df <- u_df %>%
  inner_join(demo_raw, by = demo_config$id_col) %>%
  inner_join(dx_raw, by = demo_config$id_col) %>%
  mutate(
    age_std = suppressWarnings(as.numeric(.data[[demo_config$age_col]])),
    sex_std = factor(
      standardise_sex(.data[[demo_config$sex_col]]),
      levels = c("Male", "Female", "Other")
    ),
    group_std = standardise_group(pick(everything()), demo_config),
    clinical_status = derive_clinical_status(pick(everything()), demo_config)
  )

if (nrow(df) == 0) {
  stop("Join produced zero rows. Check geometry IDs and file paths.")
}

table_df <- bind_rows(
  tibble(
    Characteristic = "Sample characteristics",
    Value = "",
    indent = 0L,
    section = TRUE
  ),
  tibble(
    Characteristic = "Age, years",
    Value = fmt_mean_sd(df$age_std, digits = 1),
    indent = 0L,
    section = FALSE
  ),
  tibble(
    Characteristic = "Sex",
    Value = "",
    indent = 0L,
    section = TRUE
  ),
  tibble(
    Characteristic = "Male",
    Value = fmt_npct(df$sex_std == "Male", digits = 0),
    indent = 1L,
    section = FALSE
  ),
  tibble(
    Characteristic = "Female",
    Value = fmt_npct(df$sex_std == "Female", digits = 0),
    indent = 1L,
    section = FALSE
  ),
  if (any(df$sex_std == "Other", na.rm = TRUE)) {
    tibble(
      Characteristic = "Other",
      Value = fmt_npct(df$sex_std == "Other", digits = 0),
      indent = 1L,
      section = FALSE
    )
  },
  tibble(
    Characteristic = "Recruitment group",
    Value = "",
    indent = 0L,
    section = TRUE
  ),
  tibble(
    Characteristic = "General population",
    Value = fmt_npct(df$group_std == "General population", digits = 0),
    indent = 1L,
    section = FALSE
  ),
  tibble(
    Characteristic = "Clinical group",
    Value = fmt_npct(df$group_std == "Clinical group", digits = 0),
    indent = 1L,
    section = FALSE
  ),
  tibble(
    Characteristic = "Clinical status",
    Value = "",
    indent = 0L,
    section = TRUE
  ),
  tibble(
    Characteristic = "No diagnosis",
    Value = fmt_npct(df$clinical_status == "No diagnosis", digits = 0),
    indent = 1L,
    section = FALSE
  ),
  tibble(
    Characteristic = "Any diagnosis",
    Value = fmt_npct(df$clinical_status == "Any diagnosis", digits = 0),
    indent = 1L,
    section = FALSE
  )
)

ft <- flextable(table_df[, c("Characteristic", "Value")]) %>%
  set_header_labels(
    Characteristic = "Characteristic",
    Value = paste0(demo_config$cohort, " (n = ", nrow(df), ")")
  ) %>%
  bold(part = "header", bold = TRUE) %>%
  bold(i = table_df$section, j = 1:2, bold = TRUE, part = "body") %>%
  padding(
    i = !table_df$section & table_df$indent == 1L,
    j = 1,
    padding.left = 18,
    part = "body"
  ) %>%
  align(j = 1, align = "left", part = "all") %>%
  align(j = 2, align = "center", part = "all") %>%
  border_remove() %>%
  hline_top(part = "header", border = fp_border(width = 1.2)) %>%
  hline(part = "header", border = fp_border(width = 1.0)) %>%
  hline(i = nrow(table_df), border = fp_border(width = 1.2), part = "body") %>%
  autofit()

doc <- read_docx() %>%
  body_add_flextable(ft) %>%
  body_add_par("Note. Values are mean (SD) or n (%).", style = "Normal")

print(doc, target = demo_config$output_path)
