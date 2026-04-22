library(tidyverse)
library(readr)
library(flextable)
library(tibble)
library(officer)

# cfg <- list(
#   cohort = "LA5c",
#   demo_path = "data/demo.csv",
#   dx_path = "data/wide_diagnoses.csv",
#   output_path = "out/Table1_Cohort_Characteristics.docx",
#
#   id_col = "participant_id",
#   age_col = "age",
#
#   sex_col = "gender",
#   sex_map = c("1" = "Male", "2" = "Female"),
#
#   nodiag_col = "NODIAG",
#
#   group_source = "id"
# )

cfg <- list(
  cohort = "TCP",
  demo_path = "data/demo.csv",
  dx_path = "data/wide_diagnoses.csv",
  output_path = "out/Table1_Cohort_Characteristics.docx",
  
  id_col = "participant_id",
  age_col = "age",
  
  sex_col = "sex",
  sex_map = c("M" = "Male", "F" = "Female", "Male" = "Male", "Female" = "Female"),
  
  nodiag_col = "NODIAG",
  ndiagnoses_col = "NDIAGNOSES",
  
  group_source = "column",
  group_col = "group",
  group_map = c(
    "GenPop" = "General population",
    "Patient" = "Clinical group"
  )
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

standardise_from_map <- function(x, map) {
  x_chr <- clean_id(x)
  out <- unname(map[x_chr])
  if_else(is.na(out), x_chr, out)
}

derive_group_from_id <- function(ids) {
  factor(
    if_else(str_detect(clean_id(ids), "^sub-1"), "General population", "Clinical group"),
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
  read_delim(path, delim = ";", show_col_types = FALSE) %>%
    mutate("{id_col}" := clean_id(.data[[id_col]]))
}

derive_clinical_status <- function(data, cfg) {
  if (cfg$cohort == "TCP") {
    if (!cfg$ndiagnoses_col %in% names(data)) {
      stop("Missing TCP diagnoses-count column: ", cfg$ndiagnoses_col)
    }
    
    ndiag <- suppressWarnings(as.numeric(data[[cfg$ndiagnoses_col]]))
    
    factor(
      case_when(
        is.na(ndiag) ~ NA_character_,
        ndiag == 0 ~ "No diagnosis",
        ndiag > 0 ~ "Any diagnosis",
        TRUE ~ NA_character_
      ),
      levels = c("No diagnosis", "Any diagnosis")
    )
  } else {
    if (!cfg$nodiag_col %in% names(data)) {
      stop("Missing LA5c no-diagnosis column: ", cfg$nodiag_col)
    }
    
    factor(
      case_when(
        normalise_binary(data[[cfg$nodiag_col]]) == 1L ~ "No diagnosis",
        normalise_binary(data[[cfg$nodiag_col]]) == 0L ~ "Any diagnosis",
        TRUE ~ NA_character_
      ),
      levels = c("No diagnosis", "Any diagnosis")
    )
  }
}

demo_raw <- read_delim_clean(cfg$demo_path, cfg$id_col)
dx_raw   <- read_delim_clean(cfg$dx_path, cfg$id_col)

u_mat <- as.data.frame(geom$U)

u_ids_from_rownames <- rownames(u_mat)
u_ids_from_rownames <- if (!is.null(u_ids_from_rownames)) clean_id(u_ids_from_rownames) else character(0)

demo_ids <- clean_id(demo_raw[[cfg$id_col]])
overlap_rownames <- sum(u_ids_from_rownames %in% demo_ids)

if (overlap_rownames > 0) {
  u_df <- u_mat %>%
    rownames_to_column(cfg$id_col) %>%
    mutate("{cfg$id_col}" := clean_id(.data[[cfg$id_col]]))
} else if (exists("ids_base", inherits = TRUE) && length(ids_base) == nrow(u_mat)) {
  u_df <- u_mat %>%
    mutate("{cfg$id_col}" := clean_id(ids_base)) %>%
    relocate(all_of(cfg$id_col))
} else {
  stop(
    "Could not match geometry IDs to demo.csv IDs. ",
    "Overlap from rownames(geom$U): ", overlap_rownames, ". ",
    "First few demo IDs: ", paste(head(demo_ids, 5), collapse = ", "), "."
  )
}

df <- u_df %>%
  inner_join(demo_raw, by = cfg$id_col) %>%
  inner_join(dx_raw, by = cfg$id_col) %>%
  mutate(
    age_std = suppressWarnings(as.numeric(.data[[cfg$age_col]])),
    sex_std = factor(
      standardise_from_map(.data[[cfg$sex_col]], cfg$sex_map),
      levels = c("Male", "Female")
    ),
    group_std = if (cfg$group_source == "column") {
      factor(
        standardise_from_map(.data[[cfg$group_col]], cfg$group_map),
        levels = c("General population", "Clinical group")
      )
    } else {
      derive_group_from_id(.data[[cfg$id_col]])
    },
    clinical_status = derive_clinical_status(cur_data(), cfg)
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
    Value = paste0(cfg$cohort, " (n = ", nrow(df), ")")
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

print(doc, target = cfg$output_path)