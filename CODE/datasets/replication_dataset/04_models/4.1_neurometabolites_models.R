
library(dplyr)
library(readr)
library(tidyr)
library(lme4)
library(tibble)


# ============================================================
# Configuration
# ============================================================

field <- "Hursts"

task_levels <- c("Rest", "0Back", "1Back", "2Back")
metabolites <- c("EIB", "GABA", "Glx")

selected_parcels <- c(26, 27, 35, 43, 45)

base_path <- normalizePath(getwd(), winslash = "/")

MRS_path <- file.path(
  base_path,
  "DATA/replication_dataset/EIB_12s_MRS_static_&_AUC_New_criteria.csv"
)

Hurst_path <- file.path(
  base_path,
  paste0(
    "RESULTS/replication_dataset/hurst/dataframes/",
    "df_hurst_parcelwise_FIN_D4_modwt_1_5_padding_0_",
    "decimate_0_residuals_w_rest_IQR_cleaned.csv"
  )
)

output_dir <- file.path(
  base_path,
  "RESULTS/replication_dataset/models"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Load data
# ============================================================

mrs_df <- read_csv(
  MRS_path,
  show_col_types = FALSE
) %>%
  rename(Task = Runs) %>%
  mutate(
    Task = recode(
      Task,
      "ARest" = "Rest",
      "B0Back" = "0Back",
      "C1Back" = "1Back",
      "D2Back" = "2Back"
    ),
    Task = factor(Task, levels = task_levels)
  )

hurst_df <- read_csv(
  Hurst_path,
  show_col_types = FALSE
) %>%
  select(-starts_with("Unnamed"))


# ============================================================
# Merge Hurst and MRS data
# ============================================================

merged_df <- hurst_df %>%
  filter(Parcel %in% selected_parcels) %>%
  left_join(
    mrs_df %>%
      select(Subject, Task, EIB, GABA, Glx),
    by = c("Subject", "Task")
  ) %>%
  mutate(
    Task = factor(Task, levels = task_levels),
    Subject = factor(Subject),
    Parcel = factor(Parcel)
  )


# ============================================================
# Define reference EIB sample
# ============================================================

eib_df <- merged_df %>%
  drop_na(
    all_of(c(field, "EIB", "Task", "Subject", "Parcel"))
  ) %>%
  droplevels()

reference_subjects <- unique(as.character(eib_df$Subject))


# ============================================================
# Fit models and export results for Python
# ============================================================

for (measure in metabolites) {
  
  cat("\nFitting model for:", measure, "\n")
  
  # Prepare data.
  df <- merged_df %>%
    drop_na(
      all_of(c(
        field, measure, "Task", "Subject", "Parcel"
      ))
    )
  
  # Restrict GABA and Glx to subjects in the EIB sample.
  if (measure != "EIB") {
    df <- df %>%
      filter(
        as.character(Subject) %in% reference_subjects
      )
  }
  
  df <- df %>%
    droplevels() %>%
    mutate(
      Task = factor(Task, levels = task_levels),
      Subject = factor(Subject),
      Parcel = factor(Parcel)
    )
  
  # Fit linear mixed-effects model.
  formula <- as.formula(
    paste0(
      field, " ~ ", measure,
      " * Task + (1 | Subject) + (1 | Parcel)"
    )
  )
  
  fit <- lmer(
    formula,
    data = df,
    REML = TRUE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  )
  
  # Output prefix.
  prefix <- paste0(
    "lmm_", measure, "_", field, "_full_task"
  )
  
  # 1. Save model dataframe.
  write_csv(
    df,
    file.path(
      output_dir,
      paste0(prefix, "_df.csv")
    )
  )
  
  # 2. Save fixed-effect coefficients.
  fixef_df <- data.frame(
    Term = names(fixef(fit)),
    Estimate = unname(fixef(fit))
  )
  
  write_csv(
    fixef_df,
    file.path(
      output_dir,
      paste0(prefix, "_fixef.csv")
    )
  )
  
  # 3. Save fixed-effect covariance matrix.
  vcov_df <- as.data.frame(as.matrix(vcov(fit))) %>%
    rownames_to_column("Term")
  
  write_csv(
    vcov_df,
    file.path(
      output_dir,
      paste0(prefix, "_vcov.csv")
    )
  )
  
  cat(
    "Observations:", nrow(df),
    "| Subjects:", n_distinct(df$Subject),
    "| Parcels:", n_distinct(df$Parcel),
    "\n"
  )
}

cat("\nDone. Plotting data saved to:\n", output_dir, "\n")
