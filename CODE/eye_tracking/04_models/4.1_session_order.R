library(dplyr)
library(readr)
library(lme4)
library(lmerTest)

#### Config ####

base_path <- normalizePath(getwd(), winslash = "/")
data_type <- "microsaccades"
decimation <- "decimated_2s"

data_path <- file.path(
  base_path,
  "RESULTS/pupillometry",
  data_type,
  "hurst/dataframes",
  paste0("df_hurst_", data_type, "_FIN_D4_modwt_1_5_", decimation, ".csv")
)

output_path <- file.path(
  base_path,
  "RESULTS/pupillometry",
  data_type,
  "/models/order/decimated_2s/range_1_5"
)

TASK_COL <- "Task"
MODALITY_COL <- "Modality"
SUBJECT_COL <- "Subject"
VALUE_COL <- "Hurst"
ORDER_COL <- "order"

MODALITY_ORDER <- c("fMRI", "fMRS")

dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

#### Load data ####

df <- read_csv(data_path, show_col_types = FALSE)

required_cols <- c(VALUE_COL, TASK_COL, MODALITY_COL, SUBJECT_COL, ORDER_COL)
missing_cols <- setdiff(required_cols, colnames(df))

if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df <- df %>%
  filter(if_all(all_of(required_cols), ~ !is.na(.))) %>%
  mutate(
    Modality = factor(.data[[MODALITY_COL]]),
    Task = factor(.data[[TASK_COL]]),
    Subject = factor(.data[[SUBJECT_COL]])
  )

task_order <- sort(unique(as.character(df[[TASK_COL]])))

df <- df %>%
  mutate(
    Task = factor(.data[[TASK_COL]], levels = task_order),
    Modality = factor(.data[[MODALITY_COL]], levels = MODALITY_ORDER),
    Subject = factor(.data[[SUBJECT_COL]])
  )

#### Model ####

formula_str <- paste0(
  VALUE_COL, " ~ ",
  TASK_COL, " + ",
  ORDER_COL, " + ",
  TASK_COL, ":", ORDER_COL,
  " + (1 | ", SUBJECT_COL, ")"
)

model_formula <- as.formula(formula_str)

fit <- lmer(
  model_formula,
  data = df,
  REML = FALSE,
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(maxfun = 2e5)
  )
)

print(summary(fit))

#### Fixed-effects table ####

coef_mat <- summary(fit)$coefficients
coef_df <- as.data.frame(coef_mat)
coef_df$term <- rownames(coef_df)
rownames(coef_df) <- NULL

p_col <- grep("Pr\\(", colnames(coef_df), value = TRUE)

if (length(p_col) == 0) {
  stop("No p-value column found. Make sure lmerTest is loaded.")
}

ci_mat <- confint(fit, parm = names(fixef(fit)), method = "Wald")

fixed_effects_df <- data.frame(
  term = coef_df$term,
  coef = coef_df$Estimate,
  std_err = coef_df$`Std. Error`,
  z = coef_df$`t value`,
  p_value = coef_df[[p_col[1]]],
  ci_low = ci_mat[coef_df$term, 1],
  ci_high = ci_mat[coef_df$term, 2],
  row.names = NULL,
  check.names = FALSE
)

#### Prepare data for plotting ####

plot_df <- df %>%
  group_by(
    .data[[SUBJECT_COL]],
    .data[[TASK_COL]],
    .data[[MODALITY_COL]]
  ) %>%
  summarise(
    Hurst = mean(.data[[VALUE_COL]], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Task = factor(.data[[TASK_COL]], levels = task_order, ordered = TRUE),
    Modality = factor(.data[[MODALITY_COL]], levels = MODALITY_ORDER, ordered = TRUE)
  )

#### Save ####

write_csv(
  plot_df,
  file.path(output_path, "plot_df.csv")
)

write_csv(
  fixed_effects_df,
  file.path(output_path, "fixed_effects.csv")
)

capture.output(
  summary(fit),
  file = file.path(output_path, "model_summary.txt")
)

analysis_params <- list(
  task_order = task_order,
  modality_order = MODALITY_ORDER,
  value_col = VALUE_COL,
  task_col = TASK_COL,
  modality_col = MODALITY_COL,
  subject_col = SUBJECT_COL,
  data_path = data_path,
  formula = formula_str
)

saveRDS(
  analysis_params,
  file.path(output_path, "analysis_params_R.rds")
)

saveRDS(
  fit,
  file.path(output_path, "mixedlm_result_R.rds")
)

fixef_df <- data.frame(
  Term = names(fixef(fit)),
  Estimate = as.numeric(fixef(fit))
)

write_csv(
  fixef_df,
  file.path(output_path, "lmm_fixef_R.csv")
)

vcov_mat <- as.matrix(vcov(fit))
vcov_df <- as.data.frame(vcov_mat)
vcov_df <- cbind(Term = rownames(vcov_mat), vcov_df)

write_csv(
  vcov_df,
  file.path(output_path, "lmm_vcov_R.csv")
)

cat("Saved results to: ", output_path, "\n", sep = "")

