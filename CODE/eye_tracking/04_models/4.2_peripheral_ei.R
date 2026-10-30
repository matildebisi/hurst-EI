library(dplyr)
library(readr)
library(lme4)
library(lmerTest)

#### Config ####
measures <- c("EIB", "GABA", "Glx")

TASK_COL <- "Task"
MODALITY_COL <- "Modality"
SUBJECT_COL <- "Subject"
VALUE_COL <- "Hurst"
ORDER_COL <- "order"

base_path <- normalizePath(getwd(), winslash = "/")
data_type <- "pupil_size"
decimation <- "decimated_2s"

Hurst_path <- file.path(
  base_path,
  "RESULTS/pupillometry",
  data_type,
  "hurst/dataframes",
  paste0("df_hurst_", data_type, "_FIN_D4_modwt_1_5_", decimation, ".csv")
)

output_dir <- file.path(
  base_path,
  "RESULTS/pupillometry",
  data_type,
  "/models/MRS/decimated_2s/range_1_5"
)


MRS_path <- file.path(
  base_path,
  "DATA/replication_dataset/EIB_12s_MRS_static_&_AUC_New_criteria.csv"
)

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#### Load MRS data ####
stat_df <- read_csv(MRS_path, show_col_types = FALSE)

stat_df <- stat_df %>%
  rename(Task = Runs) %>%
  mutate(
    Task = recode(
      Task,
      "ARest"  = "Rest",
      "B0Back" = "0Back",
      "C1Back" = "1Back",
      "D2Back" = "2Back"
    )
  )

ids <- stat_df$Subject
num_sub <- length(ids)

#### Helper functions ####

sig_stars <- function(p) {
  ifelse(
    p < 0.001, "***",
    ifelse(
      p < 0.01, "**",
      ifelse(p < 0.05, "*", "")
    )
  )
}

make_key <- function(measure) {
  paste0("measure_", measure)
}

extract_fixed_effect_table <- function(fit) {
  coef_mat <- summary(fit)$coefficients
  
  coef_df <- as.data.frame(coef_mat)
  coef_df$Term <- rownames(coef_df)
  rownames(coef_df) <- NULL
  
  p_col <- grep("Pr\\(", colnames(coef_df), value = TRUE)
  
  if (length(p_col) == 0) {
    stop("No p-value column found. Make sure lmerTest is loaded.")
  }
  
  out <- data.frame(
    Term = coef_df$Term,
    Coef = coef_df$Estimate,
    `Std Err` = coef_df$`Std. Error`,
    z = coef_df$`t value`,
    `P (raw)` = coef_df[[p_col[1]]],
    check.names = FALSE
  )
  
  return(out)
}

#### Load Hurst data ####
hurst_df <- read_csv(Hurst_path, show_col_types = FALSE)

hurst_df <- hurst_df %>%
  select(-starts_with("Unnamed")) %>%
  filter(.data[[MODALITY_COL]] == "fMRS") %>%
  filter(.data[[SUBJECT_COL]] %in% ids)

#### Merge MRS data ####
merged_base <- hurst_df %>%
  left_join(
    stat_df %>% select(Subject, Task, EIB, Glx, GABA),
    by = c("Subject", "Task")
  )

if ("Rest" %in% unique(merged_base$Task)) {
  merged_base$Task <- factor(
    merged_base$Task,
    levels = c("Rest", setdiff(sort(unique(merged_base$Task)), "Rest"))
  )
} else {
  merged_base$Task <- factor(merged_base$Task)
}

merged_base <- merged_base %>%
  mutate(
    Subject = factor(Subject)
  )

#### 1. Fit models and collect p-values ####

all_results <- list()

for (measure in measures) {
  
  model_vars <- c(VALUE_COL, measure, TASK_COL, SUBJECT_COL)
  
  merged <- merged_base %>%
    tidyr::drop_na(all_of(model_vars))
  
  formula_str <- paste0(
    VALUE_COL, " ~ ", measure, " * Task + ",
    "(1 | Subject)"
  )
  
  model_formula <- as.formula(formula_str)
  
  fit <- lmer(
    model_formula,
    data = merged,
    REML = FALSE,
    control = lmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  )
  
  cat("\nModel diagnostics for measure:", measure, "\n")
  
  cat("\nRandom effects variance:\n")
  print(VarCorr(fit), comp = "Variance")
  
  cat("\nIs singular:", isSingular(fit, tol = 1e-4), "\n")
  
  cat("\nObservations per subject:\n")
  print(
    merged %>%
      count(Subject, name = "n_obs") %>%
      summarise(
        n_subjects = n(),
        min_obs = min(n_obs),
        median_obs = median(n_obs),
        max_obs = max(n_obs)
      )
  )
  
  cat("\nTask counts:\n")
  print(
    merged %>%
      count(Task, name = "n_obs")
  )
  
  summary_table_raw <- extract_fixed_effect_table(fit)
  
  key <- make_key(measure)
  
  all_results[[key]] <- list(
    measure = measure,
    data_type = data_type,
    decimation = decimation,
    formula = formula_str,
    model = fit,
    df = merged,
    summary_table_raw = summary_table_raw,
    pvalues_raw = summary_table_raw$`P (raw)`
  )
}

#### 2. FDR correction across the 3 models ####

family_keys <- sapply(
  measures,
  function(m) make_key(m)
)

pval_concat <- unlist(
  lapply(family_keys, function(k) all_results[[k]]$pvalues_raw)
)

pval_fdr_concat <- p.adjust(pval_concat, method = "BH")

idx <- 1

for (k in family_keys) {
  n <- length(all_results[[k]]$pvalues_raw)
  
  all_results[[k]]$pvalues_fdr <- pval_fdr_concat[idx:(idx + n - 1)]
  
  idx <- idx + n
}

#### 3. Print summary tables and save results ####

for (key in names(all_results)) {
  
  r <- all_results[[key]]
  
  summary_table <- r$summary_table_raw
  summary_table$`P (FDR)` <- r$pvalues_fdr
  summary_table$`Sig (FDR)` <- sig_stars(summary_table$`P (FDR)`)
  
  r$summary_table <- summary_table
  all_results[[key]] <- r
  
  cat("\n", paste(rep("=", 70), collapse = ""), "\n", sep = "")
  cat(
    "measure=", r$measure,
    "  |  data_type=", r$data_type,
    "  |  decimation=", r$decimation,
    "\n",
    sep = ""
  )
  cat(paste(rep("=", 70), collapse = ""), "\n", sep = "")
  print(summary_table, digits = 4, row.names = FALSE)
  
  #### Model fit indices: AIC / BIC  ####
  
  model_fit_df <- data.frame(
    formula = r$formula,
    AIC = AIC(r$model),
    BIC = BIC(r$model)
  )
  
  cat("\nModel fit indices:\n")
  print(model_fit_df, row.names = FALSE, digits = 4)
  cat("\n")
  
  output_name_rds <- paste0(
    "lmm_results_",
    r$measure,
    "_R.rds"
  )
  
  output_path_rds <- file.path(output_dir, output_name_rds)
  saveRDS(r, output_path_rds)
  
  output_name_csv <- paste0(
    "lmm_summary_",
    r$measure,
    "_R.csv"
  )
  
  output_path_csv <- file.path(output_dir, output_name_csv)
  write_csv(summary_table, output_path_csv)
  
  cat("Saved RDS: ", output_path_rds, "\n", sep = "")
  cat("Saved CSV: ", output_path_csv, "\n", sep = "")
  
  df_name_csv <- paste0(
    "lmm_df_",
    r$measure,
    "_R.csv"
  )
  
  write_csv(r$df, file.path(output_dir, df_name_csv))
  
  fixef_df <- data.frame(
    Term = names(fixef(r$model)),
    Estimate = as.numeric(fixef(r$model))
  )
  
  fixef_name_csv <- paste0(
    "lmm_fixef_",
    r$measure,
    "_R.csv"
  )
  
  write_csv(fixef_df, file.path(output_dir, fixef_name_csv))
  
  vcov_mat <- as.matrix(vcov(r$model))
  vcov_df <- as.data.frame(vcov_mat)
  vcov_df <- cbind(Term = rownames(vcov_mat), vcov_df)
  
  vcov_name_csv <- paste0(
    "lmm_vcov_",
    r$measure,
    "_R.csv"
  )
  
  write_csv(vcov_df, file.path(output_dir, vcov_name_csv))
}

all_output_path <- file.path(
  output_dir,
  paste0("all_lmm_results_", data_type, "_", decimation, "_R.rds")
)

saveRDS(all_results, all_output_path)

cat("\nSaved all results to: ", all_output_path, "\n", sep = "")

