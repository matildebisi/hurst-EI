library(dplyr)
library(readr)
library(tidyr)
library(lme4)
library(emmeans)
library(lmerTest)
library(sjPlot)
library(ggeffects)

#### _____Config_____ ####

base_path <- normalizePath(getwd(), winslash = "/")
field <- "Hursts"

Hurst_path <- file.path(base_path,"RESULTS/discovery_dataset/hurst/dataframes/df_hurst_parcelwise_FIN_D4_modwt_1_5_padding_0_decimate_0_residuals_w_rest_IQR_cleaned.csv")
MRS_path <- file.path(base_path, "DATA/discovery_dataset/cleaned_MRS_df.csv")
selected_parcels_path <- file.path(base_path, "DATA/discovery_dataset/selected_parcels_overlap30.csv")

output_dir <- file.path(base_path, "RESULTS/discovery_dataset/models")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

EIB_col <-  "GlxGABA_AlphaTissCorr"
voxel_col <- "voxel"
parcel_col <- "parcel_id"

voxel_levels <- c(
  "anteriorcingulate",
  "occipital",
  "rightauditory",
  "smleft"
)

parcel_cyto_path <- file.path(base_path, "DATA/cytoarchitecture/bb_gradient_df.csv")   
gradient_col <- "Cyto_class"


#### _____Dataset preparation_____ ####

#### 1. Load Hurst parcelwise data ####

hurst_df <- read_csv(Hurst_path, show_col_types = FALSE) %>%
  select(-starts_with("Unnamed")) %>%
  mutate(
    Parcel = as.integer(Parcel)
  )

#### 2. Load selected parcels per voxel ####

selected_raw <- read_csv(selected_parcels_path, show_col_types = FALSE)

selected_parcels <- selected_raw %>%
  transmute(
    Voxel = (.data[[voxel_col]]),
    Parcel = as.integer(.data[[parcel_col]])
  ) %>%
  filter(!is.na(Voxel), Voxel != "", !is.na(Parcel)) %>%
  distinct()

cat("\nSelected parcels per voxel:\n")
print(selected_parcels %>% count(Voxel, name = "n_selected_parcels"))

#### 3. Load MRS data ####

mrs_model_df <- read_csv(MRS_path, show_col_types = FALSE)

#### 4. Keep Hurst parcelwise over selected parcels ####

group_cols <- c("Subject", "Voxel", "Parcel")

hurst_voxel_df <- hurst_df %>%
  inner_join(selected_parcels, by = "Parcel") %>%
  filter(Voxel %in% voxel_levels) %>%
  filter(!is.na(.data[[field]])) %>%
  mutate(
    Hursts = .data[[field]]
  ) %>%
  select(Subject, Voxel, Parcel, Hursts, everything())

cat("\nHurst data after selecting parcels, without averaging:\n")
cat("Rows: ", nrow(hurst_voxel_df), "\n", sep = "")
cat("Subjects: ", length(unique(hurst_voxel_df$Subject)), "\n", sep = "")
print(hurst_voxel_df %>% count(Voxel, name = "n_rows"))

cat("\nParcels per voxel actually present in Hurst data:\n")
print(
  hurst_voxel_df %>%
    distinct(Voxel, Parcel) %>%
    count(Voxel, name = "n_parcels")
)


#### 5. Merge with MRS ####

merged_before_na <- hurst_voxel_df %>%
  inner_join(
    mrs_model_df %>%
      filter(Voxel %in% voxel_levels),
    by = c("Subject", "Voxel")
  )

cat("\nRows after inner_join before NA removal = ", nrow(merged_before_na), "\n", sep = "")

merged_all <- merged_before_na %>%
  drop_na(Hursts, all_of(EIB_col)) %>%
  mutate(
    Subject = factor(Subject),
    Voxel = factor(Voxel, levels = voxel_levels),
    Parcel = factor(Parcel)
  )

cat("Rows after NA removal = ", nrow(merged_all), "\n", sep = "")

cat("\nMerged dataset:\n")
cat("N rows     = ", nrow(merged_all), "\n", sep = "")
cat("N subjects = ", length(unique(merged_all$Subject)), "\n", sep = "")
cat("N parcels  = ", length(unique(merged_all$Parcel)), "\n", sep = "")
print(merged_all %>% count(Voxel, name = "n_rows"))

cat("\nRows per subject:\n")
print(summary(merged_all %>% count(Subject) %>% pull(n)))

cat("\nRows per voxel x parcel:\n")
print(
  merged_all %>%
    count(Voxel, Parcel, name = "n_rows") %>%
    arrange(Voxel, Parcel)
)

#### 6. Add parcel gradients / cyto classes ####

parcel_cyto_df <- read_csv(parcel_cyto_path, show_col_types = FALSE) %>%
  transmute(
    Parcel = as.integer(Parcel),
    Gradient = as.numeric(.data[[gradient_col]])
  ) %>%
  distinct(Parcel, .keep_all = TRUE)

merged_all <- merged_all %>%
  mutate(
    Parcel_int = as.integer(as.character(Parcel))
  ) %>%
  left_join(
    parcel_cyto_df,
    by = c("Parcel_int" = "Parcel")
  ) %>%
  select(-Parcel_int) %>%
  relocate(Gradient, .after = Parcel)


cat("\nGradient / cyto-class diagnostics:\n")
cat("Rows with missing Gradients = ", sum(is.na(merged_all$Gradient)), "\n", sep = "")
print(
  merged_all %>%
    distinct(Parcel, Gradient) %>%
    arrange(Parcel)
)



############################### _____Models_____ ###############################

#### Final polynomial models ####

merged_all$EIB_c <- scale(merged_all$GlxGABA_AlphaTissCorr, scale = FALSE)
merged_all$Gradient_c <- scale(merged_all$Gradient, scale = FALSE)

# check best polynomial order
m_linear <- lmer(Hursts ~ EIB_c * Gradient_c + (1 | Subject) + (0 + EIB_c | Parcel), data = merged_all, REML = FALSE)
m_quad <- lmer(Hursts ~ poly(EIB_c, 2, raw = TRUE) * Gradient_c + (1 | Subject) + (0 + EIB_c | Parcel), data = merged_all, REML = FALSE)
m_cube <- lmer(Hursts ~ poly(EIB_c, 3, raw = TRUE) * Gradient_c + (1 | Subject) + (0 + EIB_c | Parcel), data = merged_all, REML = FALSE)
anova(m_linear, m_quad, m_cube)

m_correct <- lmer(Hursts ~ poly(EIB_c, 3, raw = TRUE) * Gradient_c + 
                    (1 | Subject) + (0 + EIB_c | Parcel), 
                  data = merged_all)

summary(m_correct)
plot_model(m_correct, type = "pred", terms = c("EIB_c [all]", "Gradient_c"))


g_agranular <- min(merged_all$Gradient_c, na.rm = TRUE)
g_granular  <- max(merged_all$Gradient_c, na.rm = TRUE)
end_trends_bigbrain <- emtrends(
  m_correct,
  specs = ~ Gradient_c,
  var = "EIB_c",
  at = list(
    Gradient_c = c(g_agranular, g_granular),
    EIB_c = 0
  )
)

summary(
  end_trends_bigbrain,
  infer = c(TRUE, TRUE),
  adjust = "BH"
)

pairs(
  end_trends_bigbrain,
  adjust = "BH"
)

# Selected gradient levels
gradient_levels <- c(-0.08, -0.04, 0.01, 0.05)

# Predictions at selected gradient levels
preds <- ggpredict(
  m_correct,
  terms = c(
    "EIB_c [all]",
    paste0("Gradient_c [", paste(gradient_levels, collapse = ","), "]")
  )
)

# Save predictions
output_path <- file.path(output_dir, "bigbrain_model.csv")
write.csv(as.data.frame(preds), output_path, row.names = FALSE)

raw_df <- merged_all[, c("EIB_c", "Gradient_c", "Hursts", "Subject", "Parcel")]
output_path <- file.path(output_dir, "bigbrain_raw.csv")
write.csv(raw_df, output_path, row.names = FALSE)
