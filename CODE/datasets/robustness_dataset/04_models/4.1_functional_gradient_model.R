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

Hurst_path <- file.path(base_path,"RESULTS/robustness_dataset/hurst/dataframes/df_hurst_parcelwise_FIN_D4_modwt_1_5_padding_0_decimate_0_residuals_w_rest_IQR_cleaned.csv")
gradient_path <- file.path(base_path, "DATA/gradients/Margulies_gradients_Schaefer100_7Networks.csv")

output_dir <- file.path(base_path, "RESULTS/robustness_dataset/models")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


#### _____Dataset preparation_____ ####

#### 1. Load Hurst parcelwise data ####

hurst_df_all <- read_csv(Hurst_path, show_col_types = FALSE) %>%
  select(-starts_with("Unnamed")) %>%
  mutate(
    Parcel = as.integer(Parcel),
    ordered = FALSE
  )
hurst_df <- hurst_df_all #%>%

#### 2. Load Gradient data ####

gradient_df <- read_csv(gradient_path, show_col_types = FALSE)

#### 3. Merge ####

merged <- hurst_df  %>%
  mutate(
    Parcel_int = as.integer(as.character(Parcel))
  ) %>%
  left_join(
    gradient_df,
    by = c("Parcel_int" = "Parcel")
  ) %>%
  select(-Parcel_int)


cat("\nGradient :\n")
cat("Rows with missing Gradients = ", sum(is.na(merged)), "\n", sep = "")
print(
  merged %>%
    distinct(Parcel) %>%
    arrange(Parcel)
)


############################### _____Models_____ ###############################

#### Final polynomial models ####

# Gradient1
m_linear <- lmer(Hursts ~ MarguliesGradient1 + (1 | Subject) + (1 | Parcel), data = merged, REML = FALSE)
m_quad <- lmer(Hursts ~ poly(MarguliesGradient1, 2, raw = TRUE) + (1 | Subject) + (1 | Parcel), data = merged, REML = FALSE)
m_cube <- lmer(Hursts ~ poly(MarguliesGradient1, 3, raw = TRUE) + (1 | Subject) + (1 | Parcel), data = merged, REML = FALSE)
anova(m_linear, m_quad, m_cube)

summary(m_quad)
plot_model(m_quad, type = "pred", terms = c("MarguliesGradient1 [all]"))

# save
preds <- ggpredict(m_quad, terms = c("MarguliesGradient1 [all]"))
preds_df <- as.data.frame(preds)
write.csv(preds_df, file.path(output_dir, "margulies_model_quad.csv"), row.names = FALSE)
