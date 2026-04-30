# Project: Feeding behaviour, growth performance and nutrient utilisation of abalone with dietary Ulva sp. supplementation: A meta-analysis 

## Step 4: Running meta-analysis

## Load required packages
packages_cran <- c(
  "here",         
  "readr",         
  "dplyr",         
  "tidyr",         
  "metafor",       
  "ggplot2",       
  "ggcorrplot",    
  "patchwork",     
  "car",
  "clubSandwich",
  "devtools",     
  "emmeans",      
  "R.rsp"         
)

# Install and load any missing CRAN packages
installed_packages <- rownames(installed.packages())
for (p in packages_cran) {
  if (!(p %in% installed_packages)) {
    install.packages(p, dependencies = TRUE)
  }
  library(p, character.only = TRUE)
}

# Install orchaRd from GitHub if not already installed
if (!("orchaRd" %in% installed_packages)) {
  devtools::install_github("daniel1noble/orchaRd", ref = "main", force = TRUE)
}
library(orchaRd)

### Please download GitHub repository and then run the following
here()

# Load data
clean_data <- read_csv(here("Data", "cleaned_data_for_meta_analysis.csv"))
head(clean_data)

# Calculate effect size (lnRR) and variance
clean_data <- clean_data %>%
  mutate(
    lnRR = log(treatment_mean / control_mean),
    vi_lnRR = (treatment_SD^2) / (treatment_n * treatment_mean^2) +
              (control_SD^2) / (control_n * control_mean^2)
  )

# Reverse the sign of FCR effect size (negative result indicates positive biological improvement)
clean_data <- clean_data %>%
  mutate(lnRR = ifelse(outcome == "FCR", -lnRR, lnRR))

# Scale publication year 
clean_data$publication_year_scaled <- as.numeric(scale(clean_data$publication_year))

# Create check to make sure that direction of lnRR are biologically correct
clean_data %>%
  mutate(
    expected_direction = case_when(
      outcome == "FCR" ~ control_mean - treatment_mean, 
      TRUE ~ treatment_mean - control_mean
    ),
    consistency = sign(expected_direction) == sign(lnRR)
  ) %>%
  group_by(outcome) %>%
  summarise(
    agreement = mean(consistency, na.rm = TRUE)
  )

# Check distribution of effect sizes overall and for each outcome category
lnRR_feed    <- clean_data %>% filter(outcome_category == "feed behaviour")
lnRR_growth  <- clean_data %>% filter(outcome_category == "growth performance")
lnRR_nutrient <- clean_data %>% filter(outcome_category == "nutrient utilisation")

# Create histograms
hist(clean_data$lnRR,    breaks = 30, main = "all data")  
hist(lnRR_feed$lnRR,    breaks = 30, main = "Feed behaviour")  
hist(lnRR_growth$lnRR,  breaks = 30, main = "Growth performance")  
hist(lnRR_nutrient$lnRR, breaks = 30, main = "Nutrient utilisation") 

# Check highly negative effect sizes
clean_data %>%
  arrange(desc(abs(lnRR))) %>%
  dplyr::select(ES_ID, short_citation, treatment_name, outcome_long,
         treatment_mean, control_mean,
         treatment_SD, control_SD,
         treatment_n, control_n,
         lnRR, vi_lnRR) %>%
  head(5)  # Flag S004.12, S004.13, S004.4, S001.30 and S004.5 for sensitivity analysis

# Check distribution of variances overall and for each outcome category
vi_lnRR_feed    <- clean_data %>% filter(outcome_category == "feed behaviour")
vi_lnRR_growth  <- clean_data %>% filter(outcome_category == "growth performance")
vi_lnRR_nutrient <- clean_data %>% filter(outcome_category == "nutrient utilisation")

# Create histograms
hist(clean_data$vi_lnRR,     breaks = 30, main = "all data")  
hist(vi_lnRR_feed$vi_lnRR,   breaks = 30, main = "Feed behaviour")  
hist(vi_lnRR_growth$vi_lnRR, breaks = 30, main = "Growth performance")  
hist(vi_lnRR_nutrient$vi_lnRR, breaks = 30, main = "Nutrient utilisation") 

# Check very low vi_lnRR
clean_data %>%
  arrange(vi_lnRR) %>%
  dplyr::select(ES_ID, short_citation, outcome,
         treatment_mean, control_mean,
         treatment_SD, control_SD,
         treatment_n, control_n,
         lnRR, vi_lnRR) %>%
  head(5) # Flag S004.11, S001.18, S006.23, S006.26, S008.3 for sensitivity analysis

# Many ES_ID compare to the same control (e.g., S001 compares three different Ulva doses to the same 0% control)
# Calculate variance covariance (VCV) matrix to account for this
VCV <- vcalc(
  vi = vi_lnRR,
  cluster = cohort_ID,
  obs = ES_ID,
  rho = 0.5,
  data = clean_data
)

# Run overall meta-analysis
res_3L_all <- rma.mv(yi = lnRR, V = VCV,
                     random = ~ 1 | study_ID / ES_ID,
                     test = "t",
                     data = clean_data,
                     method = "REML")
res_3L_all

# Check VCV
dim(res_3L_all$V)        # should be n x n
is.matrix(res_3L_all$V)  # should return TRUE
str(res_3L_all$V)

# Function: extract variance components (3-level model)
  calc_I2_3level <- function(model) {
  sigma <- model$sigma2
  mean_vi <- mean(diag(model$V), na.rm = TRUE)
  total_var <- sum(sigma) + mean_vi
  data.frame(
    I2_study = sigma[1] / total_var * 100,
    I2_ES    = sigma[2] / total_var * 100,
    I2_total = sum(sigma) / total_var * 100
  )
}

# Extract variance components for res_3L_all
I2_all <- calc_I2_3level(res_3L_all)
print(I2_all)

# Function: generate orchard plots
run_orchard_plot <- function(model, I2, 
                             mod = "1",
                             title = "Orchard Plot",
                             xlab = "Effect Size (lnRR)",
                             x_pos = 0.7,
                             y_pos = 0.5,
                             y_limits = c(-1.5, 1.5),
                             colour = "grey") {
  
  p <- orchaRd::orchard_plot(
    model,
    mod = mod,
    xlab = xlab,
    group = "study_ID"
  ) +
    annotate(
      geom  = "text",
      x     = x_pos,
      y     = y_pos,
      label = paste0("italic(I)^2 == ", round(I2$I2_total, 2), "*\"%\""),
      color = "black",
      parse = TRUE,
      size  = 4,
      face  = "bold"
    ) +
    ggtitle(title) +
    scale_fill_manual(values = colour) +
    scale_colour_manual(values = colour) +
    scale_y_continuous(limits = y_limits)
  return(p)
}

# Create and save orchard plot for full dataset
all_data_plot <- run_orchard_plot(model = res_3L_all, I2 = I2_all)
all_data_plot
ggsave(here("Figures", "all_data_orchard_plot.png"), width = 9, height = 8, units = "in")

# Testing effect of species 

# MLMA with species as fixed effect (no-intercept: estimate per level)
res_species_fixed <- rma.mv(
  yi     = lnRR,
  V      = VCV,
  mods   = ~ 0 + species,
  random = ~ 1 | study_ID / ES_ID,
  test   = "t",
  data   = clean_data,
  method = "REML"
)
summary(res_species_fixed)

# R² for species as fixed moderator
r2_res_species_fixed <- r2_ml(res_species_fixed, res_3L_all)
cat("R2 (Meta-regression vs base model): ", r2_res_species_fixed, "\n")

# MLMA with species as random effect
res_species_random <- rma.mv(
  yi     = lnRR,
  V      = VCV,
  random = list(
    ~ 1 | species,
    ~ 1 | study_ID / ES_ID
  ),
  test   = "t",
  data   = clean_data,
  method = "REML"
)
summary(res_species_random)

# Function: extract variance components (4-level model with species)
calc_I2_4level <- function(model) {
  sigma <- model$sigma2
  mean_vi <- mean(diag(model$V), na.rm = TRUE)
  total_var <- sum(sigma) + mean_vi
  data.frame(
    I2_species = sigma[1] / total_var * 100,
    I2_study   = sigma[2] / total_var * 100,
    I2_ES      = sigma[3] / total_var * 100,
    I2_total   = sum(sigma) / total_var * 100
  )
}

I2_species_model <- calc_I2_4level(res_species_random)
print(I2_species_model)

# Compare models — species random effect does not improve fit
anova(res_3L_all, res_species_random)

#### Adding species as a random effect does not improve model fit; continue without

# Publication bias function

# year_var should be the SCALED publication year column (publication_year_scaled)to ensure the intercept is interpretable at the mean publication year.
run_bias_models <- function(data, 
                             yi = "lnRR",
                             vi = "vi_lnRR",
                             vcv = NULL,
                             study_id = "study_ID",
                             es_id = "ES_ID",
                             n_treatment = "treatment_n",
                             n_control = "control_n",
                             year_var = "publication_year_scaled",
                             bio_moderators = NULL) {
  
  # Compute sqrt_inv_e_n
  data$sqrt_inv_e_n <- sqrt(
    (data[[n_treatment]] + data[[n_control]]) /
    (data[[n_treatment]] * data[[n_control]])
  )
  
  # Build random effects formula
  random_fx <- list(
    as.formula(paste0("~1 | ", study_id)),
    as.formula(paste0("~1 | ", es_id))
  )
  
  # Egger's test (uses scalar V, not VCV, as per standard practice)
  egger_mod <- rma.mv(
    yi     = data[[yi]],
    V      = data[[vi]],
    mods   = ~ 1 + sqrt_inv_e_n,
    random = random_fx,
    test   = "t",
    method = "REML",
    sparse = TRUE,
    data   = data
  )
  
  # Decline effect — uses scaled year and VCV
  year_formula <- as.formula(paste0("~ 1 + ", year_var))
  
  year_mod <- rma.mv(
    yi     = data[[yi]],
    V      = if (!is.null(vcv)) vcv else data[[vi]],
    mods   = year_formula,
    random = random_fx,
    test   = "t",
    method = "REML",
    sparse = TRUE,
    data   = data
  )
  
  # Multivariate bias model — biological moderators + Egger + scaled year
  if (!is.null(bio_moderators)) {
    multi_formula <- as.formula(
      paste0("~ ",
             paste(bio_moderators, collapse = " + "),
             " + sqrt_inv_e_n + ", year_var)
    )
    multi_mod <- rma.mv(
      yi     = data[[yi]],
      V      = if (!is.null(vcv)) vcv else data[[vi]],
      mods   = multi_formula,
      random = random_fx,
      test   = "t",
      method = "REML",
      sparse = TRUE,
      data   = data
    )
  } else {
    multi_mod <- NULL
  }
  
  # Bubble plots
  p_egger <- bubble_plot(egger_mod,
                         mod   = "sqrt_inv_e_n",
                         group = study_id,
                         xlab  = "Square root of inverse effective sample size",
                         ylab  = "Effect size (lnRR)")
  
  p_year <- bubble_plot(year_mod,
                        mod   = year_var,
                        group = study_id,
                        xlab  = "Publication year (scaled)",
                        ylab  = "Effect size (lnRR)")
  
  if (!is.null(multi_mod)) {
    p_multi_egger <- bubble_plot(multi_mod,
                                 mod   = "sqrt_inv_e_n",
                                 group = study_id,
                                 xlab  = "Square root of inverse effective sample size",
                                 ylab  = "Effect size (lnRR)")
    p_multi_year <- bubble_plot(multi_mod,
                                mod   = year_var,
                                group = study_id,
                                xlab  = "Publication year (scaled)",
                                ylab  = "Effect size (lnRR)")
  }
  
  list(
    egger_model      = egger_mod,
    year_model       = year_mod,
    multi_model      = multi_mod,
    plot_egger       = p_egger,
    plot_year        = p_year,
    plot_multi_egger = if (!is.null(multi_mod)) p_multi_egger else NULL,
    plot_multi_year  = if (!is.null(multi_mod)) p_multi_year  else NULL
  )
}

# Create quadratic dose term for publication bias models
clean_data$intervention_dose2 <- clean_data$intervention_dose^2

# Run publication bias models — full dataset
# Note: year_var defaults to "publication_year_scaled" inside the function
results_all_data <- run_bias_models(
  data           = clean_data,
  vcv            = VCV,
  year_var       = "publication_year_scaled",
  bio_moderators = c("intervention_dose",
                     "intervention_dose2",
                     "study_duration_days")
)

# Print results — full dataset
summary(results_all_data$egger_model)
summary(results_all_data$year_model)
summary(results_all_data$multi_model)

# Save publication bias plots — full dataset
results_all_data$plot_egger
ggsave(here("Figures", "plot_egger.png"),
       plot = results_all_data$plot_egger, dpi = 300, width = 8, height = 6, units = "in")

results_all_data$plot_year
ggsave(here("Figures", "plot_year.png"),
       plot = results_all_data$plot_year, dpi = 300, width = 8, height = 6, units = "in")

results_all_data$plot_multi_egger
ggsave(here("Figures", "plot_multi_egger.png"),
       plot = results_all_data$plot_multi_egger, dpi = 300, width = 8, height = 6, units = "in")

results_all_data$plot_multi_year
ggsave(here("Figures", "plot_multi_year.png"),
       plot = results_all_data$plot_multi_year, dpi = 300, width = 8, height = 6, units = "in")

# Stacked combined plot — full dataset
combined_plot <- results_all_data$plot_egger + results_all_data$plot_year +
  plot_annotation(
    tag_levels = list(c("A)", "B)")),
    theme = theme(plot.tag = element_text(face = "bold"))
  )
combined_plot
ggsave(here("Figures", "combined_plot.png"),
       plot = combined_plot, dpi = 300, width = 16, height = 8, units = "in")

# Funnel plot — full dataset
png(here("Figures", "funnel_plot_all.png"), width = 8, height = 6, units = "in", res = 600)
funnel(res_3L_all,
       yaxis    = "seinv",
       xlab     = "Standardised residuals (lnRR)",
       ylab     = "Precision (inverse of SE)",
       col      = alpha("black", 0.5),
       cex      = 1.0,
       grid     = FALSE,
       xaxt     = "n",
       yaxt     = "n",
       font.lab = 1)
axis(1, at = axTicks(1), labels = formatC(axTicks(1), digits = 1, format = "f"))
axis(2, at = axTicks(2), labels = formatC(axTicks(2), digits = 1, format = "f"))
dev.off()

# Sensitivity analysis 
# Leave-one-out sensitivity analysis with VCV subsetting
study_list <- unique(clean_data$study_ID)

sensitivity_results_VCV <- lapply(study_list, function(S) {
  
  dat_sub  <- subset(clean_data, study_ID != S)
  keep_idx <- which(clean_data$study_ID != S)
  V_sub    <- VCV[keep_idx, keep_idx]
  
  res_sub <- rma.mv(
    yi     = lnRR,
    V      = V_sub,
    random = ~ 1 | study_ID / ES_ID,
    data   = dat_sub,
    method = "REML"
  )
  
  data.frame(
    study_removed = S,
    est           = res_sub$beta,
    se            = res_sub$se,
    tau2_study    = res_sub$sigma2[1],
    tau2_ES       = res_sub$sigma2[2],
    pct_change    = (100 * (res_sub$beta - res_3L_all$beta) / res_3L_all$beta)
  )
  
}) %>% bind_rows()

print(sensitivity_results_VCV)

# S004 is highly influential (>100% change, sign reversal)
# Author (D. Francis) noted poor performance was linked to poor diet stability
# Re-run MLMA without S004

clean_data_sens <- clean_data %>% filter(study_ID != "S004")

# Scale publication year for sensitive dataset
clean_data_sens$publication_year_scaled <- as.numeric(scale(clean_data_sens$publication_year))

# Create quadratic dose term for sensitive dataset
clean_data_sens$intervention_dose2 <- clean_data_sens$intervention_dose^2

# VCV for sensitive dataset
VCV_sens <- vcalc(
  vi      = vi_lnRR,
  cluster = cohort_ID,
  obs     = ES_ID,
  rho     = 0.5,
  data    = clean_data_sens
)

# 3-level MLMA — sensitive dataset
res_3L_sens <- rma.mv(yi = lnRR, V = VCV_sens,
                      random = ~ 1 | study_ID / ES_ID,
                      test   = "t",
                      data   = clean_data_sens,
                      method = "REML")
res_3L_sens

# Check VCV
dim(res_3L_sens$V)
is.matrix(res_3L_sens$V)
str(res_3L_sens$V)

# Comparison table: full vs sensitive model
sensitivity_compare <- data.frame(
  Model      = c("Full", "S004 removed"),
  k          = c(res_3L_all$k,       res_3L_sens$k),
  Estimate   = c(res_3L_all$beta,    res_3L_sens$beta),
  pval       = c(res_3L_all$pval,    res_3L_sens$pval),
  SE         = c(res_3L_all$se,      res_3L_sens$se),
  CI_low     = c(res_3L_all$ci.lb,   res_3L_sens$ci.lb),
  CI_high    = c(res_3L_all$ci.ub,   res_3L_sens$ci.ub),
  tau2_study = c(res_3L_all$sigma2[1], res_3L_sens$sigma2[1]),
  tau2_ES    = c(res_3L_all$sigma2[2], res_3L_sens$sigma2[2]),
  Q          = c(res_3L_all$QE,      res_3L_sens$QE),
  Q_pval     = c(res_3L_all$QEp,     res_3L_sens$QEp)
)
sensitivity_compare

# Variance components — sensitive model
I2_sens <- calc_I2_3level(res_3L_sens)
print(I2_sens)

# Orchard plot — sensitive dataset
sens_data_plot <- run_orchard_plot(
  model = res_3L_sens,
  I2    = I2_sens,
  title = "B) All Outcome Categories (S004 removed)"
)
sens_data_plot
ggsave(here("Figures", "sens_data_orchard_plot.png"), width = 9, height = 8, units = "in")

# Publication bias — sensitive dataset
results_sens_data <- run_bias_models(
  data           = clean_data_sens,
  vcv            = VCV_sens,
  year_var       = "publication_year_scaled",
  bio_moderators = c("intervention_dose",
                     "intervention_dose2",
                     "study_duration_days")
)

# Print results — sensitive dataset
summary(results_sens_data$egger_model)
summary(results_sens_data$year_model)
summary(results_sens_data$multi_model)

# Save publication bias plots — sensitive dataset
results_sens_data$plot_egger
ggsave(here("Figures", "plot_egger_sens.png"),
       plot = results_sens_data$plot_egger, dpi = 300, width = 8, height = 6, units = "in")

results_sens_data$plot_year
ggsave(here("Figures", "plot_year_sens.png"),
       plot = results_sens_data$plot_year, dpi = 300, width = 8, height = 6, units = "in")

results_sens_data$plot_multi_egger
ggsave(here("Figures", "plot_multi_egger_sens.png"),
       plot = results_sens_data$plot_multi_egger, dpi = 300, width = 8, height = 6, units = "in")

results_sens_data$plot_multi_year
ggsave(here("Figures", "plot_multi_year_sens.png"),
       plot = results_sens_data$plot_multi_year, dpi = 300, width = 8, height = 6, units = "in")

# Stacked combined plot — sensitive dataset
combined_plot_sens <- results_sens_data$plot_egger / results_sens_data$plot_year +
  plot_annotation(
    tag_levels = list(c("A)", "B)")),
    theme = theme(plot.tag = element_text(face = "bold"))
  )
combined_plot_sens
ggsave(here("Figures", "combined_plot_sens.png"),
       plot = combined_plot_sens, dpi = 300, width = 12, height = 14, units = "in")

# Funnel plot — sensitive dataset
png(here("Figures", "funnel_plot_sens.png"), width = 8, height = 6, units = "in", res = 600)
funnel(res_3L_sens,
       yaxis    = "seinv",
       xlab     = "Standardised residuals (lnRR)",
       ylab     = "Precision (inverse of SE)",
       col      = alpha("black", 0.5),
       cex      = 1.0,
       grid     = FALSE,
       xaxt     = "n",
       yaxt     = "n",
       font.lab = 1)
axis(1, at = axTicks(1), labels = formatC(axTicks(1), digits = 1, format = "f"))
axis(2, at = axTicks(2), labels = formatC(axTicks(2), digits = 1, format = "f"))
dev.off()

# Outcome category as moderator (MLMR) 
# Two models are fitted for each dataset
#
#   Model A — no-intercept (mods = ~ 0 + outcome_category):
#     Returns a separate estimate for each level; tests each against zero
#
#   Model B — with-intercept (mods = ~ outcome_category):
#     Returns contrasts between levels; tests whether levels differ from each other.
#     Reference level = feed behaviour (alphabetical default).
#     Re-level to obtain additional pairwise contrasts

##  Full dataset 

# Model A: no-intercept — estimate per outcome category (full dataset)
res_meta_reg <- rma.mv(yi = lnRR, V = VCV,
  mods   = ~ 0 + outcome_category,
  random = list(~ 1 | study_ID / ES_ID),
  test   = "t",
  data   = clean_data,
  method = "REML"
)
res_meta_reg

# Model B: with-intercept — contrasts between outcome categories (full dataset)
# Reference level = feed behaviour
res_meta_reg_contrasts <- rma.mv(yi = lnRR, V = VCV,
  mods   = ~ outcome_category,
  random = list(~ 1 | study_ID / ES_ID),
  test   = "t",
  data   = clean_data,
  method = "REML"
)
res_meta_reg_contrasts

# Intercept                              = feed behaviour estimate
# outcome_categorygrowth performance     = growth vs feed behaviour
# outcome_categorynutrient utilisation   = nutrient utilisation vs feed behaviour

# Re-level to get growth performance as reference (nutrient vs growth contrast)
clean_data$outcome_category <- relevel(factor(clean_data$outcome_category),
                                       ref = "growth performance")

res_meta_reg_contrasts_gp <- rma.mv(yi = lnRR, V = VCV,
  mods   = ~ outcome_category,
  random = list(~ 1 | study_ID / ES_ID),
  test   = "t",
  data   = clean_data,
  method = "REML"
)
res_meta_reg_contrasts_gp

# Intercept                              = growth performance estimate
# outcome_categoryfeed behaviour         = feed behaviour vs growth performance
# outcome_categorynutrient utilisation   = nutrient utilisation vs growth performance

# Restore original factor ordering for downstream code
clean_data$outcome_category <- factor(clean_data$outcome_category,
                                      levels = c("feed behaviour",
                                                 "growth performance",
                                                 "nutrient utilisation"))

# R² and variance components for full-dataset MLMR (no-intercept model)
r2_meta_reg <- r2_ml(res_meta_reg, res_3L_all)
cat("R2 (outcome category moderator, full dataset): ", r2_meta_reg, "\n")

I2_all_data_mlmr <- calc_I2_3level(res_meta_reg)
print(I2_all_data_mlmr)

# Orchard plot — full dataset, categorical (uses with-intercept model)
all_data_mlmr_plot <- run_orchard_plot(
  model = res_meta_reg_contrasts,
  I2    = I2_all_data_mlmr,
  mod   = "outcome_category",
  title = "All data — outcome category"
)
all_data_mlmr_plot
ggsave(here("Figures", "all_data_mlmr_plot.png"), dpi = 500, width = 9, height = 8, units = "in")

## Sensitive dataset 

# Model A: no-intercept — estimate per outcome category (sensitive dataset)
res_meta_reg_sens <- rma.mv(yi = lnRR, V = VCV_sens,
  mods   = ~ 0 + outcome_category,
  random = list(~ 1 | study_ID / ES_ID),
  test   = "t",
  data   = clean_data_sens,
  method = "REML"
)
res_meta_reg_sens

# Model B: with-intercept — contrasts between outcome categories (sensitive dataset)
# Reference level = feed behaviour
res_meta_reg_sens_contrasts <- rma.mv(yi = lnRR, V = VCV_sens,
  mods   = ~ outcome_category,
  random = list(~ 1 | study_ID / ES_ID),
  test   = "t",
  data   = clean_data_sens,
  method = "REML"
)
res_meta_reg_sens_contrasts
# Intercept                              = feed behaviour estimate
# outcome_categorygrowth performance     = growth vs feed behaviour
# outcome_categorynutrient utilisation   = nutrient utilisation vs feed behaviour

# Re-level to get growth performance as reference (sensitive dataset)
clean_data_sens$outcome_category <- relevel(factor(clean_data_sens$outcome_category),
                                            ref = "growth performance")

res_meta_reg_sens_contrasts_gp <- rma.mv(yi = lnRR, V = VCV_sens,
  mods   = ~ outcome_category,
  random = list(~ 1 | study_ID / ES_ID),
  test   = "t",
  data   = clean_data_sens,
  method = "REML"
)
res_meta_reg_sens_contrasts_gp

# Intercept                              = growth performance estimate
# outcome_categoryfeed behaviour         = feed behaviour vs growth performance
# outcome_categorynutrient utilisation   = nutrient utilisation vs growth performance

# Restore original factor ordering for downstream code
clean_data_sens$outcome_category <- factor(clean_data_sens$outcome_category,
                                           levels = c("feed behaviour",
                                                      "growth performance",
                                                      "nutrient utilisation"))

# R² and variance components for sensitive MLMR (no-intercept model)
I2_sens_mlmr <- calc_I2_3level(res_meta_reg_sens)
print(I2_sens_mlmr)

# Orchard plot — sensitive dataset, categorical (uses with-intercept model)
sens_mlmr_plot <- run_orchard_plot(
  model = res_meta_reg_sens_contrasts,
  I2    = I2_sens_mlmr,
  mod   = "outcome_category",
  title = "Sensitive dataset — outcome category (S004 removed)"
)
sens_mlmr_plot
ggsave(here("Figures", "sens_mlmr_plot.png"),
       plot = sens_mlmr_plot, dpi = 300, width = 9, height = 8, units = "in")

# Subgroup MLMA for each outcome category 

# Create subgroup datasets
clean_data_sens_feed   <- clean_data_sens %>% filter(outcome_category == "feed behaviour")
clean_data_sens_growth <- clean_data_sens %>% filter(outcome_category == "growth performance")
clean_data_sens_nutr   <- clean_data_sens %>% filter(outcome_category == "nutrient utilisation")

# Check direction of effect sizes for nutrient utilisation
clean_data_sens %>%
  filter(outcome_category == "nutrient utilisation") %>%
  select(lnRR, control_mean, treatment_mean, outcome_long) %>%
  print(n = 55)

# Function: run subgroup MLMA for a given outcome category
run_mlma <- function(data, outcome_cat, rho = 0.5,
                     random_structure = "~1 | study_ID/ES_ID") {
  
  dat_sub <- data %>% filter(outcome_category == outcome_cat)
  
  if (nrow(dat_sub) < 2) {
    warning(paste("Not enough data for outcome:", outcome_cat))
    return(NULL)
  }
  
  V_sub <- vcalc(
    vi      = vi_lnRR,
    cluster = cohort_ID,
    obs     = ES_ID,
    rho     = rho,
    data    = dat_sub
  )
  
  res <- rma.mv(
    yi     = lnRR,
    V      = V_sub,
    random = as.formula(random_structure),
    test   = "t",
    data   = dat_sub,
    method = "REML"
  )
  
  list(model = res, VCV = res$V, data = dat_sub)
}

# Run subgroup MLMA for all three outcome categories
outcome_list <- c("feed behaviour", "growth performance", "nutrient utilisation")
results_mlma <- lapply(outcome_list, function(cat) run_mlma(clean_data_sens, cat))
names(results_mlma) <- outcome_list

# Print subgroup results
results_mlma[["feed behaviour"]]$model
results_mlma[["growth performance"]]$model
results_mlma[["nutrient utilisation"]]$model

# Extract variance components for all subgroups
I2_results_subgroups <- lapply(results_mlma, function(x) {
  if (is.null(x)) return(NULL)
  calc_I2_3level(x$model)
})
names(I2_results_subgroups) <- outcome_list
I2_results_subgroups

# Moderator collinearity check 
corr_cont <- round(
  cor(clean_data_sens[, c("intervention_dose", "study_duration_days",
                          "initial_size_g", "intervention_dose2")]), 2)

p_cont <- ggcorrplot(corr_cont, hc.order = TRUE, lab = TRUE,
                     outline.col = "white",
                     colors = c("#6D9EC1", "white", "#E46726"),
                     title = "(a) Continuous variables")

corr_cont_log <- round(
  cor(log(clean_data_sens[, c("intervention_dose", "study_duration_days",
                              "initial_size_g", "intervention_dose2")])), 2)

p_cont_log <- ggcorrplot(corr_cont_log, hc.order = TRUE, lab = TRUE,
                         outline.col = "white",
                         colors = c("#6D9EC1", "white", "#E46726"),
                         title = "(b) Log-transformed continuous variables")

p_cont + p_cont_log + plot_layout(guides = "collect")

# VIF check
lm_check <- lm(lnRR ~ study_duration_days + initial_size_g + intervention_dose,
               data = clean_data_sens)
car::vif(lm_check)

### initial_size_g and study_duration_days are moderately correlated
### Fit size and duration models separately, then compare

# MLMR with continuous moderators

fixed_vars_size     <- c("intervention_dose", "initial_size_g", "intervention_dose2")
fixed_vars_duration <- c("study_duration_days", "intervention_dose", "intervention_dose2")
fixed_vars_all      <- c("study_duration_days", "intervention_dose", "intervention_dose2", "initial_size_g")

# Function: run MLMR with specified fixed effects (no-intercept for continuous moderators)
run_mlmr_fe <- function(data, outcome_cat = NULL, fixed_effects = NULL,
                        rho = 0.5, random_structure = "~1 | study_ID/ES_ID") {
  
  if (!is.null(outcome_cat)) {
    dat_sub <- data %>% filter(outcome_category == outcome_cat)
  } else {
    dat_sub <- data
  }
  
  V_sub <- vcalc(
    vi      = vi_lnRR,
    cluster = cohort_ID,
    obs     = ES_ID,
    rho     = rho,
    data    = dat_sub
  )
  
  if (is.null(fixed_effects)) {
    mods_formula <- NULL
  } else {
    mods_formula <- as.formula(paste("~ 0 +", paste(fixed_effects, collapse = " + ")))
  }
  
  res <- rma.mv(
    yi     = lnRR,
    V      = V_sub,
    mods   = mods_formula,
    random = as.formula(random_structure),
    test   = "t",
    data   = dat_sub,
    method = "REML"
  )
  
  list(model = res, VCV = res$V, data = dat_sub)
}

# Overall models — size, duration, combined
size_MLMR     <- run_mlmr_fe(clean_data_sens, outcome_cat = NULL, fixed_effects = fixed_vars_size)
duration_MLMR <- run_mlmr_fe(clean_data_sens, outcome_cat = NULL, fixed_effects = fixed_vars_duration)
all_MLMR      <- run_mlmr_fe(clean_data_sens, outcome_cat = NULL, fixed_effects = fixed_vars_all)

summary(size_MLMR$model)
summary(duration_MLMR$model)
summary(all_MLMR$model)

# R² for overall MLMR models
r2_size_all     <- r2_ml(size_MLMR$model,     res_3L_sens)
r2_duration_all <- r2_ml(duration_MLMR$model, res_3L_sens)
r2_combined_all <- r2_ml(all_MLMR$model,      res_3L_sens)

cat("R2 size MLMR:        ", r2_size_all,     "\n")
cat("R2 duration MLMR:    ", r2_duration_all, "\n")
cat("R2 combined MLMR:    ", r2_combined_all, "\n")

## Dose-response estimates are robust regardless of whether initial_size_g or
## study_duration_days is retained (AIC difference < 2). Inspect by subgroup.

outcome_list_all <- c("overall", "feed behaviour", "growth performance", "nutrient utilisation")

# Duration models by subgroup
results_mlmr_duration <- lapply(outcome_list_all, function(cat) {
  if (cat == "overall") {
    run_mlmr_fe(clean_data_sens, outcome_cat = NULL,  fixed_effects = fixed_vars_duration)
  } else {
    run_mlmr_fe(clean_data_sens, outcome_cat = cat,   fixed_effects = fixed_vars_duration)
  }
})
names(results_mlmr_duration) <- outcome_list_all

results_mlmr_duration$overall$model
results_mlmr_duration$`feed behaviour`$model
results_mlmr_duration$`growth performance`$model
results_mlmr_duration$`nutrient utilisation`$model

# Size models by subgroup
results_mlmr_size <- lapply(outcome_list_all, function(cat) {
  if (cat == "overall") {
    run_mlmr_fe(clean_data_sens, outcome_cat = NULL,  fixed_effects = fixed_vars_size)
  } else {
    run_mlmr_fe(clean_data_sens, outcome_cat = cat,   fixed_effects = fixed_vars_size)
  }
})
names(results_mlmr_size) <- outcome_list_all

results_mlmr_size$overall$model
results_mlmr_size$`feed behaviour`$model
results_mlmr_size$`growth performance`$model
results_mlmr_size$`nutrient utilisation`$model

# Dose-response bubble plot
# Univariate MLMR using quadratic + linear dose terms (with intercept)
res_meta_dose <- rma.mv(yi = lnRR, V = VCV_sens,
  mods   = ~ intervention_dose + I(intervention_dose^2),
  random = list(~ 1 | study_ID / ES_ID),
  test   = "t",
  data   = clean_data_sens,
  method = "REML"
)
res_meta_dose

# Bubble plot — Ulva inclusion level vs effect size
ulva_inclusion <- bubble_plot(res_meta_dose,
                  mod      = "intervention_dose",
                  group    = "study_ID",
                  xlab     = expression(italic(Ulva) ~ "inclusion level (% w/w)"),
                  ylab     = "Effect size (lnRR)",
                  est.lwd  = 1.2,
                  est.col  = "deeppink3") +
  scale_x_continuous(breaks = seq(0, 30, by = 5), limits = c(0, 30)) +
  scale_y_continuous(breaks = seq(floor(min(clean_data_sens$lnRR)),
                                   ceiling(max(clean_data_sens$lnRR)), 0.5)) +
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(color = "black")
  )

ulva_inclusion
ggsave(here("Figures", "ulva_inclusion_relationship.png"),
       plot = ulva_inclusion, width = 10, height = 6, units = "in")
