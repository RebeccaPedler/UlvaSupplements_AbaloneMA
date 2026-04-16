# Project: Feeding behaviour, growth performance and nutrient utilisation of abalone with dietary Ulva sp. supplementation: A meta-analysis 

## Step4: Running meta-analysis

## Load required packages
# List packages
packages <- c(
  "here",         
  "readr",         
  "dplyr",         
  "tidyr",         
  "metafor",       
  "orchaRd",       
  "ggplot2",       
  "ggcorrplot",    
  "patchwork",     
  "car",
  "clubSandwich"
)

# Install and load any missing packages
installed_packages <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed_packages)) {
    install.packages(p, dependencies = TRUE)
  }
  library(p, character.only = TRUE)
}

### Please download GitHub repository and then run the following
here()
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
lnRR_feed <- clean_data %>% filter(outcome_category == "feed behaviour")
lnRR_growth <- clean_data %>% filter(clean_data$outcome_category == "growth performance")
lnRR_nutrient <- clean_data %>% filter(clean_data$outcome_category == "nutrient utilisation")

# Create histograms
hist(clean_data$lnRR, breaks = 30, main = "all data")  
hist(lnRR_feed$lnRR, breaks = 30, main = "Feed behaviour")  
hist(lnRR_growth$lnRR, breaks = 30, main = "Growth performance")  
hist(lnRR_nutrient$lnRR, breaks = 30, main = "Nutrient utilisation") 

# Check highly negative effect sizes
clean_data %>%
  arrange(desc(abs(lnRR))) %>%
  dplyr::select(ES_ID, short_citation, treatment_name, outcome_long,
         treatment_mean, control_mean,
         treatment_SD, control_SD,
         treatment_n, control_n,
         lnRR, vi_lnRR) %>%
  head(5)  #Flag S004.12, S004.13, S004.4, S001.30 and S004.5 for sensitivity analysis

# Check distribution of variances overall and for each outcome category
vi_lnRR_feed <- clean_data %>% filter(outcome_category == "feed behaviour")
vi_lnRR_growth <- clean_data %>% filter(clean_data$outcome_category == "growth performance")
vi_lnRR_nutrient <- clean_data %>% filter(clean_data$outcome_category == "nutrient utilisation")

# Create histograms
hist(clean_data$vi_lnRR, breaks = 30, main = "all data")  
hist(vi_lnRR_feed$vi_lnRR, breaks = 30, main = "Feed behaviour")  
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
  head(5) #Flag S004.11, S001.18, S006.23, S006.26, S008.3 for sensitivity analysis

# Many ES_ID compare to the same control (e.g., S001 compares three different Ulva doses to the same 0% control)
# Calculate variance covariance (VCV) matrix to account for this
VCV <- vcalc(
  vi = vi_lnRR,
  cluster = cohort_ID,
  obs = ES_ID,
  rho = 0.5,
  data = clean_data
)

## Run meta-analysis
# All clean_data (overall model)
res_3L_all <- rma.mv(yi = lnRR, V = VCV,
                     random = ~ 1 | study_ID / ES_ID,
                     test = "t",
                     data = clean_data,
                     method = "REML")
res_3L_all

# Check VCV
dim(res_3L_all$V)           # should be n x n
is.matrix(res_3L_all$V)     # should return TRUE
str(res_3L_all$V)

## Create function for extracting variance components
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

## Create function for generating orchard plots 
run_orchard_plot <- function(model, I2, 
                             title = "Orchard Plot",
                             xlab = "Effect Size (lnRR)",
                             x_pos = 0.7,
                             y_pos = 0.5,
                             y_limits = c(-1.5, 1.5),
                             colour = "grey") {
  
  p <- orchaRd::orchard_plot(
    model,
    mod = "1",
    xlab = xlab,
    group = "study_ID",
  ) +
    
    annotate(
      geom = "text",
      x = x_pos,
      y = y_pos,
      label = paste0(
      "italic(I)^2 == ", round(I2$I2_total, 2), "*\"%\""
      ),
      color = "black",
      parse = TRUE,
      size = 4
    ) +
    
    ggtitle(title) +
    theme(plot.title = element_text(face = "bold")) +
    scale_fill_manual(values = colour) +
    scale_colour_manual(values = colour) +
    scale_y_continuous(limits = y_limits)
  
  return(p)
}

# Create orchard plot for all clean_data
all_data_plot <- run_orchard_plot(
  model = res_3L_all,
  I2 = I2_all,
  title = "A) All Outcome Categories"
)

all_data_plot   # Print orchard plot
ggsave("all_data_orchard_plot.png", width = 9, height = 8, units = "in")  # Save orchard plot

## Testing effect of species
# MLMA with species as fixed effect (no intercept)
res_species_fixed <- rma.mv(
  yi   = lnRR,
  V    = VCV,
  mods = ~ 0 + species,                    
  random = ~ 1 | study_ID / ES_ID,  
  test = "t",
  data = clean_data,
  method = "REML"
)
summary(res_species_fixed)

# MLMA with species as random effect
res_species_random <- rma.mv(
yi   = lnRR,
V    = VCV,
  random = list(
    ~ 1 | species,                  
    ~ 1 | study_ID / ES_ID 
  ),
  test = "t",
  data = clean_data,
  method = "REML"
)
summary(res_species_random)

# Compare models
anova(res_3L_all, res_species_random)

#### Adding species as a random effect does not improve model fit, continue without

## Publication bias

#Function for publication bias models and plots
run_bias_models <- function(data, 
                             yi = "lnRR",
                             vi = "vi_lnRR",
                             vcv = NULL,
                             study_id = "study_ID",
                             es_id = "ES_ID",
                             n_treatment = "treatment_n",
                             n_control = "control_n",
                             year_var = "publication_year",
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
  
  # Egger's test (scalar V) 
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
  
  # Decline effect (VCV)
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
  
  # Multi-moderator bias model (VCV)
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
                         ylab = "Effect size (lnRR)")
  
  p_year <- bubble_plot(year_mod,
                        mod   = year_var,
                        group = study_id,
                        xlab  = "Publication year",
                        ylab = "Effect size (lnRR)")
  
  if (!is.null(multi_mod)) {
    p_multi_egger <- bubble_plot(multi_mod,
                                 mod   = "sqrt_inv_e_n",
                                 group = study_id,
                                 xlab  = "Square root of inverse effective sample size",                                                         
                                 ylab = "Effect size (lnRR)")
    
    p_multi_year <- bubble_plot(multi_mod,
                                mod   = year_var,
                                group = study_id,
                                xlab  = "Publication year",
                                ylab = "Effect size (lnRR)")
  }
  
  # Return everything
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

# Create quadratic term for intervention_dose
clean_data$intervention_dose2 <- clean_data$intervention_dose^2

results_all_data <- run_bias_models(
  data           = clean_data,
  vcv            = VCV,
  year_var       = "publication_year",
  bio_moderators = c("intervention_dose", 
                     "intervention_dose2",
                     "study_duration_days")
)

# Print results (all data)
summary(results_all_data$egger_model)
summary(results_all_data$year_model)
summary(results_all_data$multi_model)

# Print and save plots (all data)
results_all_data$plot_egger 
ggsave("plot_egger.png", plot = results_all_data$plot_egger, dpi = 300, width = 8, height = 6, units = "in")

results_all_data$plot_year
ggsave("plot_year.png", plot = results_all_data$plot_year, dpi = 300, width = 8, height = 6, units = "in")

results_all_data$plot_multi_egger
ggsave("plot_multi_egger.png", plot = results_all_data$plot_multi_egger, dpi = 300, width = 8, height = 6, units = "in")

results_all_data$plot_multi_year
ggsave("plot_multi_year.png", plot = results_all_data$plot_multi_year, dpi = 300, width = 8, height = 6, units = "in")

# Print funnel plot
funnel(res_3L_all, yaxis = "seinv", 
      xlab = "Standarised residuals",
      ylab = "Precision (inverse of SE)",
      # xlim = c(-4.0, 4.5), 
      # ylim = c(0.01, 60.0),
      col = c(alpha("black", 0.5)),
      cex = 1.0) 

ggsave("funnel_all_data.png", dpi = 300, width = 8, height = 6, units = "in")

## Testing if outcome category is a significant moderator
res_outcome_fixed <- rma.mv(
  yi   = lnRR,
  V    = VCV,
  mods = ~ 0 + outcome_category,                    
  random = ~ 1 | study_ID / ES_ID,  
  data = clean_data,
  method = "REML"
)
summary(res_outcome_fixed)

### Sensitivity analysis
# Leave-one-out sensitivity analysis with VCV subsetting
study_list <- unique(clean_data$study_ID)

sensitivity_results_VCV <- lapply(study_list, function(S) {
  
  # Subset data excluding study S
  dat_sub <- subset(clean_data, study_ID != S)
  
  # Identify rows/effects to keep in the VCV
  keep_idx <- which(clean_data$study_ID != S)
  
  # Subset the VCV matrix accordingly
  V_sub <- VCV[keep_idx, keep_idx]
  
  # Re-run the meta-analysis with the subsetted VCV
  res_sub <- rma.mv(
    yi = lnRR,
    V  = V_sub,
    random = ~1 | study_ID / ES_ID,
    data = dat_sub,
    method = "REML"
  )
  
  # Create summary table
  data.frame(
    study_removed = S,
    est = res_sub$beta,
    se = res_sub$se,
    tau2_study = res_sub$sigma2[1],
    tau2_ES = res_sub$sigma2[2],
    pct_change = (100 * (res_sub$beta - res_3L_all$beta) / res_3L_all$beta)
  )
  
}) %>% bind_rows()

print(sensitivity_results_VCV) #Print sensitivity analysis results

# S004 is highly influential (> 100% pct_change) and flips the sign of lnRR. Removal improves publication bias risk.
# Author (D. Francis) suggested that poor performance linked to poor diet stability
# Run MLMA without S004

# Create dataset 
clean_data_sens <- clean_data %>%
     filter(study_ID != "S004")

# Create new variance covariance (VCV) matrix for clean_data_sens
VCV_sens <- vcalc(
  vi = vi_lnRR,
  cluster = cohort_ID,
  obs = ES_ID,
  rho = 0.5,
  data = clean_data_sens
)

# Run 3-level meta-analysis
res_3L_sens <- rma.mv(yi = lnRR, V = VCV_sens,
                     random = ~ 1 | study_ID / ES_ID,
                     test = "t",
                     data = clean_data_sens,
                     method = "REML")
res_3L_sens

# Check VCV
dim(res_3L_sens$V)           # should be n x n
is.matrix(res_3L_sens$V)     # should return TRUE
str(res_3L_sens$V)

# Create summary table comparing sensitive model (S004 removed) and full model
sensitivity_compare <- data.frame(
  Model = c("Full", "S004 removed"),
  k = c(res_3L_all$k, res_3L_sens$k),
  Estimate = c(res_3L_all$beta, res_3L_sens$beta),
  pval = c(res_3L_all$pval, res_3L_sens$pval),
  SE = c(res_3L_all$se, res_3L_sens$se),
  CI_low = c(res_3L_all$ci.lb, res_3L_sens$ci.lb),
  CI_high = c(res_3L_all$ci.ub, res_3L_sens$ci.ub),
  tau2_study = c(res_3L_all$sigma2[1], res_3L_sens$sigma2[1]),
  tau2_ES = c(res_3L_all$sigma2[2], res_3L_sens$sigma2[2]),
  Q = c(res_3L_all$QE, res_3L_sens$QE),
  Q_pval = c(res_3L_all$QEp, res_3L_sens$QEp)
)

sensitivity_compare # Print comparison table

# Extract variance components for res_3L_sens
I2_sens <- calc_I2_3level(res_3L_sens)
print(I2_sens)

# Create orchard plot for sensitivity data
sens_data_plot <- run_orchard_plot(
  model = res_3L_sens,
  I2 = I2_sens,
  title = "B) All Outcome Categories (S004 removed)"
)

sens_data_plot   # Print orchard plot
ggsave("sens_data_orchard_plot.png", width = 9, height = 8, units = "in")  # Save orchard plot

# Publication bias with S004 removed
results_sens_data <- run_bias_models(
  data           = clean_data_sens,
  vcv            = VCV_sens,
  year_var       = "publication_year",
  bio_moderators = c("intervention_dose", 
                     "intervention_dose2",
                     "study_duration_days")
)

# Print results (S004 removed)
summary(results_sens_data$egger_model)
summary(results_sens_data$year_model)
summary(results_sens_data$multi_model)

# Print and save plots (all data)
results_sens_data$plot_egger 
ggsave("plot_egger_sens.png", plot = results_sens_data$plot_egger_sens, dpi = 300, width = 8, height = 6, units = "in")

results_all_data$plot_year
ggsave("plot_year_sens.png", plot = results_sens_data$plot_year_sens, dpi = 300, width = 8, height = 6, units = "in")

results_all_data$plot_multi_egger
ggsave("plot_multi_egger_sens.png", plot = results_sens_data$plot_multi_egger_sens, dpi = 300, width = 8, height = 6, units = "in")

results_all_data$plot_multi_year
ggsave("plot_multi_year_sens.png", plot = results_sens_data$plot_multi_year_sens, dpi = 300, width = 8, height = 6, units = "in")

# Create funnel plot with S004 excluded
funnel(res_3L_sens, yaxis = "seinv", 
      xlab = "Standarised residuals",
      ylab = "Precision (inverse of SE)",
      # xlim = c(-4.0, 4.5), 
      # ylim = c(0.01, 60.0),
      col = c(alpha("black", 0.5)),
      cex = 1.0) 
ggsave("funnel_sens.png", dpi = 300, width = 8, height = 6, units = "in")

## Test effect of outcome_category as moderator by running MLMR with OC as fixed effect (for both full dataset and with S004 removed)
# Run MLMR (no intercept model)

# Full dataset
res_meta_reg <- rma.mv(yi = lnRR, V  = VCV,
  mods = ~ 0 + outcome_category,        
  random = list(                
    ~ 1 | study_ID / ES_ID),
  test = "t",
  data = clean_data,
  method = "REML"
)
res_meta_reg

# Extract variance components for res_meta_reg_sens
I2_sens_all_data_mlmr <- calc_I2_3level(res_meta_reg)
print(I2_sens_all_data_mlmr)

# Create orchard plot for sensitivity data
all_data_mlmr_plot <- run_orchard_plot(
  model = res_meta_reg,
  I2 = I2_sens_all_data_mlmr,
  title = ""
)

all_data_mlmr_plot
ggsave("all_data_mlmr_plot.png", dpi = 300, width = 9, height = 8, units = "in")

# S004 removed
res_meta_reg_sens <- rma.mv(yi = lnRR, V  = VCV_sens,
  mods = ~ 0 + outcome_category,        
  random = list(                
    ~ 1 | study_ID / ES_ID),
  test = "t",
  data = clean_data_sens,
  method = "REML"
)
res_meta_reg_sens

# Extract variance components for res_meta_reg_sens
I2_sens_mlmr <- calc_I2_3level(res_meta_reg_sens)
print(I2_sens_mlmr)

# Create orchard plot for sensitivity data
sens_mlmr_plot <- run_orchard_plot(
  model = res_meta_reg_sens,
  I2 = I2_sens_mlmr,
  title = ""
)

sens_mlmr_plot
ggsave("sens_mlmr_plot.png", plot = sens_mlmr_plot, dpi = 300, width = 9, height = 8, units = "in")

### Outcome category is a significant moderator
## Create sub-groups for each outcome category and run MLMA

# Create subgroup datasets
clean_data_sens_feed <- clean_data_sens %>%
     filter(outcome_category == "feed behaviour")  # Feed behaviour dataset
clean_data_sens_growth <- clean_data_sens %>%
     filter(outcome_category == "growth performance")   # Growth performance dataset
clean_data_sens_nutr <- clean_data_sens %>%
     filter(outcome_category == "nutrient utilisation")    # Nutrient utilisation dataset

# Create function for subgroup MLMA and loop through outcome categories
run_mlma <- function(data, outcome_cat, rho = 0.5, random_structure = "~1 | study_ID/ES_ID") {
  
  # Subset data
  dat_sub <- data %>% filter(outcome_category == outcome_cat)
  
  # Check if there are enough rows
  if(nrow(dat_sub) < 2) {
    warning(paste("Not enough data for outcome:", outcome_cat))
    return(NULL)
  }
  
  # Create variance-covariance matrix
  V_sub <- vcalc(
    vi = vi_lnRR,
    cluster = cohort_ID,
    obs = ES_ID,
    rho = rho,
    data = dat_sub
  )
  
  # Run 3-level meta-analysis
  res <- rma.mv(
    yi = lnRR,
    V  = V_sub,
    random = as.formula(random_structure),
    test = "t",
    data = dat_sub,
    method = "REML"
  )
  
  # Return both the result and the VCV
  list(
    model = res,
    VCV = res$V,
    data = dat_sub
  )
}

# Define outcome categories
outcome_list <- c("feed behaviour", "growth performance", "nutrient utilisation")

# Run MLMA for all outcomes
results_mlma <- lapply(outcome_list, function(cat) run_mlma(clean_data_sens, cat))

# Name the list elements
names(results_mlma) <- outcome_list

# Print results for each outcome category subgroup
results_mlma[["feed behaviour"]]$model
results_mlma[["growth performance"]]$model
results_mlma[["nutrient utilisation"]]$model

# Extract variance components for all subgroup MLMA
I2_results_subgroups <- lapply(results_mlma, function(x) {
  
  if(is.null(x)) return(NULL)  # handle cases with insufficient data
  
  calc_I2_3level(x$model)
})

names(I2_results_subgroups) <- outcome_list

# Print results
I2_results_subgroups

## Assess correlation between continious moderators before choosing variables for inclusion in MLMA 
# Correlation between continuous moderators (code has been adapted from Ayumi et al. (https://zenodo.org/records/13147019))

# Create quadratic term for intervention_dose
clean_data_sens$intervention_dose2 <- clean_data_sens$intervention_dose^2
head(clean_data_sens[, c("intervention_dose", "intervention_dose2")]) #Check that this worked

corr_cont <- round(cor(clean_data_sens[, c("intervention_dose", "study_duration_days", "initial_size_g", "intervention_dose2")]),2)

p_cont <- ggcorrplot(corr_cont, hc.order = TRUE, lab = TRUE, 
                    outline.col = "white", colors = c("#6D9EC1", "white", "#E46726"), 
                    title = "(a) Continuous variables")

corr_cont_log <- round(cor(log(clean_data_sens[, c("intervention_dose", "study_duration_days", "initial_size_g", "intervention_dose2")])),2)

p_cont_log <- ggcorrplot(corr_cont_log, hc.order = TRUE, lab = TRUE, 
                    outline.col = "white", colors = c("#6D9EC1", "white", "#E46726"), 
                    title = "(b) Log-transormed continuous variables")

p_cont + p_cont_log + plot_layout(guides = 'collect')

# Simple linear model for VIF check
lm_check <- lm(lnRR ~ study_duration_days + initial_size_g + intervention_dose,
               data = clean_data_sens)
car::vif(lm_check)

### Initial_size_g and study_duration_days moderately correlated. 
## Proceed with multivariate MLMR including Initial_size_g and study_duration_days seperately and then run a sensitivity analysis 

# Set fixed effects
fixed_vars_size <- c("intervention_dose", "initial_size_g", "intervention_dose2")
fixed_vars_duration <- c("study_duration_days", "intervention_dose", "intervention_dose2")
fixed_vars_all <- c("study_duration_days", "intervention_dose", "intervention_dose2", "initial_size_g")

# Create function for running MLMR with whole dataset and each outcome category (no intercept model)
run_mlmr_fe <- function(data, outcome_cat = NULL, fixed_effects = NULL, 
                        rho = 0.5, random_structure = "~1 | study_ID/ES_ID") {
  
  # Subset data
  if(!is.null(outcome_cat)) {
    dat_sub <- data %>% filter(outcome_category == outcome_cat)
  } else {
    dat_sub <- data
  }
  
  # Create VCV
  V_sub <- vcalc(
    vi = vi_lnRR,
    cluster = cohort_ID,
    obs = ES_ID,
    rho = rho,
    data = dat_sub
  )
  
  # Build moderators
  if(is.null(fixed_effects)) {
    mods_formula <- NULL
  } else {
    mods_formula <- as.formula(paste("~ 0 +", paste(fixed_effects, collapse = " + ")))
  }
  
  # Run model
  res <- rma.mv(
    yi = lnRR,
    V = V_sub,
    mods = mods_formula,
    random = as.formula(random_structure),
    test = "t",
    data = dat_sub,
    method = "REML"
  )
  
  list(
    model = res,
    VCV = res$V,
    data = dat_sub
  )
}

# Run all three models and compare
# Size model
size_MLMR <- run_mlmr_fe(
  data = clean_data_sens,
  outcome_cat = NULL,   # overall model
  fixed_effects = fixed_vars_size
)

# Duration model
duration_MLMR <- run_mlmr_fe(
  data = clean_data_sens,
  outcome_cat = NULL,   # overall model
  fixed_effects = fixed_vars_duration
)

# All model
all_MLMR <- run_mlmr_fe(
  data = clean_data_sens,
  outcome_cat = NULL,   # overall model
  fixed_effects = fixed_vars_all
)

# Inspect results
summary(size_MLMR$model)
summary(duration_MLMR$model)
summary(all_MLMR$model)

## Dose-response estimates vary between size_MLMR and duration_MLMR BUT dose effects are robust
# Difference in AIC < 2 between models
# Inspect in MLMR for each individual subgroup

outcome_list <- c("overall", "feed behaviour", "growth performance", "nutrient utilisation")

results_mlmr_duration <- lapply(outcome_list, function(cat) {
  
  if(cat == "overall") {
    run_mlmr_fe(clean_data_sens, outcome_cat = NULL, fixed_effects = fixed_vars_duration)
  } else {
    run_mlmr_fe(clean_data_sens, outcome_cat = cat, fixed_effects = fixed_vars_duration)
  }
  
})

names(results_mlmr_duration) <- outcome_list

results_mlmr_duration$overall$model
results_mlmr_duration$`feed behaviour`$model
results_mlmr_duration$`growth performance`$model
results_mlmr_duration$`nutrient utilisation`$model

results_mlmr_size <- lapply(outcome_list, function(cat) {
  
  if(cat == "overall") {
    run_mlmr_fe(clean_data_sens, outcome_cat = NULL, fixed_effects = fixed_vars_size)
  } else {
    run_mlmr_fe(clean_data_sens, outcome_cat = cat, fixed_effects = fixed_vars_size)
  }
  
})

names(results_mlmr_size) <- outcome_list

results_mlmr_size$overall$model
results_mlmr_size$`feed behaviour`$model
results_mlmr_size$`growth performance`$model
results_mlmr_size$`nutrient utilisation`$model

## Determining heterogeneity explained by MLMR models
# Calculate r2 for overall model and sub-group models
r2_overall  <- r2_ml(results_mlmr_duration$overall$model, 
                     results_mlma$overall$model)

r2_feed     <- r2_ml(results_mlmr_duration$`feed behaviour`$model, 
                     results_mlma[["feed behaviour"]]$model)

r2_growth   <- r2_ml(results_mlmr_duration$`growth performance`$model, 
                     results_mlma[["growth performance"]]$model)

r2_nutrient <- r2_ml(results_mlmr_duration$`nutrient utilisation`$model, 
                     results_mlma[["nutrient utilisation"]]$model)

# Print results

cat("R2 Overall:              ", r2_overall, "\n")
cat("R2 Feed behaviour:       ", r2_feed, "\n")
cat("R2 Growth performance:   ", r2_growth, "\n")
cat("R2 Nutrient utilisation: ", r2_nutrient, "\n")

## Again, dose-response estimates vary between size_MLMR and duration_MLMR BUT dose effects are robust regardless of if initial_size_g or study duration is included as a fixed effect

## Plot the relationship between Ulva dose and effect size (lnRR)
# Run univariate MLMR (no intercept model) for clean_data_sens with intervention dose as fixed effect
res_meta_dose <- rma.mv(yi = lnRR, V  = VCV_sens,
  mods = ~ 0 + intervention_dose,        
  random = list(                
    ~ 1 | study_ID / ES_ID),
  data = clean_data_sens,
  method = "REML"
)
res_meta_dose

# Create plot
ulva_dose <- bubble_plot(res_meta_dose,           
                  mod = "intervention_dose", 
                  group = "study_ID",       
                  xlab = "Ulva dose (% w/w)", 
                  ylab = "Effect size (lnRR)",
                  est.lwd = 1.2,
                  est.col = "deeppink3") +
       scale_x_continuous(breaks = seq(min(clean_data_sens$intervention_dose),
                                       max(clean_data_sens$intervention_dose),
                                       length.out = 6)) +
       scale_y_continuous(breaks = seq(floor(min(clean_data_sens$lnRR)),
                                       ceiling(max(clean_data_sens$lnRR)), 1)) +
       theme_minimal() +
       theme(
         panel.grid.major = element_blank(),  
         panel.grid.minor = element_blank(), 
         axis.line = element_line(color = "black")  
       )

# Print plot
ulva_dose

# Save ulva_dose_plot
ggsave("ulva_dose_relationship.png", width = 8, height = 6, units = "in") 
