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
  "car"          
)

# Install any missing packages
installed_packages <- rownames(installed.packages())
for (p in packages) {
  if (!(p %in% installed_packages)) {
    install.packages(p, dependencies = TRUE)
  }
}

###Please download GitHub repository and then run the following
here()
clean_data <- read_csv(here("GitHub", "UlvaSupplements_AbaloneMA", "Data", "cleaned_data_for_meta_analysis.csv"))
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
                             y_pos = 1.5,
                             y_limits = c(-2.5, 2.5),
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
ggsave("all_data_orchard_plot.png", width = 11, height = 8, units = "in")  # Save orchard plot

## Testing effect of species
# MLMA with species as fixed effect (no intercept)
res_species_fixed <- rma.mv(
  yi   = lnRR,
  V    = VCV,
  mods = ~ 0 + species,                    
  random = ~ 1 | study_ID / ES_ID,  
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
  data = clean_data,
  method = "REML"
)
summary(res_species_random)

# Compare models
anova(res_3L_all, res_species_random)

#### Adding species as a random effect does not improve model fit. Continue without?

## Testing if outcome is a significant moderator
res_outcome_fixed <- rma.mv(
  yi   = lnRR,
  V    = VCV,
  mods = ~ 0 + outcome_category,                    
  random = ~ 1 | study_ID / ES_ID,  
  data = clean_data,
  method = "REML"
)
summary(res_outcome_fixed)

## Publication bias
# Create function for assessing publication bias
run_pub_bias <- function(data, outcome_name = NULL, random_structure = "~1 | study_ID/ES_ID",
                         se_limit = NULL, x_limits = c(-2, 2), y_limit = NULL) {
  
  # Optional: subset by outcome category
  if(!is.null(outcome_name)) {
    data <- data %>% filter(outcome_category == outcome_name)
  }
  
  # Compute SE, precision, effective sample size
  data <- data %>%
    mutate(
      se_lnRR = sqrt(vi_lnRR),
      precision = 1 / se_lnRR,
      n_eff = (treatment_n * control_n) / (treatment_n + control_n),
      inv_n_eff = 1 / n_eff
    )
  
  # Fit 3-level meta-analysis with precision as moderator
  res <- rma.mv(
    yi = lnRR,
    V  = vi_lnRR,
    mods = ~ precision,
    random = as.formula(random_structure),
    data = data,
    method = "REML"
  )
  
  # Robust small-sample corrected test for the moderator
  cr2_test <- coef_test(
    res,
    vcov = "CR2",
    cluster = data$study_ID,
    test = "Satterthwaite"
  )
  
  # Funnel plot preparation
  data <- data %>% mutate(SE = 1 / precision)
  
  precision_grid <- seq(min(data$precision, na.rm = TRUE),
                        max(data$precision, na.rm = TRUE),
                        length.out = 100)
  
  mean_lnRR <- mean(data$lnRR, na.rm = TRUE)
  
  bounds <- data.frame(
    precision = precision_grid,
    lnRR_upper = mean_lnRR + 1.96 / precision_grid,
    lnRR_lower = mean_lnRR - 1.96 / precision_grid
  )
  
  funnel_plot <- ggplot(data, aes(x = lnRR, y = precision)) +
    geom_ribbon(
      data = bounds,
      aes(y = precision, xmin = lnRR_lower, xmax = lnRR_upper),
      fill = "white",
      colour = "black",
      linetype = "dotted",
      size = 0.5,
      inherit.aes = FALSE
    ) +
    geom_point(shape = 21, fill = "grey", color = "black", size = 3) +
    geom_vline(xintercept = 0, linetype = "dotted", color = "black", size = 1) +
    scale_y_continuous(name = "Precision (1/SE)", limits = y_limit) +
    scale_x_continuous(name = "Effect Size (lnRR)", limits = x_limits) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_line(color = "white", size = 0.5),
      panel.grid.minor.y = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.8),
      axis.ticks.length = unit(0.3, "cm"),
      axis.ticks = element_line(color = "black"),
      plot.title = element_text(face = "bold", hjust = 0.01, size = 18)
    ) +
    ggtitle(ifelse(is.null(outcome_name), "Funnel Plot", paste("Funnel Plot:", outcome_name)))
  
  # Return results
  list(
    model = res,
    cr2_test = cr2_test,
    funnel_plot = funnel_plot
  )
}

# Run publication bias test on clean_data and print results
pubbias_all <- run_pub_bias(clean_data)
pubbias_all$model
pubbias_all$cr2_test
pubbias_all$funnel_plot

# Save funnel plot
ggsave("funnel_plot_all_data.png", width = 11, height = 8, units = "in") 

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

# Recreate funnel plot with S004 highlighted
# Funnel plot
funnel_plot_S004 <- ggplot(clean_data, aes(x = lnRR, y = precision)) +
  geom_ribbon(
    data = bounds,
    aes(y = precision, xmin = lnRR_lower, xmax = lnRR_upper),
    fill = "white",
    colour = "black",
    linetype = "dotted",
    linewidth = 0.5,
    inherit.aes = FALSE
  ) +

  geom_point(aes(fill = ifelse(study_ID == "S004", "red", "grey")),
             shape = 21,
             color = "black",
             size = 3) +

  geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 1) +

  scale_fill_identity() +

  scale_y_continuous(name = "Precision (1/SE)", limits = c(0, 400)) +
  scale_x_continuous(name = "Effect Size (lnRR)", limits = c(-2.0, 2.0),
                     breaks = seq(-5, 5, by = 0.5)) +

  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(color = "white", linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.8),
    axis.ticks.length = unit(0.3, "cm"),
    axis.ticks = element_line(color = "black"),
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 16),
    plot.title = element_text(face = "bold", hjust = 0.01, size = 18)
  ) +
  ggtitle("")

# Print plot
funnel_plot_S004

# Save funnel plot
ggsave("funnel_plot_sens_data.png", width = 11, height = 8, units = "in") 

# S004 is highly influential (> 100% pct_change) and flips the sign of lnRR. Removal improves publication bias risk.
# Author (D. Francis) suggested that poor performance linked to poor diet stability
# Run MLMA without S004?????

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
ggsave("sens_data_orchard_plot.png", width = 11, height = 8, units = "in")  # Save orchard plot

# Re-test publication bias
pubbias_sens <- run_pub_bias(clean_data_sens)
pubbias_sens$model
pubbias_sens$cr2_test
pubbias_sens$funnel_plot

# Removal of S004 reduces publication bias risk

# Save funnel plot
ggsave("funnel_plot_sens_data.png", width = 11, height = 8, units = "in") 

## Test effect of outcome_category as moderator by running MLMR with OC as fixed effect 
# Run MLMR (no intercept model)
res_meta_reg <- rma.mv(yi = lnRR, V  = VCV_sens,
  mods = ~ 0 + outcome_category,        
  random = list(                
    ~ 1 | study_ID / ES_ID),
  data = clean_data_sens,
  method = "REML"
)
res_meta_reg

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

# Create orchard_plots for all subgroup MLMA
orchard_plots_subgroup <- lapply(names(results_mlma), function(cat) {
  
  res <- results_mlma[[cat]]
  if(is.null(res)) return(NULL)
  
  I2 <- calc_I2_3level(res$model)
  
  run_orchard_plot(
    model = res$model,
    I2 = I2,
    title = paste(cat)
  )
})

# Assign letters for subplots and new titles
names(orchard_plots_subgroup) <- names(results_mlma)
orchard_plots_clean <- orchard_plots_subgroup[sapply(orchard_plots_subgroup, function(p) inherits(p, "gg"))]
letters <- paste0(LETTERS[1:length(orchard_plots_clean)], ")")
new_titles <- paste(letters, names(orchard_plots_subgroup))

orchard_plots_titled <- mapply(function(plot, title_new) {
  plot +
    ggtitle(title_new) +
    theme(plot.title = element_text(size = 14)) +
    scale_fill_manual(values = "grey") +
    scale_colour_manual(values = "grey")
}, orchard_plots_clean, new_titles, SIMPLIFY = FALSE)

# Create stacked plot with all outcome categories
multi_panel_orchard <- wrap_plots(orchard_plots_titled, nrow = 1)

multi_panel_orchard   # Print orchard plot
ggsave("multi_panel_orchard.png", width = 15, height = 6, units = "in")  # Save orchard plot

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
lm_check <- lm(lnRR ~ study_duration_days + initial_size_g + intervention_dose + I(intervention_dose^2),
               data = clean_data_sens)
car::vif(lm_check)

### Initial_size_g and study_duration_days moderately correlated. 
## Proceed with multivariate MLMR including Initial_size_g and study_duration_days seperately and then run a sensitivity analysis 

# Set fixed effects
fixed_vars_size <- c("intervention_dose", "initial_size_g", "intervention_dose2")
fixed_vars_duration <- c("study_duration_days", "intervention_dose", "intervention_dose2")

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
    data = dat_sub,
    method = "REML"
  )
  
  list(
    model = res,
    VCV = res$V,
    data = dat_sub
  )
}

# Run both models and compare
# Size model
size_MLMR <- run_mlma_fe(
  data = clean_data_sens,
  outcome_cat = NULL,   # overall model
  fixed_effects = fixed_vars_size
)

# Duration model
duration_MLMR <- run_mlma_fe(
  data = clean_data_sens,
  outcome_cat = NULL,   # overall model
  fixed_effects = fixed_vars_duration
)

# Inspect results
summary(size_MLMR$model)
summary(duration_MLMR$model)

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

## Again, dose-response estimates vary between size_MLMR and duration_MLMR BUT dose effects are robust regardless of if initial_size_g
or study duration is included as a fixed effect

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

# Save plot
ggsave("ulva_dose_relationship.png", width = 8, height = 6, units = "in") 

