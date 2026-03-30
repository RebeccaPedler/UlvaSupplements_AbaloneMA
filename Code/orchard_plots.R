install.packages(c("tidyr", "tidyverse", "rnaturalearth", "rnaturalearthdata", "sf", "ggplot2", "grid", "dplyr", "patchwork", "here", 
"meta", "metafor", "clubSandwich", "robumeta", "devtools", "MASS"))

# Load necessary libraries
library(dplyr)
library(tidyr)
library(tidyverse)
library(rnaturalearth)   # for map data
library(rnaturalearthdata)
library(sf)              # for spatial handling
library(ggplot2)
library(grid)
library(patchwork)
library(here)
library(metafor)
library(meta)
library(clubSandwich)
library(robumeta)
library(MASS)

#install.packages("pacman")
rm(list = ls())
devtools::install_github("daniel1noble/orchaRd", ref = "main", force = TRUE)
pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, orchaRd, emmeans, ape, phytools, flextable)


###Please download GitHub repository and then run the following
setwd("C:/Users/RebeccaPedler/OneDrive - Yumbah/Documents/R&D/Industry PhD/Trials/Meta analysis/R_datasets")
mydata <- read.csv("ulva_meta_analysis_data.csv")

#Add another column to assign short citation
citations <- data.frame(
  study_ID = c("S001", "S002", "S003", "S004", "S005", "S006", "S007", "S008","S009", "S010", "S011"),
  short_citation = c(
    "Bansemer et al. (2016)",
    "Lange et al. (2014)",
    "Stone et al. (2022)",
    "Francis et al. (2021)",
    "Mwangudza (2024)",
    "Falade (2023)",
    "Bates et al. (2017)",
    "Chi et al. (2018)",
    "Angell et al. (2012)",
    "Searle et al. (2025)", 
    "Duong et al. (2021)" 
  ))

mydata <- mydata %>%
  left_join(citations, by = "study_ID")

#Add long outcome names
mydata <- mydata %>%
  mutate(
    outcome_long = case_when(
      outcome == "FI"  ~ "Feed intake",
      outcome == "BG"  ~ "Biomass gain",
      outcome == "CF"  ~ "Condition factor",
      outcome == "FW"  ~ "Final weight",
      outcome == "FSL" ~ "Final shell length",
      outcome == "SG"  ~ "Shell gain",
      outcome == "SGR" ~ "Specific growth rate",
      outcome == "WG"  ~ "Weight gain",
      outcome == "EER" ~ "Energy efficiency ratio",
      outcome == "ED"  ~ "Energy deposition",
      outcome == "FCR" ~ "Feed conversion ratio",
      outcome == "PER" ~ "Protein efficiency ratio",
      outcome == "PD"  ~ "Protein deposition",
      outcome == "LG"  ~ "Length gain",
      outcome == "I" ~ "Ingested feed energy", 
      outcome == "Pg" ~ "Somatic feed energy", 
      outcome == "S" ~ "Shell growth energy",
      TRUE ~ outcome
    )
  )

# Calculate effect size (Hedges g) and pooled SD
all_data <- mydata %>%
  mutate(
    s_pooled = sqrt(((treatment_n - 1) * treatment_SD^2 + (control_n - 1) * control_SD^2) /
                    (treatment_n + control_n - 2)),
    J_all = 1 - (3 / (4 * (treatment_n + control_n) - 9)),
    g_all = J_all * ((treatment_mean - control_mean) / s_pooled),
    vi_g_all = (treatment_n + control_n) / (treatment_n * control_n) + (g_all^2 / (2 * (treatment_n + control_n - 2)))
  )

#Reverse the sign of FCR effect size
all_data <- all_data %>%
  mutate(g_all = ifelse(outcome == "FCR", -g_all, g_all))

# MLMA with species as random effect
res_species_random <- rma.mv(
yi   = g_all,
V    = vi_g_all,
  random = list(
    ~ 1 | species,                  
    ~ 1 | study_ID / experiment_ID 
  ),
  data = all_data,
  method = "REML"
)

# Extract variance components
sigma_species_all <- res_species_random$sigma2[1]
sigma_study_all   <- res_species_random$sigma2[2]
sigma_exp_all    <- res_species_random$sigma2[3]

# Mean sampling variance
mean_vi_all <- mean(all_data$vi_g_all)

# Total variance
total_var_all <- sigma_species_all + sigma_study_all + sigma_exp_all + mean_vi_all

# I² calculations
I2_species_all <- sigma_species_all / total_var_all * 100
I2_study_all   <- sigma_study_all   / total_var_all * 100
I2_exp_all    <- sigma_exp_all     / total_var_all * 100
I2_total_all   <- (sigma_species_all + sigma_study_all + sigma_exp_all) / total_var_all * 100

I2_species_all; I2_study_all; I2_exp_all; I2_total_all

#Plot orchard plot
all_plot <- orchaRd::orchard_plot(
    res_species_random,
    mod = "1",
    xlab = "Effect Size (Hedges g)",
    group = "study_ID"
) +
  annotate(
    geom = "text",
    x = 0.7,
    y = 12,
    label = paste0("italic(I)^{2} == ", round(I2_total_all[1], 2), "*\"%\""),
    color = "black",
    parse = TRUE,
    size = 3
  ) +
  scale_fill_manual(values = "grey") +
  scale_colour_manual(values = "grey")

all_data_plot <- orchaRd::orchard_plot(
  res_species_random,
  mod = "1",
  xlab = "Effect Size (Hedges g)",
  group = "study_ID"
) +
  annotate(
    geom = "text",
    x = 0.7,
    y = 12,
    label = paste0("italic(I)^{2} == ", round(I2_total_all[1], 2), "*\"%\""),
    color = "black",
    parse = TRUE,
    size = 3
  ) +
  ggtitle("A) All Outcome Categories") +
  theme(plot.title = element_text(face = "bold")) +
  scale_fill_manual(values = "grey") +
  scale_colour_manual(values = "grey")

#########Feed behaviour
#Feed behaviour
# Filter for feed behaviour
feed_data <- mydata %>%
  filter(outcome_category == "feed behaviour")
feed_data$intervention_preparation <- as.factor(feed_data$intervention_preparation)
print(feed_data)

# Calculate effect size (Hedges g) and pooled SD
feed_data <- feed_data %>%
  mutate(
    s_pooled_feed = sqrt(((treatment_n - 1) * treatment_SD^2 + (control_n - 1) * control_SD^2) /
                    (treatment_n + control_n - 2)),
    J_feed = 1 - (3 / (4 * (treatment_n + control_n) - 9)),
    g_feed = J_feed * ((treatment_mean - control_mean) / s_pooled_feed),
    vi_g_feed = (treatment_n + control_n) / (treatment_n * control_n) + (g_feed^2 / (2 * (treatment_n + control_n - 2)))
  )

#Run MLMA with species as additional random effect
res_3L_feed_species <- rma.mv(yi = g_feed, V = vi_g_feed,
  random = list(
    ~ 1 | species,                  
    ~ 1 | study_ID / experiment_ID),
data = feed_data,
method = "REML")
res_3L_feed_species

I2_feed <- orchaRd::i2_ml(res_3L_feed_species)

feed_data_plot <- orchaRd::orchard_plot(
  res_3L_feed_species,
  mod = "1",
  xlab = "Effect Size (Hedges g)",
  group = "study_ID"
) +
  annotate(
    geom = "text",
    x = 0.7,
    y = 10,
    label = paste0("italic(I)^{2} == ", round(I2_feed[1], 2), "*\"%\""),
    color = "black",
    parse = TRUE,
    size = 3
  ) +
  scale_fill_manual(values = "pink") +
  scale_colour_manual(values = "pink") +
  ggtitle("B) Feed Behaviour") +
  theme(plot.title = element_text(face = "bold"))

######### Growth Performance
# Filter for growth performance
growth_data <- mydata %>%
  filter(outcome_category == "growth performance")
growth_data$intervention_preparation <- as.factor(growth_data$intervention_preparation)
print(growth_data)

#Summarise number of studies and effect sizes
growth_dataset <- growth_data %>%
  summarise(
    n_studies = n_distinct(study_ID),
    n_effects = n_distinct(experiment_ID),
    .groups = "drop"
  )
growth_dataset

# Calculate effect size (Hedges g) and pooled SD
growth_data <- growth_data %>%
  mutate(
    s_pooled_growth = sqrt(((treatment_n - 1) * treatment_SD^2 + (control_n - 1) * control_SD^2) /
                    (treatment_n + control_n - 2)),
    J_growth = 1 - (3 / (4 * (treatment_n + control_n) - 9)),
    g_growth = J_growth * ((treatment_mean - control_mean) / s_pooled_growth),
    vi_g_growth = (treatment_n + control_n) / (treatment_n * control_n) + (g_growth^2 / (2 * (treatment_n + control_n - 2)))
  )

# Run MLMA with species as additional random effect
res_3L_growth_species <- rma.mv(yi = g_growth, V = vi_g_growth,
  random = list(
    ~ 1 | species,                  
    ~ 1 | study_ID / experiment_ID
  ),
  data = growth_data,
  method = "REML"
)
res_3L_growth_species

I2_growth <- orchaRd::i2_ml(res_3L_growth_species)

growth_data_plot <- orchaRd::orchard_plot(
   res_3L_growth_species,
    mod = "1",
    xlab = "Effect Size (Hedges g)",
    group = "study_ID"
) +
  annotate(
    geom = "text",
    x = 0.7,
    y = 7,
    label = paste0("italic(I)^{2} == ", round(I2_growth[1], 2), "*\"%\""),
    color = "black",
    parse = TRUE,
    size = 3
  ) +
  scale_fill_manual(values = "blue") +
  scale_colour_manual(values = "blue") + ggtitle("C) Growth Performance") + theme(plot.title = element_text(face = "bold"))

######### Nutrient Utilisation
# Filter for nutrient utilisation
nutrient_data <- mydata %>%
  filter(outcome_category == "nutrient utilisation")
nutrient_data$intervention_preparation <- as.factor(nutrient_data$intervention_preparation)
print(nutrient_data)

# Summarise number of studies and effect sizes
nutrient_dataset <- nutrient_data %>%
  summarise(
    n_studies = n_distinct(study_ID),
    n_effects = n_distinct(experiment_ID),
    .groups = "drop"
  )
nutrient_dataset

# Calculate effect size (Hedges g) and pooled SD
nutrient_data <- nutrient_data %>%
  mutate(
    s_pooled_nutrient = sqrt(((treatment_n - 1) * treatment_SD^2 + (control_n - 1) * control_SD^2) /
                    (treatment_n + control_n - 2)),
    J_nutrient = 1 - (3 / (4 * (treatment_n + control_n) - 9)),
    g_nutrient = J_nutrient * ((treatment_mean - control_mean) / s_pooled_nutrient),
    vi_g_nutrient = (treatment_n + control_n) / (treatment_n * control_n) + (g_nutrient^2 / (2 * (treatment_n + control_n - 2)))
  )

# Reverse sign for outcomes where lower values are better (e.g., FCR)
nutrient_data <- nutrient_data %>%
  mutate(
    g_nutrient = ifelse(outcome == "FCR", -g_nutrient, g_nutrient),
    vi_g_nutrient = ifelse(outcome == "FCR", vi_g_nutrient, vi_g_nutrient))

# Run MLMA with species as additional random effect
res_3L_nutrient_species <- rma.mv(yi = g_nutrient, V = vi_g_nutrient,
  random = list(
    ~ 1 | species,                  
    ~ 1 | study_ID / experiment_ID
  ),
  data = nutrient_data,
  method = "REML"
)
res_3L_nutrient_species

I2_nutrient <- orchaRd::i2_ml(res_3L_nutrient_species)

nutrient_data_plot <- orchaRd::orchard_plot(
   res_3L_nutrient_species,
    mod = "1",
    xlab = "Effect Size (Hedges g)",
    group = "study_ID"
) +
  annotate(
    geom = "text",
    x = 0.7,
    y = 12,
    label = paste0("italic(I)^{2} == ", round(I2_nutrient[1], 2), "*\"%\""),
    color = "black",
    parse = TRUE,
    size = 3
  ) +
  scale_fill_manual(values = "green") +
  scale_colour_manual(values = "green")  +
  ggtitle("D) Nutrient Utilisation") + theme(
  plot.title = element_text(face = "bold")
)

(all_data_plot | feed_data_plot) /
(growth_data_plot | nutrient_data_plot)
