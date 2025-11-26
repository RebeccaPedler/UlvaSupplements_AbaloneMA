install.packages(c("metafor","meta","tidyr", "tidyverse", "rnaturalearth", "rnaturalearthdata", "sf", "ggplot2", "grid", "dplyr", "patchwork"))

# Load necessary libraries
library(dplyr)
library(metafor)
library(meta)
library(tidyr)
library(tidyverse)
library(rnaturalearth)   # for map data
library(rnaturalearthdata)
library(sf)              # for spatial handling
library(ggplot2)
library(grid)
library(patchwork)

# Read the CSV file and check structure
mydata <- read.csv("C:/Users/RebeccaPedler/OneDrive - Yumbah/Documents/R&D/Industry PhD/Trials/Meta analysis/Ulva inclusion in Haliotis sp. diets_ A Meta-analysis.CSV")
str(mydata)

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
      outcome == "Ab" ~ "Absorbed feed energy", 
      outcome == "S" ~ "Shell growth energy", 
      outcome == "Pg" ~ "Somatic growth energy",
      TRUE ~ outcome
    )
  )

data <- dplyr::left_join(mydata, citations, by = "study_ID")
head(data)

#####Study Characteristics

# Get unique study IDs
unique_ids <- unique(mydata$study_ID)
print(unique_ids)
total_studies <- length(unique(data$study_ID))
print(total_studies)

# Step 1: Summarise number of records per species–outcome combination
bubble_data <- mydata %>%
  group_by(species, outcome, outcome_category) %>%
  summarise(
    count = n(),        # number of times each combo appears
    .groups = "drop"
  )

# Wrap long species names onto two lines
bubble_plot_species <- ggplot(bubble_data, aes(x = species, y = outcome)) +
  geom_point(
    aes(size = count, fill = outcome_category),
    shape = 21,
    color = "black",
    alpha = 0.9
  ) +
  scale_size_continuous(
    range = c(3, 15),
    breaks = c(2, 5, 10),
    name = "Number of Observations"
  ) +
  scale_fill_manual(
    values = c(
      "growth performance" = "steelblue",
      "nutrient utilisation" = "grey70",
      "feed behaviour" = "white"
    )
  ) +
  labs(
    title = "",
    x = "Species",
    y = "Outcome",
    size = "Number of Observations",
    fill = "Outcome Category"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(hjust = 0.5),  # centered text
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.05),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.25, "cm"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks.outside = TRUE,
    panel.border = element_blank(),
    legend.position = "right"
  ) +
  coord_cartesian(clip = "off")

bubble_plot_species <- bubble_plot_species +
  scale_x_discrete(labels = function(x) {
    sapply(x, function(label) {
      # check if label contains " x " (hybrid)
      if (grepl(" x ", label)) {
        parts <- strsplit(label, " ")[[1]]  # split by space
        # join first two words as first line, "x" as second line, last two words as third line
        paste0(parts[1], " ", parts[2], "\n", parts[3], "\n", parts[4], " ", parts[5])
      } else {
        # normal species: keep genus + species on one line
        paste(label)
      }
    })
  })

bubble_plot_species <- bubble_plot_species +
  guides(
    size = guide_legend(
      override.aes = list(size = 8)   # increase legend bubble size
    ),
    fill = guide_legend(
      override.aes = list(size = 8)   # optional: also affect fill legend if needed
    )
  )
# Bubble plot for species with multi-line labels and larger legend bubbles
bubble_plot_species <- ggplot(bubble_data, aes(x = species, y = outcome)) +
  geom_point(
    aes(size = count, fill = outcome_category),
    shape = 21,         # filled circle with border
    color = "black",
    alpha = 0.9
  ) +
  scale_size_continuous(
    range = c(3, 15),
    breaks = c(2, 5, 10),
    name = "Number of Observations"
  ) +
  scale_fill_manual(
    values = c(
      "growth performance" = "steelblue",
      "nutrient utilisation" = "grey70",
      "feed behaviour" = "white"
    )
  ) +
  labs(
    title = "A",
    x = "Species",
    y = "Outcome",
    size = "Number of Observations",
    fill = "Outcome Category"
  ) +
  scale_x_discrete(labels = function(x) {
    sapply(x, function(label) {
      if (grepl(" x ", label)) {
        # split hybrid label into three lines
        parts <- strsplit(label, " ")[[1]]
        paste0(parts[1], " ", parts[2], "\n", parts[3], "\n", parts[4], " ", parts[5])
      } else {
        # keep normal species on one line
        label
      }
    })
  }) +
  guides(
    size = guide_legend(override.aes = list(size = 8)),  # larger legend bubbles
    fill = guide_legend(override.aes = list(size = 8))   # optional for fill legend
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(hjust = 0.5, face = "italic"),          # center multi-line labels
    axis.text.y = element_text(size = 10),
    plot.title = element_text(face = "bold", hjust = 0.05),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank(),
    axis.ticks.length = unit(0.25, "cm"),
    axis.ticks = element_line(color = "black"),
    axis.line = element_line(color = "black"),
    axis.ticks.outside = TRUE,
    panel.border = element_blank(),
    legend.position = "right"
  ) +
  coord_cartesian(clip = "off")
bubble_plot_species

# Recode country names to match rnaturalearth conventions
mydata <- mydata %>%
  mutate(country_standardised = case_when(
    country == "Australia"     ~ "Australia",
    country == "China"         ~ "China",
    country == "Hawaii"        ~ "United States of America",  # Hawaii is part of USA
    country == "Indonesia"     ~ "Indonesia",
    country == "Ireland"       ~ "Ireland",
    country == "Japan"         ~ "Japan",
    country == "Korea"         ~ "South Korea",               # or "Republic of Korea"
    country == "New Zealand"   ~ "New Zealand",
    country == "South Africa"  ~ "South Africa",
    country == "Taiwan"        ~ "Taiwan",                    # not always present in shapefile
    TRUE                       ~ country                      # default: keep original
  ))

country_summary <- mydata %>%
  group_by(country_standardised) %>%
  summarise(
    observation_count = n(),                                  # count all rows
    percentage = 100 * observation_count / nrow(mydata),      # proportion of total rows
    .groups = "drop"
  ) %>%
  rename(name = country_standardised)
print(country_summary)

# Get world map shapefile
world <- ne_countries(scale = "medium", returnclass = "sf")

# Join country summary with world shapefile
world_data <- world %>%
  left_join(country_summary, by = "name")

# Step 4: Plot heat map
world_map <- ggplot(world_data) +
  geom_sf(aes(fill = observation_count)) +
  scale_fill_gradient(
    low = "#deebf7",
    high = "#08306b",
    na.value = "grey90",
    name = "Number of\nObservations"   # line break in legend title
  ) +
  labs(
    title = "B",
    caption = "Grey areas = no studies"
  ) +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(face = "italic"),
    panel.grid = element_blank()           # remove all gridlines
  )

# Combine side by side with one shared legend
species_and_map_combined <- bubble_plot_species + world_map +
  plot_annotation(tag_levels = 'A') &          # & applies theme to all plots
  theme(plot.tag = element_text(face = "bold", size = 14))
species_and_map_combined

# List validity metrics
validity_metrics <- c(
  "random_assignment",
  "blind_measurement",
  "protocol_registration",
  "raw_data_availability",
  "source_code_availability",
  "funding_disclosed",
  "funding_potential_COI",
  "COI_disclosed"
)

# Loop through each metric and summarise
validity_summary <- data %>%
  select(study_ID, all_of(validity_metrics)) %>%
  distinct() %>%
  pivot_longer(
    cols = all_of(validity_metrics),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(metric, value) %>%
  summarise(
    study_count = n(),
    percentage = 100 * study_count / total_studies,
    .groups = "drop"
  ) %>%
  arrange(metric, desc(study_count))

print(validity_summary)

# Summarise number and percentage of unique studies per publication year
publication_summary <- data %>%
  distinct(study_ID, publication_year) %>%
  group_by(publication_year) %>%
  summarise(unique_study_count = n()) %>%
  mutate(
    percentage = 100 * unique_study_count / total_studies,
    label_percent = paste0(round(percentage, 1), "%"),
    label_count = unique_study_count
  )
print(publication_summary)


#########Feed behaviour

#Feed behaviour
# Filter for feed behaviour
feed_data <- data %>%
  filter(outcome_category == "feed behaviour", study_ID != "S004")
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

# Run meta-analysis
res_HG_feed_behaviour <- rma(yi = g_feed, vi = vi_g_feed, data = feed_data, method = "REML")
summary(res_HG_feed_behaviour)
forest(
  res_HG_feed_behaviour,
  slab = feed_data$short_citation,
  xlab = "Hedges g",
  main = "Feed Behaviour - All Data"
)

#3 level MA
res_3L_feed <- rma.mv(yi = g_feed, V = vi_g_feed,
                 random = ~ 1 | study_ID / experiment_ID,
                 data = feed_data,
                 method = "REML")
res_3L_feed

# Extract variance components
sigma_study <- res_3L_feed$sigma2[1]   # between-study
sigma_exp   <- res_3L_feed$sigma2[2]   # within-study / experimental

# Mean sampling variance
mean_vi <- mean(feed_data$vi_g_feed)

# Calculate I²
I2_study <- sigma_study / (sigma_study + sigma_exp + mean_vi) * 100
I2_exp   <- sigma_exp / (sigma_study + sigma_exp + mean_vi) * 100
I2_total <- (sigma_study + sigma_exp) / (sigma_study + sigma_exp + mean_vi) * 100
I2_study; I2_exp; I2_total

#Adding intervention_dose as a moderator
#Suspected quadratic relationship
res_3L_feed_dose_quad <- rma.mv(
  yi = g_feed,
  V  = vi_g_feed,
  mods = ~ intervention_dose + I(intervention_dose^2),
  random = ~ 1 | study_ID/experiment_ID,
  data = feed_data,
  method = "REML"
)
summary(res_3L_feed_dose_quad)

# variance components from model with and without moderator
sigma_null <- res_3L_feed$sigma2
sigma_mod  <- res_3L_feed_dose_quad$sigma2

# pseudo-R² for each level
R2_study <- (sigma_null[1] - sigma_mod[1]) / sigma_null[1]
R2_exp   <- (sigma_null[2] - sigma_mod[2]) / sigma_null[2]

R2_study; R2_exp

# Create a sequence of doses for prediction
dose_seq <- seq(min(feed_data$intervention_dose), max(feed_data$intervention_dose), length.out = 100)

# Predicted effect sizes using your quadratic model
pred_g <- res_3L_feed_dose_quad$beta[1] +
          res_3L_feed_dose_quad$beta[2] * dose_seq +
          res_3L_feed_dose_quad$beta[3] * dose_seq^2

# Approximate standard errors for the predictions (delta method)
# Here we use just fixed effects; for full CI accounting for random effects you can use predict()
pred_se <- sqrt(res_3L_feed_dose_quad$vb[1,1] +
                dose_seq^2 * res_3L_feed_dose_quad$vb[2,2] +
                (dose_seq^2)^2 * res_3L_feed_dose_quad$vb[3,3] +
                2 * dose_seq * res_3L_feed_dose_quad$vb[1,2] +
                2 * dose_seq^2 * res_3L_feed_dose_quad$vb[1,3] +
                2 * dose_seq^3 * res_3L_feed_dose_quad$vb[2,3])

# Confidence intervals
ci_lb <- pred_g - 1.96 * pred_se
ci_ub <- pred_g + 1.96 * pred_se

# Optimum dose (from quadratic formula)
opt_dose <- -res_3L_feed_dose_quad$beta[2] / (2 * res_3L_feed_dose_quad$beta[3])

# Create prediction dataframe
df_plot <- data.frame(dose = dose_seq, pred = pred_g, ci.lb = ci_lb, ci.ub = ci_ub)

# Plot
ggplot(df_plot, aes(x = dose, y = pred)) +
  geom_line(color = "steelblue", size = 1) +
  geom_ribbon(aes(ymin = ci.lb, ymax = ci.ub), fill = "steelblue", alpha = 0.2) +
  geom_point(data = feed_data, aes(x = intervention_dose, y = g_feed),
             color = "black", fill = "white", shape = 21, size = 3) +
  geom_vline(xintercept = opt_dose, linetype = "dashed", color = "red", size = 1) +
  annotate("text", x = opt_dose, y = max(pred_g)+0.2, 
           label = paste0("Optimum ~", round(opt_dose,1), "%"), 
           color = "red", hjust = 0.5) +
  labs(
    title = "Optimum Ulva Inclusion for Growth Performance",
    x = "Ulva Inclusion (%)",
    y = "Predicted Hedges' g (Effect Size)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid = element_blank(),  # remove all internal gridlines
    axis.line = element_line(color = "black", size = 0.8), # black border lines
    axis.ticks.length = unit(0.3, "cm"),
    axis.ticks = element_line(color = "black"),
    axis.title.x = element_text(margin = margin(t = 10)),
    axis.title.y = element_text(margin = margin(r = 10)),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )
