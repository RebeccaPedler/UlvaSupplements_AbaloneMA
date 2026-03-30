# Project: Feeding behaviour, growth performance and nutrient utilisation of abalone with dietary Ulva sp. supplementation: A meta-analysis 

## Step1: Example script for deduplicating references downloaded from Scopus and Web of Science (WoS) sources 
### Note: Deduplication is based on Title field, after standardizing case and removing punctuation and extra spaces

## Install.packages(c("tidyverse", "here"))

## Load required libraries
library(tidyverse)
library(here)

## Load WoS and Scopus data
df <- read.csv(here("Searches/primary_literature", "WOS_and_SCOPUS_combined_11092025.csv"), header = TRUE, stringsAsFactors = FALSE)
dim(df) #997 rows

names(df)

## Clean the titles  remove duplicates
df_unique <- df %>%
  mutate(
    Title_clean = Title %>%
      stringr::str_to_lower() %>%                     
      stringr::str_replace_all("[[:punct:]]", "") %>% 
      stringr::str_replace_all("\\s+", " ") %>%       
      stringr::str_trim()                             
  ) %>%
  distinct(Title_clean, .keep_all = TRUE)

dim(df_unique) #718 rows

## Report how many duplicates were removed
cat("Original rows:", dim(df)[1], "\n")
cat("Rows after removing duplicates:", dim(df_unique)[1], "\n")
cat("Duplicates removed:", dim(df)[1] - dim(df_unique)[1], "\n")

##   Write deduplicated records to CSV file
write.csv(df_unique, here("Searches/primary_literature", "WOS_and_SCOPUS_combined_11092025_cleaned.csv"), row.names = FALSE)
