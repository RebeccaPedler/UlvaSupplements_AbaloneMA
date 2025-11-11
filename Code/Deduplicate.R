install.packages(c("stringr","janitor", "tidyverse", "here"))

# Load required libraries
library(tidyverse)
library(janitor)
library(stringr)
library(here)

###Please download GitHub repository and then run the following
here()
df <- read_csv(here("GitHub", "UlvaSupplements_AbaloneMA","Searches", "primary_literature", "WOS AND SCOPUS 05062025.csv"))
head(data)

# Define the file name
df <- read.csv("WOS and SCOPUS 11092025.csv", header = TRUE, stringsAsFactors = FALSE)
summary(df)

# Step 1: Count original rows
original_n <- nrow(df)
print(original_n)

# Step 2: Clean and remove duplicates
df_unique <- df %>%
  mutate(
    Title_clean = Title %>%
      str_to_lower() %>%                     
      str_replace_all("[[:punct:]]", "") %>% 
      str_replace_all("\\s+", " ") %>%       
      str_trim()                             
  ) %>%
  distinct(Title_clean, .keep_all = TRUE) %>%
  select(-Title_clean)

# Step 3: Count cleaned rows
cleaned_n <- nrow(df_unique)

# Step 4: Report how many duplicates were removed
cat("Original rows:", original_n, "\n")
cat("Rows after removing duplicates:", cleaned_n, "\n")
cat("Duplicates removed:", original_n - cleaned_n, "\n")

# Step 5: Write to CSV
# Set working directory
setwd("C:/Users/RebeccaPedler/OneDrive - Yumbah")
write.csv(df_unique, "cleaned_references.csv", row.names = FALSE)

> cat("Original rows:", original_n, "\n")
Original rows: 3068 
> cat("Rows after removing duplicates:", cleaned_n, "\n")
Rows after removing duplicates: 2132 
> cat("Duplicates removed:", original_n - cleaned_n, "\n")
Duplicates removed: 936 
> # Step 5: Write to CSV
> write.csv(df_unique, "cleaned_references.csv", row.names = FALSE)
>


