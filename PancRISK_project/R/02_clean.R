# Load libraries ----------------------------------------------------------
library("tidyverse")
library("stringr")


# Load data ---------------------------------------------------------------
my_data <- read_tsv(file = "PancRISK_project/data/01_my_data.tsv")


# Wrangle data ------------------------------------------------------------
my_data_clean <- my_data %>% 
  as_tibble %>% 
  bind_cols %>% 
  mutate(diagnosis = case_when(diagnosis == 1 ~ "control",
                               diagnosis == 2 ~ "benign",
                               diagnosis == 3 ~ "malignant"),
         sex = case_when(sex == "F" ~ "Female",
                         sex == "M" ~ "Male"),
         sampling_country = case_when(sample_origin == "BPTB" ~ "United Kingdom",
                                      sample_origin == "UCL" ~ "United Kingdom",
                                      sample_origin == "LIV" ~ "United Kingdom",
                                      sample_origin == "ESP" ~ "Spain"),
         sample_id = gsub('^.', '', sample_id)) %>% 
  select(.,
         -sample_origin,
         -benign_sample_diagnosis,
         -stage) %>% 
  drop_na(plasma_CA19_9)

# Write data --------------------------------------------------------------
write_tsv(x = my_data_clean,
          file = "PancRISK_project/data/02_my_data_clean.tsv")