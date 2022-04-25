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
         patient_cohort = case_when(patient_cohort == "Cohort1" ~ "Cohort_1",
                                    patient_cohort == "Cohort2" ~ "Cohort_2"),
         sampling_country = case_when(sample_origin == "BPTB" ~ "United Kingdom",
                                      sample_origin == "UCL" ~ "United Kingdom",
                                      sample_origin == "LIV" ~ "United Kingdom",
                                      sample_origin == "ESP" ~ "Spain"),
         sample_id = gsub('^.', '', sample_id)) %>% 
  select(.,
         -sample_origin,
         -benign_sample_diagnosis)

# Write data --------------------------------------------------------------
write_tsv(x = my_data_clean,
          file = "PancRISK_project/data/02_my_data_clean.tsv")