## Modelling data augment

my_data_clean <- read_tsv(file = "PancRISK_project/data/02_my_data_clean.tsv")

my_data_cleaner <-
  my_data_clean %>% 
  mutate(diagnosis = factor(diagnosis)) %>%
  mutate(sex = factor(sex)) %>%
  mutate(diagnosis = relevel(diagnosis, ref = "control"))

my_data_REG1A_NA <-
  my_data_cleaner %>% 
  filter(is.na(REG1A)) %>% 
  mutate_at(c("REG1B", "LYVE1", "TFF1"), log)

my_data_REG1A_zero <-
  my_data_cleaner %>% 
  filter(REG1A == 0 & !is.na(REG1A)) %>% 
  mutate_at(c("REG1B", "LYVE1", "TFF1"), log)

my_data_cleaner <-
  my_data_cleaner %>% 
  filter(!is.na(plasma_CA19_9) & !is.na(REG1A) & REG1A != 0 & plasma_CA19_9 != 0) %>%
  mutate_at(c("REG1A", "REG1B", "LYVE1", "TFF1"), log) %>%
  full_join(my_data_REG1A_zero) %>%
  mutate_at(c("REG1A"), scale, scale = FALSE, center = TRUE) %>% 
  full_join(my_data_REG1A_NA) %>%
  mutate_at(c("REG1B", "LYVE1", "TFF1"), scale, scale = FALSE, center = TRUE) %>%
  mutate(REG1A = REG1A[,1],
         REG1B = REG1B[,1],
         LYVE1 = LYVE1[,1],
         TFF1 = TFF1[,1]) %>%
  mutate(cutoff_plasma = ifelse(plasma_CA19_9 > 37, 1, 0)) %>% 
  mutate(diagnosis = as.factor(diagnosis)) %>%
  group_by(diagnosis)

write_tsv(x = my_data_cleaner,
          file = "PancRISK_project/data/03_my_data_cleaner.tsv")




## model_ii augmented data 

## PCA data

my_data_clean_PCA <- select(my_data_clean, patient_cohort, age, sex, 
                            diagnosis, stage, plasma_CA19_9, creatinine, 
                            LYVE1, REG1B, TFF1) %>%
  mutate(patient_cohort = case_when(
    patient_cohort == "Cohort_1" ~ 0,
    patient_cohort == "Cohort_2" ~ 1),
    diagnosis = case_when(
    diagnosis == "control" ~ 0,
    diagnosis == "benign" ~ 1,
    diagnosis == "malignant" ~ 2),
    sex = case_when(
      sex == "Male" ~ 0,
      sex == "Female" ~ 1),
    stage = case_when(
      stage  == "IA" ~ 1,
      stage  == "IB" ~ 1,
      stage  == "IIA" ~ 2,
      stage  == "IIB" ~ 2,
      stage  == "III" ~ 3,
      stage  == "IV" ~ 4),
    stage = replace_na(stage, 0)) %>%
  drop_na(REG1B, LYVE1, TFF1, plasma_CA19_9)

write_tsv(x = my_data_clean_PCA,
          file = "PancRISK_project/data/04_my_data_clean_PCA.tsv")

