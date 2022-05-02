## Modelling data augment

my_data_clean <- read_tsv(file = "PancRISK_project/data/02_my_data_clean.tsv")

my_data_cleaner <-
  my_data_clean %>% 
  mutate(diagnosis = factor(diagnosis)) %>%
  mutate(sex = factor(sex)) %>%
  mutate(diagnosis = relevel(diagnosis, ref = "control"))

my_data_REG1A_NA <-
  my_data_cleaner %>% 
  filter(is.na(REG1A))

my_data_REG1A_zero <-
  my_data_cleaner %>% 
  filter(REG1A == 0 & !is.na(REG1A))

my_data_cleaner <-
  my_data_cleaner %>% 
  filter(!is.na(plasma_CA19_9) & !is.na(REG1A) & REG1A != 0 & plasma_CA19_9 != 0) %>%
  mutate_at(c("REG1A", "REG1B", "LYVE1", "TFF1"), log) %>%
  full_join(my_data_REG1A_zero) %>%
  mutate_at(c("REG1A", "REG1B", "LYVE1", "TFF1"), scale, scale = FALSE) %>% 
  full_join(my_data_REG1A_NA) %>%
  mutate(REG1A = REG1A[,1],
         REG1B = REG1B[,1],
         LYVE1 = LYVE1[,1],
         TFF1 = TFF1[,1]) %>%
  mutate(cutoff_plasma = ifelse(plasma_CA19_9 > 37, 1, 0)) %>% 
  group_by(diagnosis)

write_tsv(x = my_data_cleaner,
          file = "PancRISK_project/data/03_my_data_cleaner.tsv")




## model_ii augmented data 

## PCA data

my_data_clean_PCA <- select(my_data_clean , patient_cohort,age,sex, diagnosis ,stage , plasma_CA19_9,creatinine , LYVE1, REG1B, TFF1) %>%
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
    stage = replace_na(stage,0)) %>%
  drop_na(REG1B, LYVE1, TFF1, plasma_CA19_9)


## Augmenting the clean data with the data fitted on each PC (will be used in for the K-means )

data_aug_pca <- pca_fit %>% 
  augment(my_data_clean_PCA)


# This model includes all the data features 

## K-means : 

## K-means on the whole data : 

## Full model including plasma_CA19_9 + creatinine + LYVE1 + REG1B + TFF1 

data_aug_pca_aug_k_means <- data_aug_pca %>%
  select(patient_cohort,age,sex, plasma_CA19_9,creatinine , LYVE1, REG1B, TFF1) %>%
  kmeans(centers = 3)%>%
  augment(data_aug_pca) %>% 
  rename(clusters_org = .cluster)%>%
  mutate(clusters_org = case_when(
    clusters_org == 1 ~ "1",
    clusters_org == 2 ~ "2",
    clusters_org == 3 ~ "3"))

## K-means on the PC1 and PC2

k_means_onto_PC_1_2 <- data_aug_pca_aug_k_means %>%
  select(.fittedPC1, .fittedPC2)  %>%
  kmeans(centers = 3) %>%
  augment(data_aug_pca_aug_k_means) %>% 
  rename(clusters_pca = .cluster)%>%
  mutate(clusters_pca = case_when(
    clusters_pca == 1 ~ "1",
    clusters_pca == 2 ~ "2",
    clusters_pca == 3 ~ "3"),
    diagnosis = case_when(
      diagnosis == 0 ~ "control",
      diagnosis == 1 ~ "benign",
      diagnosis == 2 ~ "malignant"))

## Model including only plasma_CA19_9 

plasma_CA19_9_alone_k_means <- data_aug_pca %>%
  select(patient_cohort,age,sex,plasma_CA19_9 ) %>%
  kmeans(centers = 3)%>%
  augment(data_aug_pca) %>% 
  rename(clusters_org = .cluster)

## Model including only plasma_CA19_9 

biomarkers_only_k_means <- data_aug_pca %>%
  select(patient_cohort,age,sex, LYVE1, REG1B, TFF1)  %>%
  kmeans(centers = 3)%>%
  augment(data_aug_pca) %>% 
  rename(clusters_org = .cluster)

## Full model including plasma_CA19_9 + creatinine + LYVE1 + REG1B + TFF1 

complete_model_k_means <- data_aug_pca %>%
  select(patient_cohort,age,sex, LYVE1, REG1B, TFF1,plasma_CA19_9)  %>%
  kmeans(centers = 3)%>%
  augment(data_aug_pca) %>% 
  rename(clusters_org = .cluster)


