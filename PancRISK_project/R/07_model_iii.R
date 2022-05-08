data_aug_pca <- pca_fit %>% 
  augment(my_data_clean_PCA) %>%
  filter(.fittedPC1> -6) # removing an outlier for better K-means 


malignant_data_aug_pca <- malignant_PCA_fit %>% 
  augment(malignant_PCA)


# Model data -------------------------------------------------------------------

Kmeans_blood_diagnosis_ <- data_aug_pca %>%
  select(patient_cohort,
         age,
         diagnosis,
         sex,
         plasma_CA19_9) %>%
  kmeans(centers = 3) %>%
  augment(data_aug_pca) %>% 
  rename(clusters = .cluster) %>% 
  mutate(diagnosis = case_when(diagnosis == 0 ~ "control",
                               diagnosis == 1 ~ "benign",
                               diagnosis == 2 ~ "malignant"))


Kmeans_blood_diagnosis_control <- Kmeans_blood_diagnosis_ %>% 
  filter(diagnosis == "control") %>%
  filter(.fittedPC1 > -6) 


Kmeans_blood_diagnosis_benign <- Kmeans_blood_diagnosis_ %>% 
  filter(diagnosis == "benign") %>%
  filter(.fittedPC1 > -6) 

Kmeans_blood_diagnosis_malignant <- Kmeans_blood_diagnosis_ %>% 
  filter(diagnosis == "malignant") %>%
  filter(.fittedPC1 > -6) 



blood_control_statistics <- Kmeans_blood_diagnosis_control %>%
  group_by(clusters) %>%
  summarise(count = n())

blood_benign_statistics <- Kmeans_blood_diagnosis_benign %>%
  group_by(clusters) %>%
  summarise(count = n())

blood_malignant_statistics <- Kmeans_blood_diagnosis_malignant %>%
  group_by(clusters) %>%
  summarise(count = n())



Kmeans_urinary_diagnosis_ <- data_aug_pca %>%
  select(patient_cohort,
         age,
         diagnosis,
         sex, 
         LYVE1, 
         REG1B, 
         TFF1)  %>%
  kmeans(centers = 3) %>%
  augment(data_aug_pca) %>% 
  rename(clusters = .cluster) %>% 
  mutate(diagnosis = case_when(diagnosis == 0 ~ "control",
                               diagnosis == 1 ~ "benign",
                               diagnosis == 2 ~ "malignant"))


Kmeans_urinary_diagnosis_control <- Kmeans_urinary_diagnosis_ %>% 
  filter(diagnosis == "control")

Kmeans_urinary_diagnosis_benign <- Kmeans_urinary_diagnosis_ %>% 
  filter(diagnosis == "benign")

Kmeans_urinary_diagnosis_malignant <- Kmeans_urinary_diagnosis_ %>% 
  filter(diagnosis == "malignant")

control_statistics <- Kmeans_urinary_diagnosis_control %>%
  group_by(clusters) %>%
  summarise(count = n())

benign_statistics <- Kmeans_urinary_diagnosis_benign %>%
  group_by(clusters) %>%
  summarise(count = n())

malignant_statistics <- Kmeans_urinary_diagnosis_malignant %>%
  group_by(clusters) %>%
  summarise(count = n())


Kmeans_complete_diagnosis_ <- data_aug_pca %>%
  select(patient_cohort,
         age,
         diagnosis,
         sex,
         LYVE1, 
         REG1B, 
         TFF1,
         plasma_CA19_9) %>%
  kmeans(centers = 3) %>%
  augment(data_aug_pca) %>% 
  rename(clusters = .cluster) %>% 
  mutate(diagnosis = case_when(diagnosis == 0 ~ "control",
                               diagnosis == 1 ~ "benign",
                               diagnosis == 2 ~ "malignant"))

Kmeans_complete_diagnosis_control <- Kmeans_complete_diagnosis_ %>% 
  filter(diagnosis == "control")

Kmeans_complete_diagnosis_benign <- Kmeans_complete_diagnosis_ %>% 
  filter(diagnosis == "benign")

Kmeans_complete_diagnosis_malignant <- Kmeans_complete_diagnosis_ %>% 
  filter(diagnosis == "malignant")

complete_control_statistics <- Kmeans_complete_diagnosis_control %>%
  group_by(clusters) %>%
  summarise(count = n())

complete_benign_statistics <- Kmeans_complete_diagnosis_benign %>%
  group_by(clusters) %>%
  summarise(count = n())

complete_malignant_statistics <- Kmeans_complete_diagnosis_malignant %>%
  group_by(clusters) %>%
  summarise(count = n())



# Visualise data ---------------------------------------------------------------


# How plasma CA19_9 (blood diagnosis ) can discriminate between control benign and malignant individuals


blood_diagnosis_plot <- ggplot() +
  geom_point(data = Kmeans_blood_diagnosis_control,
             mapping = aes(x = .fittedPC1,
                           y = .fittedPC2 ,
                           color = clusters, 
                           shape = diagnosis), 
                           size = 2) +
  geom_point(data = Kmeans_blood_diagnosis_benign, 
             mapping = aes(x = .fittedPC1, 
                           y = .fittedPC2 , 
                           color = clusters, 
                           shape = diagnosis),
                           size = 2) +
  geom_point(data = Kmeans_blood_diagnosis_malignant, 
             mapping = aes(x = .fittedPC1, 
                           y = .fittedPC2, 
                           color = clusters, 
                           shape = diagnosis), 
                           size = 2) +
  scale_shape_manual(values =c(1, 2, 3)) +
  theme(plot.title = element_text(size = 10)) +
  theme_minimal() +
  background_grid() +
  labs(title = "plasma_CA19_9",
       x = "PC1",
       y = "PC2") +
  xlim(-6, 3.5) +
  ylim(-5, 5)


blood_diagnosis_plot_ <- blood_diagnosis_plot +
  guides(size = 'none') +
  theme(legend.position = "bottom",
        legend.box = "vertical", 
        legend.margin = margin())



# How the 3 biomarkers (urinary diagnosis ) can discriminate between control benign and malignant individuals


urinary_diagnosis_plot <- ggplot() +
  geom_point(data = Kmeans_urinary_diagnosis_control, 
             mapping = aes(x = .fittedPC1, 
                           y = .fittedPC2, 
                          color = clusters, 
                          shape = diagnosis), 
                          size = 2) +
  geom_point(data= Kmeans_urinary_diagnosis_benign, 
             mapping = aes(x = .fittedPC1, 
                           y = .fittedPC2, 
                           color = clusters,
                           shape = diagnosis), 
                           size = 2) +
  geom_point(data = Kmeans_urinary_diagnosis_malignant, 
             mapping = aes(x = .fittedPC1, 
                           y = .fittedPC2 , 
                           color = clusters , 
                           shape = diagnosis), 
                           size = 2) +
  scale_shape_manual(values = c(1 ,2 ,3)) +
  theme(plot.title = element_text(size = 10)) +
  theme_minimal() +
  background_grid() +
  labs(title = "LYVE1 + REG1B + TFF1",
       x = "PC1",
       y = "PC2") +
  xlim(-6, 3.5) +
  ylim(-5, 5)

urinary_diagnosis_plot_ <- urinary_diagnosis_plot +
  guides(size = 'none') +
  theme(legend.box = "vertical", 
        legend.margin = margin(),
        legend.position = "bottom")


kmeans_plot_ <- (blood_diagnosis_plot_ | urinary_diagnosis_plot_) +
  plot_annotation(title = "Comparison between Blood vs Urinary biomarkers in classification of diagnosis and cancer stages",
                  caption = "Data from Silvana Debernardi et al")+ 
  theme(legend.position = "left")