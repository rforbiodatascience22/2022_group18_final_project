# Wrangle data -----------------------------------------------------------------

data_aug_pca <- pca_fit %>% 
  augment(my_data_clean_PCA)

malignant_data_aug_pca <- malignant_PCA_fit %>% 
  augment(malignant_PCA)

# Model data -------------------------------------------------------------------

Kmeans_blood_diagnosis <- data_aug_pca %>%
  select(patient_cohort,
         age,
         sex,
         diagnosis,
         plasma_CA19_9) %>%
  kmeans(centers = 3) %>%
  augment(data_aug_pca) %>% 
  mutate(clusters_org = case_when(cluster_org == 0 ~ "control",
                                  cluster_org == 1 ~ "benign",
                                  cluster_org == 2 ~ "malignant"))

Kmeans_blood_malignant <- malignant_data_aug_pca %>%
  select(patient_cohort,
         age,
         sex,
         stage, 
         plasma_CA19_9) %>%
  kmeans(centers = 4) %>%
  augment(malignant_data_aug_pca) %>% 
  rename(clusters_org = .cluster) %>%
  mutate(clusters_org = case_when(clusters_org == 1 ~ "stage I",
                                  clusters_org == 2 ~ "stage II",
                                  clusters_org == 3 ~ "stage III",
                                  clusters_org == 4 ~ "stage IV"))

Kmeans_urinary_diagnosis <- data_aug_pca %>%
  select(patient_cohort,
         age,
         sex, 
         LYVE1, 
         REG1B, 
         TFF1)  %>%
  kmeans(centers = 4) %>%
  augment(data_aug_pca) %>% 
  rename(clusters_org = .cluster) %>% 
  mutate(clusters_org = case_when(cluster_org == 0 ~ "control",
                                  cluster_org == 1 ~ "benign",
                                  cluster_org == 2 ~ "malignant"))

Kmeans_urinary_malignant <- malignant_data_aug_pca %>%
  select(patient_cohort,
         age,
         sex,
         stage,
         LYVE1, 
         REG1B, 
         TFF1) %>%
  kmeans(centers = 4) %>%
  augment(malignant_data_aug_pca) %>% 
  rename(clusters_org = .cluster) %>%
  mutate(clusters_org = case_when(clusters_org == 1 ~ "stage I",
                                  clusters_org == 2 ~ "stage II",
                                  clusters_org == 3 ~ "stage III",
                                  clusters_org == 4 ~ "stage IV"))

# Visualise data ---------------------------------------------------------------

plasma_CA19_9_alone_plot_diagnosis <- Kmeans_blood_diagnosis %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = clusters_org)) +
  geom_point() +
  theme_minimal(base_size = 10,
                base_family = "Avenir") +
  theme_minimal(8) + 
  theme(legend.position = "bottom",
        plot.title = element_text(size = 8)) +
  labs(title = "Plasma_CA_19_9",
       x = "PC1",
       y = "PC2") +
  ylim(-5, 5) +
  xlim(-6, 6)

plasma_CA19_9_alone_plot_stages <- Kmeans_blood_malignant %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = clusters_org)) +
  geom_point() +
  theme_minimal(base_size = 10,
                base_family = "Avenir") +
  theme_minimal(8) +
  theme(legend.position = "bottom",
        plot.title = element_text(size = 8)) +
  labs(x = "PC1",
       y = "PC2") +
  ylim(-5, 5) +
  xlim(-6, 6)

urinary_alone_diagnosis <- Kmeans_urinary_diagnosis %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = clusters_org)) +
  geom_point() +
  theme_minimal(base_size = 10,
                base_family = "Avenir") +
  theme_minimal(8) + 
  theme(legend.position = "none",
        plot.title = element_text(size = 8)) +
  labs(title = "LYVE1 + REG1B + TFF1",
       x = "PC1",
       y = "PC2") +
  ylim(-5, 5) +
  xlim(-6, 6)

urinary_only_stages <- Kmeans_urinary_malignant %>%
  ggplot(aes(x = .fittedPC1, 
             y = .fittedPC2, 
             colour = clusters_org)) +
  geom_point() +
  theme_minimal(base_size = 10,
                base_family = "Avenir") +
  theme_minimal(8) +
  theme(legend.position = "none",
        plot.title = element_text(size = 8)) +
  labs(x = "PC1",
       y = "PC2") +
  ylim(-5, 5) +
  xlim(-6, 6)

Final_kmeans_plot <- (plasma_CA19_9_alone_plot_diagnosis + urinary_alone_diagnosis) / (plasma_CA19_9_alone_plot_stages + urinary_only_stages)
Final_kmeans_plot + plot_annotation(title = "Comparison between Blood vs Urinary biomarkers in classification of diagnosis and cancer stages",
                                    caption = "Data from Silvana Debernardi et al")