data_aug_pca <- pca_fit %>% 
  augment(my_data_clean_PCA)%>%
  filter(.fittedPC1> -6) # removing an outtlier for better K-means 


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
  mutate(diagnosis = case_when( diagnosis == 0 ~ "control",
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

blood_malignant_statistics <-Kmeans_blood_diagnosis_malignant %>%
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
  filter(diagnosis=="control")

Kmeans_urinary_diagnosis_benign <- Kmeans_urinary_diagnosis_ %>% 
  filter(diagnosis=="benign")

Kmeans_urinary_diagnosis_malignant <- Kmeans_urinary_diagnosis_ %>% 
  filter(diagnosis=="malignant")

control_statistics <-Kmeans_urinary_diagnosis_control %>%
  group_by(clusters) %>%
  summarise(count = n())

benign_statistics <-Kmeans_urinary_diagnosis_benign %>%
  group_by(clusters) %>%
  summarise(count = n())

malignant_statistics <-Kmeans_urinary_diagnosis_malignant %>%
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
  rename(clusters= .cluster) %>% 
  mutate(diagnosis = case_when( diagnosis == 0 ~ "control",
                                diagnosis == 1 ~ "benign",
                                diagnosis == 2 ~ "malignant"))

Kmeans_complete_diagnosis_control <- Kmeans_complete_diagnosis_ %>% 
  filter(diagnosis=="control")

Kmeans_complete_diagnosis_benign <- Kmeans_complete_diagnosis_ %>% 
  filter(diagnosis == "benign")

Kmeans_complete_diagnosis_malignant <- Kmeans_complete_diagnosis_ %>% 
  filter(diagnosis == "malignant")

complete_control_statistics <-Kmeans_complete_diagnosis_control %>%
  group_by(clusters) %>%
  summarise(count = n())

complete_benign_statistics <-Kmeans_complete_diagnosis_benign %>%
  group_by(clusters) %>%
  summarise(count = n())

complete_malignant_statistics <-Kmeans_complete_diagnosis_malignant %>%
  group_by(clusters) %>%
  summarise(count = n())



# Visualise data ---------------------------------------------------------------


# How plasma CA19_9 (blood diagnosis ) can discriminate between control benign and malignant individuals


blood_diagnosis_plot <- ggplot() +
  geom_point( data= Kmeans_blood_diagnosis_control , 
              aes(x = .fittedPC1, 
                  y = .fittedPC2 ,
                  color = clusters, 
                  shape = diagnosis), 
              size = 2) +
  geom_point( data= Kmeans_blood_diagnosis_benign , 
              aes(x = .fittedPC1, 
                  y = .fittedPC2 , 
                  color = clusters, 
                  shape = diagnosis ),
              size = 2 ) +
  geom_point( data = Kmeans_blood_diagnosis_malignant, 
              aes(x = .fittedPC1, 
                  y = .fittedPC2, 
                  color = clusters, 
                  shape = diagnosis), 
              size = 2) +
  scale_shape_manual(values=c(1, 2, 3)) +
  xlim(-6, 3.5) +
  ylim(-5, 5) +
  labs(title = "plasma_CA19_9") +
  theme_minimal(base_size = 10,
                base_family = "Avenir") +
  theme(plot.title = element_text(size=8)) +
  theme_minimal(8)


blood_diagnosis_plot_<-blood_diagnosis_plot +
  guides(size='none') +
  theme(legend.position = "bottom ",
        legend.box="vertical", 
        legend.margin=margin())



# How the 3 biomarkers (urinary diagnosis ) can discriminate between control benign and malignant individuals


urinary_diagnosis_plot <- ggplot() +
  geom_point(data = Kmeans_urinary_diagnosis_control, 
              aes(x = .fittedPC1, 
                  y = .fittedPC2, 
                  color = clusters, 
                  shape = diagnosis), 
              size = 2) +
  geom_point( data= Kmeans_urinary_diagnosis_benign, 
              aes(x = .fittedPC1, 
                  y = .fittedPC2, 
                  color = clusters,
                  shape = diagnosis ), 
              size = 2 ) +
  geom_point( data= Kmeans_urinary_diagnosis_malignant , 
              aes(x = .fittedPC1, 
                  y = .fittedPC2 , 
                  color = clusters , 
                  shape = diagnosis ), 
              size = 2 ) +
  xlim(-6, 3.5)+
  ylim(-5, 5)+
  scale_shape_manual(values=c(1 ,2 ,3))+
  labs( title ="LYVE1 + REG1B + TFF1") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(plot.title = element_text(size=8))+
  theme_minimal(8)

urinary_diagnosis_plot_<-urinary_diagnosis_plot+
  guides(size='none')+
  theme(legend.box="vertical", 
        legend.margin=margin(),
        legend.position ="bottom")



kmeans_plot_ <- ( blood_diagnosis_plot_ | urinary_diagnosis_plot_)+
  plot_annotation(title = "Comparison between Blood vs Urinary biomarkers in classification of diagnosis and cancer stages",
                  caption = "Data from Silvana Debernardi et al")+ 
  theme( legend.position = "left")





# Model data -------------------------------------------------------------------

Kmeans_blood_diagnosis <- data_aug_pca %>%
  select(patient_cohort,
         age,
         sex,
         diagnosis,
         plasma_CA19_9) %>%
  kmeans(centers = 3) %>%
  augment(data_aug_pca) %>% 
  rename(clusters_org = .cluster) %>% 
  mutate(clusters_org = case_when(diagnosis == 0 ~ "control",
                                  diagnosis == 1 ~ "benign",
                                  diagnosis == 2 ~ "malignant"))

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
  mutate(clusters_org = case_when(diagnosis == 0 ~ "control",
                                  diagnosis == 1 ~ "benign",
                                  diagnosis == 2 ~ "malignant"))

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
Final_kmeans_plot <- Final_kmeans_plot + plot_annotation(title = "Comparison between Blood vs Urinary biomarkers in classification of diagnosis and cancer stages",
                                                         caption = "Data from Silvana Debernardi et al")

top_plot <- (Data_onto_PCA_plot + cumulative_variance_plot) + (plasma_CA19_9_alone_plot_diagnosis + urinary_alone_diagnosis)
bottom_plot <- (malignant_onto_PCA_plot + malignant_cumulative_variance_plot) + (plasma_CA19_9_alone_plot_stages + urinary_only_stages)






