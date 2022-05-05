# Wrangle data -----------------------------------------------------------------

PCA <- my_data_clean_PCA %>% 
  filter(diagnosis != 0) %>%
  select(patient_cohort,
         age,
         sex, 
         stage ,
         plasma_CA19_9,
         creatinine, 
         LYVE1, 
         REG1B, 
         TFF1)

# Model data -------------------------------------------------------------------

pca_fit_ <- PCA  %>%
  select(where(is.numeric)) %>%
  prcomp(scale = TRUE)


data_pca <- pca_fit_ %>% 
  augment(PCA)


# Here i'm testing how the 3 biomarkers can allow a better distinction between the malignant samples and the benign ones 

k_means_PC_1_2 <- data_pca %>%
  select(.fittedPC1, .fittedPC2)  %>%
  kmeans(centers = 2) %>%
  augment(data_pca) %>%
  rename(clusters_pca = .cluster)%>%
  mutate(clusters_pca = case_when(
    clusters_pca == 1 ~ "1",
    clusters_pca == 2 ~ "2"),
    diagnosis = case_when(
      diagnosis == 1 ~ "benign",
      diagnosis == 2 ~ "malignant"))



## Model including only plasma_CA19_9

plasma_CA19_9_alone <- k_means_PC_1_2 %>%
  select(patient_cohort,age,sex,plasma_CA19_9 ) %>%
  kmeans(centers = 2)%>%
  augment(k_means_PC_1_2) %>% 
  rename(clusters = .cluster)



plasma_only_evaluation <- plasma_CA19_9_alone %>% 
  as_tibble %>% 
  bind_cols%>%
  mutate(eval = case_when(
    clusters == 1 & diagnosis == "malignant"~ "TP",
    clusters == 2 & diagnosis == "benign" ~ "TN",
    clusters == 1 & diagnosis == "benign" ~ "FP",
    clusters == 2 & diagnosis == "malignant" ~ "FN"))%>% 
  group_by(eval)%>%
  summarise(count = n())

plasma_only_evaluation
## Model including only LYVE1, REG1B, TFF1 biomarkers


biomarkers_only <- k_means_PC_1_2 %>%
  select( patient_cohort,age,sex, LYVE1, REG1B, TFF1)  %>%
  kmeans(centers = 2)%>%
  augment(k_means_PC_1_2) %>% 
  rename(clusters= .cluster)


biomarkers_only_evaluation <- biomarkers_only %>% 
  as_tibble %>% 
  bind_cols%>%
  mutate(eval = case_when(
    clusters == 2 & diagnosis == "malignant"~ "TP",
    clusters == 1 & diagnosis == "benign" ~ "TN",
    clusters == 2 & diagnosis == "benign" ~ "FP",
    clusters == 1 & diagnosis == "malignant" ~ "FN"))%>% 
  group_by(eval)%>%
  summarise(count = n())
  
View(biomarkers_only_evaluation)

## Full model including plasma_CA19_9 + creatinine + LYVE1 + REG1B + TFF1 

complete_model <- k_means_PC_1_2 %>%
  select( patient_cohort,age,sex, LYVE1, REG1B, TFF1,plasma_CA19_9)  %>%
  kmeans(centers = 2)%>%
  augment(k_means_PC_1_2) %>% 
  rename(clusters = .cluster)


plasma_CA19_9_alone <- plasma_CA19_9_alone %>%  
  mutate(clusters = as.factor(clusters),
         diagnosis = as.factor(diagnosis))


plasma_CA19_9_plot <- ggplot() +
  geom_point( data= plasma_CA19_9_alone , aes(x = .fittedPC1, y = .fittedPC2 , color = diagnosis , size = 1.5  , alpha=1/1000) ) +
  geom_point( data = plasma_CA19_9_alone , aes(x = .fittedPC1, y = .fittedPC2 , color = clusters) ,size = 1   ) +
  scale_color_manual(values=c("blue", "red","red", "blue" ))+
  labs( title = "plasma_CA19_91") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(plot.title = element_text(size=8))+
  theme_minimal(8)

plasma_CA19_9_plot<-plasma_CA19_9_plot+
  guides(size='none',
         alpha='none')+
  theme(legend.position = "bottom")


biomarkers <- biomarkers_only %>%  
  mutate(clusters = as.factor(clusters),
         diagnosis = as.factor(diagnosis))

  biomarkers_plot <- ggplot() +
  geom_point( data= biomarkers , aes(x = .fittedPC1, y = .fittedPC2 , color = diagnosis , size = 1.5  , alpha=1/1000) ) +
  geom_point( data = biomarkers , aes(x = .fittedPC1, y = .fittedPC2 , color = clusters) ,size = 1   ) +
  scale_color_manual(values=c("red", "blue","red", "blue" ))+
  labs( title = "LYVE1+REG1B+TFF1") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(plot.title = element_text(size=8))+
  theme_minimal(8)

  biomarkers_plot<-biomarkers_plot+
  guides(size='none',
         alpha='none')+
  theme(legend.position = "bottom")

  
  complete_model <- complete_model %>%  
    mutate(clusters = as.factor(clusters),
           diagnosis = as.factor(diagnosis))

  complete_modelplot<- ggplot() +
    geom_point( data= complete_model , aes(x = .fittedPC1, y = .fittedPC2 , color = diagnosis , size = 1.5  , alpha=1/1000) ) +
    geom_point( data = complete_model , aes(x = .fittedPC1, y = .fittedPC2 , color = clusters) ,size = 1   ) +
    scale_color_manual(values=c("blue", "red","blue", "red" ))+
    labs( title = "LYVE1+REG1B+TFF1+plasma_CA19_9") +
    theme_minimal(base_size = 10,
                  base_family = "Avenir")+
    theme(plot.title = element_text(size=8))+
    theme_minimal(8)
  
  complete_modelplot<-complete_modelplot+
    guides(size='none',
           alpha='none')+
    theme(legend.position = "bottom")
  
patchwork_diag <-  (plasma_CA19_9_plot |biomarkers_plot)+ 
  plot_layout(guides = 'collect') 
patchwork_diag + plot_annotation(
  title = "K-means results " ,
  caption= "Data from Silvana Debernardi et al" )&
  theme(legend.position = "bottom")



data_aug_pca <- pca_fit %>% 
  augment(my_data_clean_PCA)

malignant_data_aug_pca <- malignant_PCA_fit %>% 
  augment(malignant_PCA)




# Model data -------------------------------------------------------------------



Kmeans_blood_diagnosis <- data_aug_pca %>%
  select(patient_cohort,
         age,
         sex,
         #diagnosis,# the K-means need to guess the diagnosis , i think we are not allowed to give it as input
         plasma_CA19_9) %>%
  kmeans(centers = 3) %>%
  augment(data_aug_pca) %>% 
  rename(clusters= .cluster) %>% 
  mutate(diagnosis = case_when( diagnosis == 0 ~ "control",
                                diagnosis == 1 ~ "benign",
                                diagnosis == 2 ~ "malignant"))

Kmeans_urinary_diagnosis <- data_aug_pca %>%
  select(patient_cohort,
         age,
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

# Visualise data ---------------------------------------------------------------


# How plasma CA19_9 (blood diagnosis ) can discriminate between control benign and malignant individuals


Kmeans_blood_diagnosis_control <-Kmeans_blood_diagnosis %>% 
  filter(diagnosis=="control")

Kmeans_blood_diagnosis_benign <-Kmeans_blood_diagnosis %>% 
  filter(diagnosis=="benign")

Kmeans_blood_diagnosis__malignant <-Kmeans_blood_diagnosis %>% 
  filter(diagnosis=="malignant")

blood_diagnosis_plot <- ggplot() +
  geom_point( data= Kmeans_blood_diagnosis_control , aes(x = .fittedPC1, y = .fittedPC2 , color = clusters, shape = diagnosis ), size = 2) +
  geom_point( data= Kmeans_blood_diagnosis_benign , aes(x = .fittedPC1, y = .fittedPC2 , color = clusters, shape = diagnosis ), size = 2 ) +
  geom_point( data= Kmeans_blood_diagnosis_malignant , aes(x = .fittedPC1, y = .fittedPC2 , color = clusters , shape = diagnosis ), size = 2 ) +
  scale_shape_manual(values=c(1 ,2 ,3))+
  labs( title = "plasma_CA19_9") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(plot.title = element_text(size=8))+
  theme_minimal(8)

blood_diagnosis_plot<-blood_diagnosis_plot+
  guides(size='none')+
  theme(legend.position = "bottom")

# How the 3 biomarkers (urinary diagnosis ) can discriminate between control benign and malignant individuals


Kmeans_urinary_diagnosis_control <-Kmeans_urinary_diagnosis %>% 
  filter(diagnosis=="control")

Kmeans_urinary_diagnosis_benign <- Kmeans_urinary_diagnosis %>% 
  filter(diagnosis=="benign")

Kmeans_urinary_diagnosis_malignant <- Kmeans_urinary_diagnosis %>% 
  filter(diagnosis=="malignant")

urinary_diagnosis_plot <- ggplot() +
  geom_point( data= Kmeans_urinary_diagnosis_control , aes(x = .fittedPC1, y = .fittedPC2 , color = clusters, shape = diagnosis ), size = 2) +
  geom_point( data= Kmeans_urinary_diagnosis_benign , aes(x = .fittedPC1, y = .fittedPC2 , color = clusters, shape = diagnosis ), size = 2 ) +
  geom_point( data= Kmeans_urinary_diagnosis_malignant , aes(x = .fittedPC1, y = .fittedPC2 , color = clusters , shape = diagnosis ), size = 2 ) +
  scale_shape_manual(values=c(1 ,2 ,3))+
  labs( title ="LYVE1 + REG1B + TFF1") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(plot.title = element_text(size=8))+
  theme_minimal(8)

urinary_diagnosis_plot<-urinary_diagnosis_plot+
  guides(size='none')+
  theme(legend.position = "bottom")


kmeans_plot <- ( blood_diagnosis_plot | urinary_diagnosis_plot)
kmeans_plot <- kmeans_plot + plot_annotation(title = "Comparison between Blood vs Urinary biomarkers in classification of diagnosis and cancer stages",
                                                         caption = "Data from Silvana Debernardi et al")





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






