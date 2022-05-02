

## plot  showing the real clusters on PC1 and PC2 space

pl1 <- k_means_onto_PC_1_2 %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = diagnosis)) +
  geom_point() +
  labs( title = "The true clusters ") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "bottom")+
  theme_minimal(8)

#The  clusters predicted by the K means on the whole data plot into PC1 and PC2 space

pl2 <- data_aug_pca_aug_k_means %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = clusters_org)) +
  geom_point() +
  labs( subtitle = "The clusters predicted by K-means performed on the original data plotted into the PC1 and PC2 space ") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "bottom")+
  theme_minimal(8)


#The  clusters predicted by the K means on the data fitted into PC1 and PC2  plotted  into PC1 and PC2 space

pl3 <- k_means_onto_PC_1_2 %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = clusters_pca)) +
  geom_point() +
  labs( subtitle = "The clusters predicted by K-means performed on the data fitted onto PC1 and PC2 plotted into the PC1 and PC2 space ")+
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "bottom")+
  theme_minimal(8)


patchwork <-  (pl1 /pl2 / pl3) + plot_layout(guides = 'collect')
patchwork + plot_annotation(
  title = "K-means results " ,
  caption= "Data from Silvana Debernardi et al" )& 
  theme(legend.position = 'right')



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

## Model including only LYVE1, REG1B, TFF1 biomarkers

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



plasma_CA19_9_alone_plot <- plasma_CA19_9_alone_k_means %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = clusters_org)) +
  geom_point() +
  labs( title = " plasma_CA19_9 ") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "none")+
  theme_minimal(8)

biomarkers_only_plot <- biomarkers_only_k_means %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = clusters_org)) +
  geom_point() +
  labs( title = "LYVE1+REG1B+ TFF1 ") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "none")+
  theme_minimal(8)

complete_model_plot <- complete_model_k_means %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = clusters_org)) +
  geom_point() +
  labs( title = "  LYVE1+REG1B+TFF1+plasma_CA19_9 ") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "none")+
  theme_minimal(8)

patchwork_diagnosis <-  (pl1|plasma_CA19_9_alone_plot |biomarkers_only_plot | complete_model_plot)+ plot_layout(guides = 'collect') 
patchwork_diagnosis + plot_annotation(
  title = "K-means results " ,
  caption= "Data from Silvana Debernardi et al" )&
  theme(legend.position = "bottom")




#K-means to distinguish the 4 stages of malignant patients


malignant_data_aug_pca <- malignant_PCA_fit %>% 
  augment(malignant_PCA)


malignant_aug_pca_aug_k_means <- malignant_data_aug_pca %>%
  select(patient_cohort,age,sex ,stage , plasma_CA19_9,creatinine , LYVE1, REG1B, TFF1) %>%
  kmeans(centers = 4)%>%
  augment(malignant_data_aug_pca) %>% 
  rename(clusters_org = .cluster)%>%
  mutate(clusters_org = case_when(
    clusters_org == 1 ~ "stage I",
    clusters_org == 2 ~ "stage II",
    clusters_org == 3 ~ "stage III", 
    clusters_org == 4 ~ "stage IV"))


#K-means onto the PC1 and PC2

malignant_k_means_onto_PC_1_2 <- malignant_aug_pca_aug_k_means %>%
  select(.fittedPC1, .fittedPC2)  %>%
  kmeans(centers = 4) %>%
  augment(malignant_aug_pca_aug_k_means) %>% 
  rename(clusters_pca = .cluster)%>%
  mutate(clusters_pca = case_when(
    clusters_pca == 1 ~ "stage I",
    clusters_pca == 2 ~ "stage II",
    clusters_pca == 3 ~ "stage III", 
    clusters_pca == 4 ~ "stage IV"),
    stage = case_when(
      stage == 1 ~ "stage I",
      stage == 2 ~ "stage II",
      stage == 3 ~ "stage III", 
      stage == 4 ~ "stage IV"))


#The real clusters plot into the  PC1 and PC2 space
pl1 <- malignant_k_means_onto_PC_1_2 %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = stage)) +
  geom_point() +
  labs( title = "The true clusters ploted into the PC1 and PC2 space ") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "bottom")+
  theme_minimal(8)

#The  clusters predicted by the K means on the whole data plot into PC1 and PC2 space

pl2 <- malignant_aug_pca_aug_k_means %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = clusters_org)) +
  geom_point() +
  labs( subtitle = "The clusters predicted by K-means performed on the original data plotted into the PC1 and PC2 space ") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "bottom")+
  theme_minimal(8)


#The  clusters predicted by the K means on the data fitted into PC1 and PC2  plotted  into PC1 and PC2 space

pl3 <- malignant_k_means_onto_PC_1_2 %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = clusters_pca)) +
  geom_point() +
  labs( subtitle = "The clusters predicted by K-means performed on the data fitted onto PC1 and PC2 plotted into the PC1 and PC2 space ")+
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "bottom")+
  theme_minimal(8)


patchwork <-  (pl1 /pl2 / pl3) 
patchwork + plot_annotation(
  title = "K-means results " ,
  caption= "Data from Silvana Debernardi et al" )& 
  theme(legend.position = 'right')


