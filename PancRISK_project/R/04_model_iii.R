data_aug_pca <- pca_fit %>% 
  augment(my_data_clean_PCA)

#K-means onto the data 
data_aug_pca_aug_k_means <- data_aug_pca %>%
  select(patient_cohort,age,sex, diagnosis ,stage , plasma_CA19_9,creatinine , LYVE1, REG1B, TFF1) %>%
  kmeans(centers = 3)%>%
  augment(data_aug_pca) %>% 
  rename(clusters_org = .cluster)%>%
  mutate(clusters_org = case_when(
    clusters_org == 3 ~ "benign",
    clusters_org == 1 ~ "malignant",
    clusters_org == 2 ~ "control"))

#K-means onto the PC1 and PC2
k_means_onto_PC_1_2 <- data_aug_pca_aug_k_means %>%
  select(.fittedPC1, .fittedPC2)  %>%
  kmeans(centers = 3) %>%
  augment(data_aug_pca_aug_k_means) %>% 
  rename(clusters_pca = .cluster)%>%
  mutate(clusters_pca = case_when(
    clusters_pca == 2 ~ "benign",
    clusters_pca == 1 ~ "malignant",
    clusters_pca == 3 ~ "control"),
    diagnosis = case_when(
      diagnosis == 0 ~ "control",
      diagnosis == 1 ~ "benign",
      diagnosis == 2 ~ "malignant"))

#The real clusters plot into the  PC1 and PC2 space
pl1 <- k_means_onto_PC_1_2 %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = diagnosis)) +
  geom_point() +
  labs( title = "The true clusters ploted into the PC1 and PC2 space ") +
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


patchwork <-  (pl1 /pl2 / pl3) 
patchwork + plot_annotation(
  title = "K-means results " ,
  caption= "Data from Silvana Debernardi et al" )& 
  theme(legend.position = 'right')



patchwork





