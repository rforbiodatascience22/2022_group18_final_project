data_aug_pca <- pca_fit %>% 
  augment(my_data_clean_PCA)


# This model includes all the data features 

#K-means onto the data 
data_aug_pca_aug_k_means <- data_aug_pca %>%
  select(patient_cohort,age,sex, plasma_CA19_9,creatinine , LYVE1, REG1B, TFF1) %>%
  kmeans(centers = 3)%>%
  augment(data_aug_pca) %>% 
  rename(clusters_org = .cluster)%>%
  mutate(clusters_org = case_when(
    clusters_org == 1 ~ "1",
    clusters_org == 2 ~ "2",
    clusters_org == 3 ~ "3"))

#K-means onto the PC1 and PC2
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



# Let's try to perform the K-means including only a subset of the features 
#plasma_CA19_9 alone
#LYVE1, REG1B, TFF1
#and the complete model which have already been done above 



plasma_CA19_9_alone_k_means <- data_aug_pca %>%
  select(patient_cohort,age,sex,plasma_CA19_9 ) %>%
  kmeans(centers = 3)%>%
  augment(data_aug_pca) %>% 
  rename(clusters_org = .cluster)

biomarkers_only_k_means <- data_aug_pca %>%
  select(patient_cohort,age,sex, LYVE1, REG1B, TFF1)  %>%
  kmeans(centers = 3)%>%
  augment(data_aug_pca) %>% 
  rename(clusters_org = .cluster)

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
  theme(legend.position = "bottom")+
  theme_minimal(8)

biomarkers_only_plot <- biomarkers_only_k_means %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = clusters_org)) +
  geom_point() +
  labs( title = "LYVE1+REG1B+ TFF1 ") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "bottom")+
  theme_minimal(8)

complete_model_plot <- complete_model_k_means %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = clusters_org)) +
  geom_point() +
  labs( title = "  LYVE1+REG1B+TFF1+plasma_CA19_9 ") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "bottom")+
  theme_minimal(8)

patchwork_ <-  (pl1|plasma_CA19_9_alone_plot |biomarkers_only_plot | complete_model_plot) 
patchwork_ + plot_annotation(
  title = "K-means results " ,
  caption= "Data from Silvana Debernardi et al" )& 
  theme(legend.position = 'bottom')




#PCA only for the malignant samples 

malignant_PCA <- my_data_clean_PCA %>% 
  filter(diagnosis==2)%>%
  select(patient_cohort,
         age,
         sex, 
         stage ,
         plasma_CA19_9,
         creatinine , 
         LYVE1, 
         REG1B, 
         TFF1)

malignant_PCA_fit <- malignant_PCA %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  prcomp(scale = TRUE) # do PCA on scaled data


malignant_pcs <- malignant_PCA_fit %>% 
  broom::tidy("pcs")

malignant_pcs

malignant_onto_PCA_plot <- malignant_PCA_fit %>%
  augment(malignant_PCA)
View(malignant_onto_PCA_plot)


malignant_onto_PCA_plot <- malignant_PCA_fit %>%
  augment(malignant_PCA) %>%# add original dataset back in
  mutate(stage=as.factor(stage))%>%
  ggplot(aes(.fittedPC1, .fittedPC2, color = stage)) + 
  geom_point(size = 1.5) +
  scale_color_manual(
    values = c("1" ="#009E73" , "2" = "#0072B2", "3" = "#D55E00" , "4" = "#000000")) + 
  theme_half_open(12) + 
  background_grid()+
  theme(legend.position = "right")+
  labs( title = "Data projected onto the first two PCs ",
        subtitle = "The data points have been colored according to the PDAC stage ") 



#  rotation matrix plot
malignant_PCA_fit %>%
  tidy(matrix = "rotation")

arrow_style <- arrow(
  angle = 20, ends = "first", type = "closed", length = grid::unit(8, "pt")
)

malignant_rotation_matrix_plot <- malignant_PCA_fit %>%
  tidy(matrix = "rotation") %>%
  pivot_wider(names_from = "PC", names_prefix = "PC", values_from = "value") %>%
  ggplot(aes(PC1, PC2)) +
  geom_segment(xend = 0, yend = 0, arrow = arrow_style) +
  geom_text(
    aes(label = column),
    hjust = 1, nudge_x = -0.02, 
    color = "#904C2F"
  ) +
  xlim(-1.25, .5) + ylim(-.5, 1) +
  coord_fixed() + # fix aspect ratio to 1:1
  theme_minimal_grid(12)

malignant_rotation_matrix_plot

# Each PC explained variance scree plot

malignant_PC_contribution_plot <- malignant_PCA_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.06))) +
  labs( title = "The variance explained by each PC ",
        subtitle = " the first  PC explains nearly 30 % of the total data variance") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "none")+
  theme_minimal_hgrid(12)



# Scree plot =  cumulative explained variance plot 
malignant_cumulative_variance_plot <- malignant_PCA_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, cumulative)) +
  geom_col(fill = "#0072B2", alpha = 0.8) +
  geom_hline(yintercept= 0.574, 
             linetype="dashed", 
             color = "black")+
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(
    limits = c(0,1),
    breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.1))) +
  labs( title = "The cumulative explained variance ",
        subtitle = "Together the first three PCs explains around 60% of the total data variance") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "none")+
  theme_minimal_hgrid(12)



malignant_patchwork <-  (malignant_onto_PCA_plot |  malignant_cumulative_variance_plot )
malignant_patchwork + plot_annotation(
  title = "PCA results " ,
  subtitle = "The first two PCs can't really discriminate between the 4 PDAC stages ",
  caption= "Data from Silvana Debernardi et al" )& 
  theme(legend.position = 'bottom')

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


