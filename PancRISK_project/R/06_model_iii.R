## Augmenting the clean data with the data fitted on each PC (will be used in for the K-means )

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

pca_fit <- PCA  %>%
  select(where(is.numeric)) %>%
  prcomp(scale = TRUE)


data_pca <- pca_fit %>% 
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

write_tsv(x = data_aug_pca_aug_k_means,
          file = "PancRISK_project/data/05_k_means_on_whole_data.tsv")
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

write_tsv(x = k_means_onto_PC_1_2,
          file = "PancRISK_project/data/06_k_means_onto_PC_1_2.tsv")
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


## plot  showing the real clusters on PC1 and PC2 space

pl1 <- k_means_onto_PC_1_2 %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = diagnosis)) +
  geom_point() +
  labs( title = "The true clusters ") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "bottom",
                plot.title = element_text(size=8))+
  theme_minimal(8)

#The  clusters predicted by the K means on the whole data plot into PC1 and PC2 space

pl2 <- data_aug_pca_aug_k_means %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = clusters_org)) +
  geom_point() +
  labs( subtitle = "The clusters predicted by K-means performed on the original data plotted into the PC1 and PC2 space ") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "bottom",
        plot.title = element_text(size=8))+
  theme_minimal(8)


#The  clusters predicted by the K means on the data fitted into PC1 and PC2  plotted  into PC1 and PC2 space

pl3 <- k_means_onto_PC_1_2 %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = clusters_pca)) +
  geom_point() +
  labs( subtitle = "The clusters predicted by K-means performed on the data fitted onto PC1 and PC2 plotted into the PC1 and PC2 space ")+
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "bottom",
        plot.title = element_text(size=8))+
  theme_minimal(8)


patchwork <-  (pl1 /pl2 / pl3) + plot_layout(guides = 'collect')
patchwork + plot_annotation(
  title = "K-means results " ,
  caption= "Data from Silvana Debernardi et al" )& 
  theme(legend.position = 'right')



plasma_CA19_9_alone_plot <- plasma_CA19_9_alone_k_means %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = clusters_org)) +
  geom_point() +
  labs( title = " plasma_CA19_9 ") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "none",
        plot.title = element_text(size=8))+
  theme_minimal(8)

biomarkers_only_plot <- biomarkers_only_k_means %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = clusters_org)) +
  geom_point() +
  labs( title = "LYVE1+REG1B+ TFF1 ") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "none",
        plot.title = element_text(size=8))+
  theme_minimal(8)

complete_model_plot <- complete_model_k_means %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = clusters_org)) +
  geom_point() +
  labs( title = "  LYVE1+REG1B+TFF1+plasma_CA19_9 ") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "none",
        plot.title = element_text(size=8))+
  theme_minimal(8)

patchwork_diagnosis <-  (pl1|plasma_CA19_9_alone_plot |biomarkers_only_plot | complete_model_plot)+ 
  plot_layout(guides = 'collect') 
patchwork_diagnosis + plot_annotation(
  title = "K-means results " ,
  caption= "Data from Silvana Debernardi et al" )&
  theme(legend.position = "bottom")




#--------------------



#--------------------

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
pl1_malignant<- malignant_k_means_onto_PC_1_2 %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = stage)) +
  geom_point() +
  labs( title = "The true clusters ploted into the PC1 and PC2 space ") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "bottom")+
  theme_minimal(8)

#The  clusters predicted by the K means on the whole data plot into PC1 and PC2 space

pl2_malignant <- malignant_aug_pca_aug_k_means %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = clusters_org)) +
  geom_point() +
  labs( subtitle = "The clusters predicted by K-means performed on the original data plotted into the PC1 and PC2 space ") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "bottom")+
  theme_minimal(8)


#The  clusters predicted by the K means on the data fitted into PC1 and PC2  plotted  into PC1 and PC2 space

pl3_malignant <- malignant_k_means_onto_PC_1_2 %>%
  ggplot(aes(x = .fittedPC1, y = .fittedPC2, colour = clusters_pca)) +
  geom_point() +
  labs( subtitle = "The clusters predicted by K-means performed on the data fitted onto PC1 and PC2 plotted into the PC1 and PC2 space ")+
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "bottom")+
  theme_minimal(8)


patchwork <-  (pl1_malignant /pl2_malignant / pl3_malignant) 
patchwork + plot_annotation(
  title = "K-means results " ,
  caption= "Data from Silvana Debernardi et al" )& 
  theme(legend.position = 'right')


