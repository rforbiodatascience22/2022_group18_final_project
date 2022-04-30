#PCA 

library("broom")
library("cowplot")
library("tidyverse")

my_data_clean_PCA <- select(my_data_clean , patient_cohort,age,sex, diagnosis ,stage , plasma_CA19_9,creatinine , LYVE1, REG1B, TFF1) %>%
  mutate(patient_cohort = case_when(
    patient_cohort == "Cohort1" ~ 0,
    patient_cohort == "Cohort2" ~ 1),
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


pca_fit <- my_data_clean_PCA %>% 
  select(where(is.numeric)) %>% # retain only numeric columns
  prcomp(scale = TRUE) # do PCA on scaled data


Data_onto_PCA_plot <- pca_fit %>%
  augment(my_data_clean) %>% # add original dataset back in
  ggplot(aes(.fittedPC1, .fittedPC2, color = diagnosis)) + 
  geom_point(size = 1.5) +
  scale_color_manual(
    values = c(control ="#009E73" , malignant = "#0072B2", benign = "#D55E00")) + 
  theme_half_open(12) + 
  background_grid()+
  theme(legend.position = "right",
        #plot.subtitle =  element_text(size=8),
        plot.title = element_text(size=10))+
  labs( title = "Data projected onto the first two PCs ")
        #subtitle = "The data points have been colored according to the diagnosis status") 


#  rotation matrix plot


arrow_style <- arrow(
  angle = 20, ends = "first", type = "closed", length = grid::unit(8, "pt")
)

rotation_matrix_plot <- pca_fit %>%
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

# Each PC explained variance scree plot

PC_contribution_plot <- pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, percent)) +
  geom_col(fill = "#56B4E9", alpha = 0.8) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.06))) +
  labs( title = "The variance explained by each PC ",
        subtitle = " the first  PC explains nearly 40 % of the total data variance") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "none")+
  theme_minimal_hgrid(12)

# Scree plot =  cumulative explained variance plot 
cumulative_variance_plot <- pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(aes(PC, cumulative)) +
  geom_col(fill = "#0072B2", alpha = 0.8) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.1))) +
    geom_hline(yintercept= 0.64, 
             linetype="dashed", 
             color = "black")+
  labs( title = "The cumulative explained variance ")+
        #subtitle = "Together the first three PCs explains more than 60% of the total data variance") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "none",
        plot.subtitle =  element_text(size=8),
        plot.title = element_text(size=10))+
  theme_minimal_hgrid(12)

patchwork <-  (Data_onto_PCA_plot |  cumulative_variance_plot )
patchwork + plot_annotation(
  title = "PCA results " ,
  subtitle = "The first two PCs can discriminate between patients suffering from PC from the control and benign diagnosed ones  ",
  caption= "Data from Silvana Debernardi et al" )& 
  theme(legend.position = 'bottom')



malignant_PCA <- my_data_clean_PCA %>% 
  filter(diagnosis==2)%>%
  select(
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


malignant_onto_PCA_plot <- malignant_PCA_fit %>%
  augment(malignant_PCA)



malignant_onto_PCA_plot <- malignant_PCA_fit %>%
  augment(malignant_PCA) %>%# add original dataset back in
  mutate(stage=as.factor(stage))%>%
  ggplot(aes(.fittedPC1, .fittedPC2, color = stage)) + 
  geom_point(size = 1.5) +
  scale_color_manual(
    values = c("1" ="#009E73" , "2" = "#0072B2", "3" = "#D55E00" , "4" = "#000000")) + 
  theme_half_open(12) + 
  background_grid()+
  theme(legend.position = "right",
        #plot.subtitle =  element_text(size=8),
        plot.title = element_text(size=10))+
  labs( title = "Data projected onto the first two PCs ")
        #subtitle = "The data points have been colored according to the PDAC stage ") 

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
  geom_hline(yintercept= 0.598, 
             linetype="dashed", 
             color = "black")+
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(
    limits = c(0,1),
    breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1),
    labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.1))) +
  labs( title = "The cumulative explained variance ")+
        #subtitle = "Together the first three PCs explains around 60% of the total data variance") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "none",
        #plot.subtitle =  element_text(size=8),
        plot.title = element_text(size=10))+
  theme_minimal_hgrid(12)



malignant_patchwork <-  (malignant_onto_PCA_plot |  malignant_cumulative_variance_plot )
malignant_patchwork + plot_annotation(
  title = "PCA results " ,
  subtitle = "The first two PCs can't really discriminate between the 4 PDAC stages ",
  caption= "Data from Silvana Debernardi et al" )& 
  theme(legend.position = 'bottom')















