# Wrangle data -----------------------------------------------------------------

malignant_PCA <- my_data_clean_PCA %>% 
  filter(diagnosis == 2) %>%
  select(patient_cohort,
         age,
         sex, 
         stage ,
         plasma_CA19_9,
         creatinine, 
         LYVE1, 
         REG1B, 
         TFF1)


arrow_style <- arrow(angle = 20, 
                     ends = "first", 
                     type = "closed", 
                     length = grid::unit(8, "pt"))

# Model data -------------------------------------------------------------------

pca_fit <- my_data_clean_PCA  %>%
  select(where(is.numeric)) %>%
  prcomp(scale = TRUE)

malignant_PCA_fit <- malignant_PCA %>% 
  select(where(is.numeric)) %>%
  prcomp(scale = TRUE)

malignant_pcs <- malignant_PCA_fit %>% 
  broom::tidy("pcs")

# Data visualization -----------------------------------------------------------

Data_onto_PCA_plot <- pca_fit %>%
  augment(my_data_clean_PCA) %>%
  mutate(diagnosis = case_when(
    diagnosis == 0 ~ "control",
    diagnosis == 1 ~ "benign",
    diagnosis == 2 ~ "malignant")) %>%
  ggplot(mapping = aes(x = .fittedPC1, 
                       y = .fittedPC2, 
                       color = diagnosis)) + 
  geom_point(size = 1.5) +
  background_grid() +
  theme(legend.position = "right",
        plot.title = element_text(size = 10)) +
  theme_minimal() + 
  labs(x = "PC1", 
       y = "PC2",
       title = "Biplot : data onto PC1 and PC2") +
  xlim(-6, 3.5) +
  ylim(-5, 5)


rotation_matrix_plot <- pca_fit %>%
  tidy(matrix = "rotation") %>%
  pivot_wider(names_from = "PC", 
              names_prefix = "PC", 
              values_from = "value") %>%
  ggplot(mapping = aes(x = PC1, 
                       y = PC2)) +
  geom_segment(xend = 0, 
               yend = 0, 
               arrow = arrow_style) +
  geom_text(aes(label = column), 
            hjust = 1, 
            nudge_x = -0.02,
            color = "#904C2F") +
  xlim(-1.25, .5) + 
  ylim(-.5, 1) +
  coord_fixed() +
  theme_minimal_grid(12)

PC_contribution_plot <- pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(mapping = aes(x = PC, 
                       y = percent)) +
  geom_col(fill = "#56B4E9", 
           alpha = 0.8) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(labels = scales::percent_format(),
    expand = expansion(mult = c(0, 0.06))) +
  labs(title = "The variance explained by each PC",
        subtitle = "the first  PC explains nearly 40 % of the total data variance") +
  theme_minimal() +
  theme(legend.position = "none") +
  theme_minimal_hgrid(12)

cumulative_variance_plot <- pca_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(mapping = aes(x = PC, 
                       y = cumulative)) +
  geom_col(fill = "#0072B2", 
           alpha = 0.8) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.1))) +
    geom_hline(yintercept = 0.64, 
             linetype = "dashed", 
             color = "black") +
  theme(legend.position = "none",
        plot.title = element_text(size = 10)) +
  theme_minimal() +
  background_grid() +
  labs(title = "Scree plot : Cumulative variance ")


malignant_onto_PCA_plot <- malignant_PCA_fit %>%
  augment(malignant_PCA) %>% 
  mutate(stage = as.factor(stage)) %>%
  ggplot(mapping = aes(x = .fittedPC1,
                       y = .fittedPC2, 
                     color = stage)) + 
  geom_point(size = 1.5) +
  theme(legend.position = "none",
        plot.title = element_text(size = 10)) +
  theme_minimal() + 
  background_grid() +
  labs(x = "PC1", 
       y = "PC2") +
  xlim(-6, 3.5) +
  ylim(-5, 5)


malignant_rotation_matrix_plot <- malignant_PCA_fit %>%
  tidy(matrix = "rotation") %>%
  pivot_wider(names_from = "PC", 
              names_prefix = "PC", 
              values_from = "value") %>%
  ggplot(aes(PC1, PC2)) +
  geom_segment(xend = 0, 
               yend = 0, 
               arrow = arrow_style) +
  geom_text(aes(label = column),
            hjust = 1, 
            nudge_x = -0.02, 
            color = "#904C2F") +
  xlim(-1.25, .5) + 
  ylim(-.5, 1) +
  coord_fixed() +
  theme_minimal_grid(12)


malignant_PC_contribution_plot <- malignant_PCA_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(mapping = aes(x = PC,
                       y = percent)) +
  geom_col(fill = "#56B4E9",
           alpha = 0.8) +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.06))) +
  theme_minimal() +
  theme(legend.position = "none") +
  theme_minimal_hgrid(12) +
  labs(title = "The variance explained by each PC",
       subtitle = "the first  PC explains nearly 30 % of the total data variance")


malignant_cumulative_variance_plot <- malignant_PCA_fit %>%
  tidy(matrix = "eigenvalues") %>%
  ggplot(mapping = aes(x = PC, 
                       y = cumulative)) +
  geom_col(fill = "#0072B2", 
           alpha = 0.8) +
  geom_hline(yintercept = 0.598, 
             linetype = "dashed", 
             color = "black") +
  scale_x_continuous(breaks = 1:10) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.1))) +
  theme(legend.position = "none",
        plot.title = element_text(size = 10)) +
  theme_minimal() +
  background_grid()
  


diagnosis_PCA <- (Data_onto_PCA_plot + cumulative_variance_plot)+ 
  plot_annotation(title = "Principal Component Analysis of PancRISK diagnosis")

diagnosis_PCA

PCA_Final_plot <- (Data_onto_PCA_plot + cumulative_variance_plot) / (malignant_onto_PCA_plot + malignant_cumulative_variance_plot)
PCA_Final_plot <- PCA_Final_plot + plot_annotation(title = "Principal Component Analysis of PancRISK diagnosis")

