library("tidyverse")
library("patchwork")
library("ggridges")


# Loading of data 
clean_data <- read.csv("../data/02_my_data_clean.tsv", sep = "\t")




# Age & Sex distribution 

age_distribution <- clean_data %>% 
  ggplot(aes(x = age, 
             color = diagnosis)) +
  geom_density() + 
  theme_minimal(base_size = 10) + 
  labs(title = "Age Distribution Across Different Diagnosis Groups", 
       x = "Age", 
       y = "Density", 
       color = "Diagnosis Group") + 
  scale_color_discrete(labels = c("Begnin", "Control", "Malignant"))


sex_distribution <- ggplot(data = clean_data, 
                           mapping = aes(x = diagnosis, 
                                         fill = sex)) + 
  geom_bar(stat="count", 
           position = position_dodge()) + 
  geom_text(stat = "count", 
            aes(label= ..count..), 
            position = position_dodge(1), 
            vjust = -0.3, 
            size = 2) + 
  theme_minimal(base_size = 10) + 
  theme(legend.position = "right") + 
  labs(title = "Sex Distribution Across Different Diagnosis Groups", 
       x = "Diagnosis Group", 
       y = "Count", 
       fill = "Sex") + 
  scale_x_discrete(labels = c("Benign", "Control", "PDAC"))



# Stage Distribution 

diagnosis_percentage <- clean_data %>%
  count(diagnosis) %>% 
  mutate(Diagnosis = "Diagnosis Group", 
         percentage = round(n/sum(n)*100), digits = 2) %>% 
  ggplot(mapping = aes(fill = diagnosis, 
                       x = n, 
                       y = Diagnosis)) + 
  geom_col() + 
  scale_fill_brewer(palette = "Pastel2", 
                    labels = c("Begnin", "Control", "Malignant")) + 
  geom_text(aes(label = paste0(percentage, "%")),
            position = position_stack(vjust = 0.5), 
            size = 4) +
  theme_minimal(base_size = 10) + 
  labs(x = NULL, 
       y = NULL, 
       title = "Diagnosis Group Distribution", 
       fill = "Diagnosis Group") + 
  theme(legend.key.height = unit(0.01, 'cm'), 
        legend.position = "right") 
#scale_fill_discrete(labels = c("Begnin", "Control", "Malignant"))


# Stacked + percent
stage_precentage <- clean_data %>% 
  filter(!is.na(stage)) %>% 
  count(stage) %>% 
  mutate(Stage = "Stage", 
         proportion = round(n/sum(n)*100), digits = 2) %>% 
  ggplot(mapping = aes(fill = stage, 
                       x = n, 
                       y = Stage)) + 
  geom_col() + 
  scale_fill_brewer(palette = "Set2") + 
  geom_text(aes(label = paste0(proportion, "%")),
            position = position_stack(vjust = 0.5), 
            size = 3) +
  theme_minimal(base_size = 10) + 
  labs(x = "Number of Patients", 
       y = NULL, 
       title = "Malignant Patients Stages Distribution", 
       fill = "Stage") + 
  theme(legend.key.height = unit(0.01, 'cm'), 
        legend.position = "right")



all_distribution <- (age_distribution / sex_distribution / diagnosis_percentage / stage_precentage) + 
  plot_layout(heights = unit(c(0.2, 0.1, 0.03, 0.03), c('npc', 'npc', 'npc', 'npc')))

all_distribution



#==================================================================================

# Sample distribution between cohorts regarding protein samples taken

data_long <- clean_data %>% 
  select(patient_cohort, 
         plasma_CA19_9, 
         creatinine, 
         LYVE1, 
         REG1A, 
         REG1B, 
         TFF1) %>% 
  pivot_longer(!patient_cohort, 
               names_to = "protein", 
               values_to = "expression")

data_long %>%
  group_by(patient_cohort) %>%
  ggplot(mapping = aes(y = patient_cohort, 
                       fill = patient_cohort)) +
  facet_wrap(~protein, 
             ncol = 2, 
             scales = "free") +
  geom_bar() +
  xlab("Occurence of measured protein levels per Cohort ") +
  ylab("Cohort of patients") +
  labs(title = "Sample distribution between two cohorts of patients") +
  theme_minimal() +
  theme(legend.position="none")


#  Levels of the 3 biomarkers in control, benign, and  (PDAC) samples

LYVEvscreatinine <- select(clean_data, diagnosis , stage , LYVE1, creatinine) %>%
  mutate(diagnosis = case_when(
    diagnosis == "control" ~ "control",
    diagnosis == "benign" ~ "benign",
    stage  == "IA" ~ "PDAC I-II",
    stage  == "IB" ~ "PDAC I-II",
    stage  == "IIA" ~ "PDAC I-II",
    stage  == "IIB" ~ "PDAC I-II",
    stage  == "III" ~ "PDAC III-IV",
    stage  == "IV" ~ "PDAC III-IV") , 
    LYVE1 = LYVE1 / creatinine  )


LYVEvscreatinine_plot <-ggplot(data = LYVEvscreatinine, 
                               mapping = aes(x = diagnosis, 
                                             y = LYVE1, 
                                             color = diagnosis)) + 
  geom_violin(trim=FALSE)+ 
  geom_boxplot(width=0.2)+
  labs(y = "LYVE1 ng/mg creatinine" , 
       title = "LYVE1") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "none")


REG1Bvscreatinine <- select(clean_data, diagnosis , stage , REG1B, creatinine) %>%
  mutate(diagnosis = case_when(
    diagnosis == "control" ~ "control",
    diagnosis == "benign" ~ "benign",
    stage  == "IA" ~ "PDAC I-II",
    stage  == "IB" ~ "PDAC I-II",
    stage  == "IIA" ~ "PDAC I-II",
    stage  == "IIB" ~ "PDAC I-II",
    stage  == "III" ~ "PDAC III-IV",
    stage  == "IV" ~ "PDAC III-IV"), 
    REG1B = REG1B/ creatinine)


REG1Bvscreatinine_plot <- ggplot(data = REG1Bvscreatinine, 
                                 mapping = aes(x = diagnosis, 
                                               y = REG1B, 
                                               color = diagnosis)) + 
  geom_violin(trim=FALSE)+ 
  geom_boxplot(width=0.2)+
  labs(y = "REG1B ng/mg creatinine" , 
       title = "REG1B") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "none")




TFF1vscreatinine <- select(clean_data, diagnosis , stage , TFF1, creatinine) %>%
  mutate(diagnosis = case_when(
    diagnosis == "control" ~ "control",
    diagnosis == "benign" ~ "benign",
    stage  == "IA" ~ "PDAC I-II",
    stage  == "IB" ~ "PDAC I-II",
    stage  == "IIA" ~ "PDAC I-II",
    stage  == "IIB" ~ "PDAC I-II",
    stage  == "III" ~ "PDAC III-IV",
    stage  == "IV" ~ "PDAC III-IV"),
    TFF1 = TFF1/ creatinine)





TFF1vscreatinine_plot <-ggplot(data = TFF1vscreatinine, 
                               mapping = aes(x = diagnosis, 
                                             y = TFF1, 
                                             color = diagnosis)) + 
  geom_violin(trim=FALSE)+ 
  geom_boxplot(width=0.2)+
  labs(y = "TFF1 ng/mg creatinine" , 
       title = "TFF1") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "none")





violin_plot <- (LYVEvscreatinine_plot |REG1Bvscreatinine_plot |TFF1vscreatinine_plot )
violin_plot + plot_annotation(
  title = "The levels of the 3 biomarkers in control, benign, and pancreatic ductal adenocarcinoma (PDAC) samples. " ,
  subtitle = "Violin plots are shown for each protein. All data were creatinine normalised ",
  caption= "Data from Silvana Debernardi et al" ) +
  theme(legend.position = 'none')


