library("tidyverse")
library("patchwork")
library("ggridges")

# Loading of data 
clean_data <- read.csv("PancRISK_project/data/02_my_data_clean.tsv", sep = "\t")

# Age & Sex distribution 

age_distribution <- 
  clean_data %>% 
  ggplot(aes(x = age, 
             color = fct_rev(diagnosis),
             fill = fct_rev(diagnosis))) +
  geom_density(alpha = 0.2) + 
  theme_minimal(base_size = 10) + 
  labs(title = "Age Distribution Across Different Diagnosis Groups", 
       x = "Age", 
       y = "Density", 
       color = "Diagnosis",
       fill = ""
  ) + 
  guides(fill = FALSE) +
  scale_color_brewer(palette = "Pastel2",
                     guide = guide_legend(reverse = TRUE),
                     labels = c("Malignant", "Control", "Benign")) +
  scale_fill_brewer(palette = "Pastel2",
                    guide = guide_legend(reverse = TRUE),
                    labels = c("Malignant", "Control", "Benign"))

sex_distribution <- 
  clean_data %>%
  ggplot(mapping = aes(x = diagnosis, 
                       fill = sex)) + 
  geom_bar(stat="count", 
           position = position_dodge()) + 
  geom_text(stat = "count", 
            aes(label= ..count..), 
            position = position_dodge(1), 
            vjust = -0.3, 
            size = 3) + 
  ylim(0, 140) +
  theme_minimal(base_size = 10) + 
  theme(legend.position = "right") + 
  labs(title = "Sex Distribution Across Different Diagnosis Groups", 
       x = "Diagnosis Group", 
       y = "Count", 
       fill = "Sex") + 
  scale_x_discrete(labels = c("Benign", "Control", "Malignant"))

# Stage Distribution 

diagnosis_percentage <- 
  clean_data %>%
  count(diagnosis) %>% 
  mutate(Diagnosis = "", 
         percentage = round(n/sum(n)*100), digits = 2) %>% 
  ggplot(mapping = aes(fill = fct_rev(diagnosis), 
                       x = n, 
                       y = Diagnosis)) + 
  geom_col() + 
  scale_fill_discrete(labels = c("Benign", "Control", "Malignant")) +
  scale_fill_brewer(palette = "Pastel2") +
  geom_text(aes(label = paste0(percentage, "%", " - ", str_to_sentence(diagnosis))),
            position = position_stack(vjust = 0.5), 
            size = 3) +
  theme_minimal(base_size = 10) + 
  labs(x = "", 
       y = "", 
       title = "Diagnosis Group Distribution", 
       fill = "") + 
  theme(legend.position = "none")

# Stacked + percent
stage_precentage <- 
  clean_data %>% 
  filter(!is.na(stage)) %>% 
  count(stage) %>% 
  mutate(Stage = "", 
         proportion = round(n/sum(n)*100), digits = 2) %>% 
  ggplot(mapping = aes(fill = stage, 
                       x = n, 
                       y = Stage)) + 
  geom_col(position = position_stack(reverse = TRUE)) + 
  scale_fill_brewer(palette = "Pastel2") + 
  geom_text(aes(label = ifelse(proportion >= 4, paste0(proportion, "%"), "")),
            position = position_stack(vjust = 0.5,
                                      reverse = TRUE), 
                                      size = 3) +
  theme_minimal(base_size = 10) + 
  labs(x = "Number of Patients", 
       y = "Stage", 
       title = "Malignant Patients Stages Distribution", 
       fill = "Stage") + 
  theme(legend.key.height = unit(0.001, 'cm'), 
        legend.position = "right") +
  guides(fill = guide_legend(ncol = 2))

all_distribution <- (age_distribution / sex_distribution / plot_spacer() / diagnosis_percentage / plot_spacer() / stage_precentage) + 
  plot_layout(heights = unit(c(0.15, 0.15, 0.02, 0.03, 0.005, 0.03), c('npc', 'npc', 'npc', 'npc')))

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

# data_long %>%
#   group_by(patient_cohort) %>%
#   ggplot(mapping = aes(y = patient_cohort, 
#                        fill = patient_cohort)) +
#   facet_wrap(~protein, 
#              ncol = 2, 
#              scales = "free") +
#   geom_bar() +
#   xlab("Occurence of measured protein levels per Cohort ") +
#   ylab("Cohort of patients") +
#   labs(title = "Sample distribution between two cohorts of patients") +
#   theme_minimal() +
#   theme(legend.position="none")

#  Levels of the 3 biomarkers in control, benign, and  (PDAC) samples

LYVEvscreatinine <- clean_data%>%
  select( diagnosis , stage , LYVE1, creatinine) %>%
  drop_na(diagnosis)  %>%
  mutate(diagnosis = case_when(
    diagnosis == "control" ~ "control",
    diagnosis == "benign" ~ "benign",
    stage  == "I" ~ "PDAC I-II",
    stage  == "IA" ~ "PDAC I-II",
    stage  == "IB" ~ "PDAC I-II",
    stage  == "II" ~ "PDAC I-II",
    stage  == "IIA" ~ "PDAC I-II",
    stage  == "IIB" ~ "PDAC I-II",
    stage  == "III" ~ "PDAC III-IV",
    stage  == "IV" ~ "PDAC III-IV") , 
    LYVE1 = LYVE1 / creatinine)


LYVEvscreatinine_plot <-
  ggplot(data = LYVEvscreatinine, 
         mapping = aes(x = diagnosis, 
                       y = LYVE1, 
                       color = diagnosis)) + 
  geom_violin(trim = FALSE)+ 
  geom_boxplot(width = 0.2)+
  labs(x = "",
       y = "LYVE1 ng/mg creatinine" , 
       title = "LYVE1") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45))


REG1Bvscreatinine <- clean_data %>%
  select(diagnosis , stage , REG1B, creatinine) %>%
  mutate(diagnosis = case_when(
    diagnosis == "control" ~ "control",
    diagnosis == "benign" ~ "benign",
    stage  == "I" ~ "PDAC I-II",
    stage  == "IA" ~ "PDAC I-II",
    stage  == "IB" ~ "PDAC I-II",
    stage  == "II" ~ "PDAC I-II",
    stage  == "IIA" ~ "PDAC I-II",
    stage  == "IIB" ~ "PDAC I-II",
    stage  == "III" ~ "PDAC III-IV",
    stage  == "IV" ~ "PDAC III-IV") , 
    REG1B = REG1B/ creatinine)

REG1Bvscreatinine_plot <- 
  ggplot(data = REG1Bvscreatinine, 
         mapping = aes(x = diagnosis, 
                       y = REG1B, 
                       color = diagnosis)) + 
  geom_violin(trim = FALSE) + 
  geom_boxplot(width = 0.2) +
  labs(x = "Diagnosis",
       y = "REG1B ng/mg creatinine" , 
       title = "REG1B") +
  theme_minimal(base_size = 10,
                base_family = "Avenir") +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45))

TFF1vscreatinine <- 
  select(clean_data, diagnosis , stage , TFF1, creatinine) %>%
  mutate(diagnosis = case_when(
    diagnosis == "control" ~ "control",
    diagnosis == "benign" ~ "benign",
    stage  == "I" ~ "PDAC I-II",
    stage  == "IA" ~ "PDAC I-II",
    stage  == "IB" ~ "PDAC I-II",
    stage  == "II" ~ "PDAC I-II",
    stage  == "IIA" ~ "PDAC I-II",
    stage  == "IIB" ~ "PDAC I-II",
    stage  == "III" ~ "PDAC III-IV",
    stage  == "IV" ~ "PDAC III-IV") ,
    TFF1 = TFF1/ creatinine)

TFF1vscreatinine_plot <-
  ggplot(data = TFF1vscreatinine, 
        mapping = aes(x = diagnosis, 
                      y = TFF1, 
                      color = diagnosis)) + 
  geom_violin(trim = FALSE)+ 
  geom_boxplot(width = 0.2)+
  labs(x = "",
       y = "TFF1 ng/mg creatinine" , 
       title = "TFF1") +
  theme_minimal(base_size = 10,
                base_family = "Avenir")+
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45))

violin_plot <- (LYVEvscreatinine_plot | REG1Bvscreatinine_plot | TFF1vscreatinine_plot)
violin_plot