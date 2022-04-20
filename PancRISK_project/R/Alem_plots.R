library("tidyverse")
data_long <- my_data_clean %>% 
  select(patient_cohort, plasma_CA19_9, creatinine, LYVE1, REG1A, REG1B, TFF1) %>% 
  pivot_longer(!patient_cohort, names_to = "protein", values_to = "expression")

data_long %>%
  group_by(patient_cohort) %>%
  ggplot(mapping = aes(y = patient_cohort, fill = patient_cohort)) +
  facet_wrap(~protein, ncol = 2, scales = "free") +
  geom_bar() +
  xlab("Occurence of measured protein levels per Cohort ") +
  ylab("Cohort of patients") +
  theme_minimal() +
  theme(legend.position="none")
 