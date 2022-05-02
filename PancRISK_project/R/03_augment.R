## Modelling data augment

my_data_cleaner <-
  my_data_clean %>% 
  mutate(diagnosis = factor(diagnosis)) %>%
  mutate(sex = factor(sex)) %>%
  mutate(diagnosis = relevel(diagnosis, ref = "control"))

my_data_REG1A_NA <-
  my_data_cleaner %>% 
  filter(is.na(REG1A))

my_data_REG1A_zero <-
  my_data_cleaner %>% 
  filter(REG1A == 0 & !is.na(REG1A))

my_data_cleaner <-
  my_data_cleaner %>% 
  filter(!is.na(plasma_CA19_9) & !is.na(REG1A) & REG1A != 0 & plasma_CA19_9 != 0) %>%
  mutate_at(c("REG1A", "REG1B", "LYVE1", "TFF1"), log) %>%
  full_join(my_data_REG1A_zero) %>%
  mutate_at(c("REG1A", "REG1B", "LYVE1", "TFF1"), scale, scale = FALSE) %>% 
  full_join(my_data_REG1A_NA) %>%
  mutate(REG1A = REG1A[,1],
         REG1B = REG1B[,1],
         LYVE1 = LYVE1[,1],
         TFF1 = TFF1[,1]) %>%
  mutate(cutoff_plasma = ifelse(plasma_CA19_9 > 37, 1, 0)) %>% 
  group_by(diagnosis)

write_tsv(x = my_data_cleaner,
          file = "PancRISK_project/data/03_my_data_cleaner.tsv")