library("broom")
library("cowplot")
library("scales")
library("nnet")
library("pROC")

my_data_clean <-
  my_data_clean %>% 
  mutate(diagnosis = factor(diagnosis)) %>%
  mutate(sex = factor(sex)) %>%
  mutate(diagnosis = relevel(diagnosis, ref = "control")) %>%
  drop_na(REG1A, REG1B, LYVE1, TFF1, plasma_CA19_9) %>%
  filter(plasma_CA19_9 != 0 & REG1A != 0) %>%
  mutate_at(c("REG1A", "REG1B", "LYVE1", "TFF1", "plasma_CA19_9"), log) %>%
  mutate_at(c("REG1A", "REG1B", "LYVE1", "TFF1", "plasma_CA19_9"), scale, scale = FALSE) %>% 
  mutate(REG1A = REG1A[,1],
         REG1B = REG1B[,1],
         LYVE1 = LYVE1[,1],
         TFF1 = TFF1[,1],
         plasma_CA19_9 = plasma_CA19_9[,1])

my_data_train <-
  my_data_clean %>% 
  sample_n(85)

my_data_test <-
  my_data_clean %>% 
  setdiff(my_data_train)

## PancRISK
form_risk = formula(diagnosis ~ REG1A + REG1B + LYVE1 + TFF1 + plasma_CA19_9 + age + creatinine)
pancrisk_model <- 
  my_data_train %>%
  multinom(formula = form_risk,
           data = .,
           model = TRUE,
           na.action = na.omit)

pancrisk_model %>%
  tidy() %>%
  select(y.level, term, p.value) %>%
  group_by(y.level) %>%
  mutate(significant = case_when(p.value < 0.05 ~ TRUE,
                                 p.value >= 0.05 ~ FALSE))

my_data_test %>%
  ungroup() %>%
  select(diagnosis, REG1A, REG1B, LYVE1, TFF1, plasma_CA19_9, age, creatinine) %>%
  mutate(prediction = as.numeric(predict(pancrisk_model, my_data_test))) %>%
  mutate(diagnosis = as.numeric(diagnosis)) %>%
  select(diagnosis, prediction) %>%
  drop_na() %>%
  cor()

my_data_train %>%
  multinom(formula = diagnosis ~ plasma_CA19_9,
           data = .) %>%
  tidy() %>%
  select(y.level, term, p.value) %>%
  group_by(y.level) %>%
  mutate(significant = case_when(p.value < 0.05 ~ TRUE,
                                 p.value >= 0.05 ~ FALSE))

my_data_train %>%
  mutate(stage = factor(stage)) %>%
  multinom(formula = stage ~ REG1A + REG1B + LYVE1 + TFF1 + plasma_CA19_9 + age + creatinine,
           data = .) %>%
  tidy() %>%
  select(y.level, term, p.value) %>%
  group_by(y.level) %>%
  mutate(significant = case_when(p.value < 0.05 ~ TRUE,
                                 p.value >= 0.05 ~ FALSE))

my_data_test <-
  my_data_test %>%
  ungroup() %>%
  select(diagnosis, REG1A, REG1B, LYVE1, TFF1, plasma_CA19_9, age, creatinine) %>%
  mutate(pred = as.numeric(predict(pancrisk_model, my_data_test))) %>%
  mutate(diagnosis = as.numeric(diagnosis)) %>%
  select(diagnosis, pred) %>%
  drop_na()
multiclass.roc(my_data_test$diagnosis, my_data_test$pred,
               plot = TRUE)