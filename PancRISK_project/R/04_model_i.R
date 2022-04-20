library("broom")
library("cowplot")
library("scales")
library("nnet")
library("pROC")
library("caret")

my_data_clean <-
  my_data_clean %>% 
  mutate(diagnosis = factor(diagnosis)) %>%
  mutate(sex = factor(sex)) %>%
  mutate(diagnosis = relevel(diagnosis, ref = "control"))
  
my_data_plasma_CA19_9_NA <-
  my_data_clean %>% 
  filter(is.na(plasma_CA19_9))

my_data_REG1A_NA <-
  my_data_clean %>% 
  filter(is.na(REG1A))

my_data_plasma_zero <-
  my_data_clean %>% 
  filter(plasma_CA19_9 == 0 & !is.na(plasma_CA19_9))

my_data_REG1A_zero <-
  my_data_clean %>% 
  filter(REG1A == 0 & !is.na(REG1A))

my_data_clean <-
  my_data_clean %>% 
  filter(!is.na(plasma_CA19_9) & !is.na(REG1A) & REG1A != 0 & plasma_CA19_9 != 0) %>%
  mutate_at(c("REG1A", "REG1B", "LYVE1", "TFF1", "plasma_CA19_9"), log) %>%
  full_join(my_data_plasma_zero) %>%
  full_join(my_data_REG1A_zero) %>%
  mutate_at(c("REG1A", "REG1B", "LYVE1", "TFF1", "plasma_CA19_9"), scale, scale = FALSE) %>% 
  full_join(my_data_plasma_CA19_9_NA) %>%
  full_join(my_data_REG1A_NA) %>%
  mutate(REG1A = REG1A[,1],
         REG1B = REG1B[,1],
         LYVE1 = LYVE1[,1],
         TFF1 = TFF1[,1],
         plasma_CA19_9 = plasma_CA19_9[,1]) %>%
  group_by(diagnosis)

sampleGuide <- data_frame(
  diagnosis = c("control","benign","malignant"),
  amount = c(92,104,100)
)

my_data_train <-
  my_data_clean %>% 
  nest(-diagnosis) %>% 
  left_join(sampleGuide, by = "diagnosis") %>%
  mutate(Sample = map2(data, amount, sample_n))  %>% 
  unnest(Sample) %>%
  select(-data, -amount) %>%
  mutate(diagnosis = factor(diagnosis)) %>%
  mutate(diagnosis = relevel(diagnosis, ref = "control"))

my_data_test <-
  my_data_clean %>% 
  setdiff(my_data_train)

## PancRISK
form_risk = formula(diagnosis ~ REG1A + REG1B + LYVE1 + TFF1 + plasma_CA19_9 + age + creatinine)
pancrisk_model <- 
  my_data_train %>%
  multinom(formula = diagnosis ~ REG1A + REG1B + LYVE1 + TFF1 + plasma_CA19_9 + age + creatinine,
           data = .,
           model = TRUE)

pancrisk_simple_model <- 
  my_data_train %>%
  multinom(formula = diagnosis ~ REG1A + LYVE1 + TFF1 + age + creatinine,
           data = .,
           model = TRUE)

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
  multinom(formula = diagnosis ~ REG1A + LYVE1 + TFF1,
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

my_data_val <-
  my_data_test %>%
  ungroup() %>%
  select(diagnosis, REG1A, REG1B, LYVE1, TFF1, plasma_CA19_9, age, creatinine) %>%
  mutate(pred = as.numeric(predict(pancrisk_model, my_data_test))) %>%
  mutate(diagnosis = as.numeric(diagnosis)) %>%
  select(diagnosis, pred)

roc(my_data_val$diagnosis, my_data_val$pred,
              plot = TRUE, conf.level = 0.95, levels = levels(as.factor(c(1, 2)))) %>%
  ci.auc()

roc(my_data_val$diagnosis, my_data_val$pred,
    plot = TRUE, conf.level = 0.95, levels = levels(as.factor(c(1, 3)))) %>%
  ci.auc()

my_data_val_simple <-
  my_data_test %>%
  ungroup() %>%
  select(diagnosis, REG1A, LYVE1, TFF1) %>%
  mutate(pred = as.numeric(predict(pancrisk_simple_model, my_data_test))) %>%
  mutate(diagnosis = as.numeric(diagnosis)) %>%
  select(diagnosis, pred)

roc(my_data_val_simple$diagnosis, my_data_val_simple$pred,
    plot = TRUE, conf.level = 0.95, levels = levels(as.factor(c(1, 2)))) %>%
  ggroc(colour = 'steelblue', size = 2) +
  ggtitle(paste0('ROC Curve')) +
  theme_minimal()

roc(my_data_val_simple$diagnosis, my_data_val_simple$pred,
    plot = TRUE, conf.level = 0.95, levels = levels(as.factor(c(1, 3))))

conf_matrix <- table(my_data_val$pred, my_data_val$diagnosis)
confusionMatrix(conf_matrix)# %>% 
#  tidy()

conf_matrix_simple <- table(my_data_val_simple$pred, my_data_val_simple$diagnosis)
confusionMatrix(conf_matrix_simple)
