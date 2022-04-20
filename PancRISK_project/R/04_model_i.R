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

my_data_REG1A_NA <-
  my_data_clean %>% 
  filter(is.na(REG1A))

my_data_REG1A_zero <-
  my_data_clean %>% 
  filter(REG1A == 0 & !is.na(REG1A))

my_data_clean <-
  my_data_clean %>% 
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

sampleGuide <- data_frame(
  diagnosis = c("control", "benign", "malignant"),
  amount = c(92, 104, 100) #Half of each group
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

## Plasma


## PancRISK
form_risk = formula(diagnosis ~  REG1B + LYVE1 + TFF1 + cutoff_plasma + age + creatinine)
pancrisk_model_plasma <- 
  my_data_train %>%
  multinom(formula = diagnosis ~ REG1B + LYVE1 + TFF1 + cutoff_plasma + age + creatinine,
           data = .,
           model = TRUE,
           na.action = na.omit)

pancrisk_model <- 
  my_data_train %>%
  multinom(formula = diagnosis ~ REG1B + LYVE1 + TFF1 + age + creatinine,
           data = .,
           model = TRUE,
           na.action = na.omit)

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

## Only plasma
my_data_train %>%
  multinom(formula = diagnosis ~ cutoff_plasma,
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

my_data_val <-
  my_data_test %>%
  ungroup() %>%
  select(diagnosis, REG1A, REG1B, LYVE1, TFF1, cutoff_plasma, age, creatinine) %>%
  mutate(pred = as.numeric(predict(pancrisk_model_plasma, my_data_test))) %>%
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

my_data_val_train <-
  my_data_train %>%
  ungroup() %>%
  select(diagnosis, REG1A, REG1B, LYVE1, TFF1, cutoff_plasma, age, creatinine) %>%
  mutate(pred = as.numeric(predict(pancrisk_model, my_data_train))) %>%
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

plotx <- rev(roc_test$specificities)
ploty <- rev(roc_test$sensitivities)
ggplot(NULL, aes(x = plotx, y = ploty)) +
  geom_segment(aes(x = 0, y = 1, xend = 1, yend = 0), alpha = 0.5) + 
  geom_step() +
  scale_x_reverse(name = "Specificity",limits = c(1,0), breaks = seq(0, 1, 0.2), expand = c(0.001,0.001)) + 
  scale_y_continuous(name = "Sensitivity", limits = c(0,1), breaks = seq(0, 1, 0.2), expand = c(0.001, 0.001)) +
  theme_minimal() + 
  coord_equal() + 
  ggtitle("text", x = 0.2/2, y = 0.2/2, vjust = 0, label = paste("AUC =",sprintf("%.3f",roc_test$auc)))
