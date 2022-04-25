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

my_data_plasma <-
  my_data_train %>%
  filter(plasma_CA19_9 >= 0)

my_data_test <-
  my_data_clean %>% 
  setdiff(my_data_train)

pancrisk_model_plasma <- 
  my_data_train %>%
  multinom(formula = diagnosis ~ REG1B + LYVE1 + TFF1 + cutoff_plasma + age + creatinine,
           data = .,
           model = TRUE,
           na.action = na.omit)

pancrisk_model_1 <- fit(data = my_data_train, 
                        form = formula(diagnosis ~ REG1B + LYVE1 + TFF1 + age + creatinine), 
                        contrast_1 = "control", 
                        contrast_2 = "benign")

pancrisk_model_2 <- fit(data = my_data_train, 
                        form = formula(diagnosis ~ REG1B + LYVE1 + TFF1 + age + creatinine), 
                        contrast_1 = "control", 
                        contrast_2 = "malignant")

pancrisk_simple_model_1 <- fit(data = my_data_train, 
                               form = formula(diagnosis ~ REG1A + LYVE1 + TFF1 + age + creatinine), 
                               contrast_1 = "control", 
                               contrast_2 = "benign")

pancrisk_simple_model_2 <- fit(data = my_data_train, 
                               form = formula(diagnosis ~ REG1A + LYVE1 + TFF1 + age + creatinine), 
                               contrast_1 = "control", 
                               contrast_2 = "malignant")

pancrisk_just_plasma_model_1 <- fit(data = my_data_plasma, 
                                    form = formula(diagnosis ~ cutoff_plasma + age), 
                                    contrast_1 = "control", 
                                    contrast_2 = "benign")

pancrisk_just_plasma_model_2 <- fit(data = my_data_plasma, 
                                    form = formula(diagnosis ~ cutoff_plasma + age), 
                                    contrast_1 = "control", 
                                    contrast_2 = "malignant")

pancrisk_logistic_1 <- fit(data = my_data_plasma, 
                           form = formula(diagnosis ~ REG1B + LYVE1 + TFF1 + cutoff_plasma + age + creatinine), 
                           contrast_1 = "control", 
                           contrast_2 = "benign")

pancrisk_logistic_2 <- fit(data = my_data_plasma, 
                           form = formula(diagnosis ~ REG1B + LYVE1 + TFF1 + cutoff_plasma + age + creatinine), 
                           contrast_1 = "control", 
                           contrast_2 = "malignant")

my_data_val_just_plasma_1 <-
  my_data_plasma %>%
  ungroup() %>%
  select(diagnosis, cutoff_plasma, age) %>%
  filter(diagnosis == "control" | diagnosis == "malignant") %>%
  drop_na() %>%
  mutate(pred = as.numeric(predict(pancrisk_just_plasma_model_1, .))) %>%
  mutate(diagnosis = as.numeric(diagnosis)) %>%
  select(diagnosis, pred)

my_data_val_just_plasma_2 <-
  my_data_plasma %>%
  ungroup() %>%
  select(diagnosis, cutoff_plasma, age) %>%
  filter(diagnosis == "benign" | diagnosis == "malignant") %>%
  drop_na() %>%
  mutate(pred = as.numeric(predict(pancrisk_just_plasma_model_2, .))) %>%
  mutate(diagnosis = as.numeric(diagnosis)) %>%
  select(diagnosis, pred)

my_data_val_simple_1 <-
  my_data_test %>%
  ungroup() %>%
  select(diagnosis, REG1A, LYVE1, TFF1, age, creatinine) %>%
  filter(diagnosis == "control" | diagnosis == "malignant") %>%
  drop_na() %>%
  mutate(pred = as.numeric(predict(pancrisk_simple_model_1, .))) %>%
  mutate(diagnosis = as.numeric(diagnosis)) %>%
  select(diagnosis, pred)

my_data_val_simple_2 <-
  my_data_test %>%
  ungroup() %>%
  select(diagnosis, REG1A, LYVE1, TFF1, age, creatinine) %>%
  filter(diagnosis == "benign" | diagnosis == "malignant") %>%
  drop_na() %>%
  mutate(pred = as.numeric(predict(pancrisk_simple_model_2, .))) %>%
  mutate(diagnosis = as.numeric(diagnosis)) %>%
  select(diagnosis, pred)

my_data_val_1 <-
  my_data_test %>%
  ungroup() %>%
  select(diagnosis, REG1B, LYVE1, TFF1, age, creatinine) %>%
  filter(diagnosis == "control" | diagnosis == "malignant") %>%
  drop_na() %>%
  mutate(pred = as.numeric(predict(pancrisk_model_2, .))) %>%
  mutate(diagnosis = as.numeric(diagnosis)) %>%
  select(diagnosis, pred)

my_data_val_2 <-
  my_data_test %>%
  ungroup() %>%
  select(diagnosis, REG1B, LYVE1, TFF1, age, creatinine) %>%
  filter(diagnosis == "benign" | diagnosis == "malignant") %>%
  drop_na() %>%
  mutate(pred = as.numeric(predict(pancrisk_model_2, .))) %>%
  mutate(diagnosis = as.numeric(diagnosis)) %>%
  select(diagnosis, pred)

my_data_val_plasma <-
  my_data_test %>%
  ungroup() %>%
  select(diagnosis, REG1B, LYVE1, TFF1, cutoff_plasma, age, creatinine) %>%
  drop_na() %>%
  mutate(pred = as.numeric(predict(pancrisk_model_plasma, .))) %>%
  mutate(diagnosis = as.numeric(diagnosis)) %>%
  select(diagnosis, pred)

my_data_val_log_1 <-
  my_data_test %>%
  ungroup() %>%
  filter(diagnosis == "control" | diagnosis == "malignant") %>%
  select(diagnosis, REG1B, LYVE1, TFF1, cutoff_plasma, age, creatinine) %>%
  drop_na() %>%
  mutate(pred = as.numeric(predict(pancrisk_logistic_1, .))) %>%
  mutate(diagnosis = as.numeric(diagnosis)) %>%
  select(diagnosis, pred)

my_data_val_log_2 <-
  my_data_test %>%
  ungroup() %>%
  filter(diagnosis == "benign" | diagnosis == "malignant") %>%
  select(diagnosis, REG1B, LYVE1, TFF1, cutoff_plasma, age, creatinine) %>%
  drop_na() %>%
  mutate(pred = as.numeric(predict(pancrisk_logistic_2, .))) %>%
  mutate(diagnosis = as.numeric(diagnosis)) %>%
  select(diagnosis, pred)

roc1 <- my_data_val_just_plasma_1 %>%
        pROC::roc(diagnosis, pred,
            plot = FALSE, conf.level = 0.95)

roc2 <- my_data_val_just_plasma_2 %>%
        pROC::roc(diagnosis, pred,
            plot = FALSE, conf.level = 0.95) 

roc3 <- my_data_val_simple_1 %>%
        pROC::roc(diagnosis, pred,
            plot = FALSE, conf.level = 0.95)

roc4 <- my_data_val_simple_2 %>%
        pROC::roc(diagnosis, pred,
            plot = FALSE, conf.level = 0.95) 

roc5 <- my_data_val_1 %>%
        pROC::roc(diagnosis, pred,
            plot = FALSE, conf.level = 0.95)

roc6 <- my_data_val_2 %>%
        pROC::roc(diagnosis, pred,
            plot = FALSE, conf.level = 0.95) 

roc7 <- my_data_val_plasma %>%
        pROC::roc(diagnosis, pred,
            plot = FALSE, conf.level = 0.95, levels = c(1, 2))

roc8 <- my_data_val_plasma %>%
        pROC::roc(diagnosis, pred,
            plot = FALSE, conf.level = 0.95, levels = c(2, 3)) 

roc9 <- my_data_val_log_1 %>%
        pROC::roc(diagnosis, pred,
            plot = FALSE, conf.level = 0.95)

roc10 <- my_data_val_log_2 %>%
         pROC::roc(diagnosis, pred,
            plot = FALSE, conf.level = 0.95) 

my_data_val_just_plasma_1 <- find_cutoff(data = my_data_val_just_plasma_1, 
                                         const_1 = 1, 
                                         const_2 = 3)

my_data_val_just_plasma_2 <- find_cutoff(data = my_data_val_just_plasma_2, 
                                         const_1 = 2, 
                                         const_2 = 3)

my_data_val_simple_1 <- find_cutoff(data = my_data_val_simple_1, 
                                    const_1 = 1, 
                                    const_2 = 3)

my_data_val_simple_2 <- find_cutoff(data = my_data_val_simple_2, 
                                    const_1 = 2, 
                                    const_2 = 3)

my_data_val_1 <- find_cutoff(data = my_data_val_1, 
                             const_1 = 1, 
                             const_2 = 3)

my_data_val_2 <- find_cutoff(data = my_data_val_2, 
                             const_1 = 2, 
                             const_2 = 3)

my_data_val_log_1 <- find_cutoff(data = my_data_val_log_1, 
                                 const_1 = 1, 
                                 const_2 = 3)

my_data_val_log_2 <- find_cutoff(data = my_data_val_log_2, 
                                 const_1 = 2, 
                                 const_2 = 3)

acc_1 <- 
  print_acc(data_1 = table(my_data_val_just_plasma_1$pred, my_data_val_just_plasma_1$diagnosis) %>%
              confusionMatrix(),
            data_2 = table(my_data_val_just_plasma_2$pred, my_data_val_just_plasma_2$diagnosis) %>%
              confusionMatrix(),
            title = "Binary CA19_9")

acc_2 <- 
  print_acc(data_1 = table(my_data_val_simple_1$pred, my_data_val_simple_1$diagnosis) %>%
              confusionMatrix(),
            data_2 = table(my_data_val_simple_2$pred, my_data_val_simple_2$diagnosis) %>%
              confusionMatrix(),
            title = "Old panel")

acc_3 <-
  print_acc(data_1 = table(my_data_val_1$pred, my_data_val_1$diagnosis) %>%
              confusionMatrix(),
            data_2 = table(my_data_val_2$pred, my_data_val_2$diagnosis) %>%
              confusionMatrix(),
            title = "PancRISK without Plasma")

acc_4 <- 
  print_acc(data_1 = table(my_data_val_log_1$pred, my_data_val_log_1$diagnosis) %>%
      confusionMatrix(),
      data_2 = table(my_data_val_log_2$pred, my_data_val_log_2$diagnosis) %>%
      confusionMatrix(),
    title = "PancRISK (Sep. regression)")