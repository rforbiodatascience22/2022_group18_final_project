library("broom")
library("cowplot")
library("scales")
library("nnet")
library("pROC")
library("caret")
library("broom")
library("knitr")
library("patchwork")
library("kableExtra")
library("cutoff")

find_cutoff <- function(data, const_1, const_2) {
  cutoff_val <- cutoff::roc(data$pred, data$diagnosis)$cutoff
  
  data <-
    data %>%
    mutate(pred = ifelse(pred < cutoff_val, const_1, const_2))
}

fit <- function(data, form, contrast) {
  data %>%
  filter(diagnosis == contrast[1] | diagnosis == contrast[2]) %>%
  glm(formula = form,
      data = .,
      na.action = na.omit,
      family = binomial(link = "logit"))
}

print_acc <- function(data_1, data_2, title = "Group") {
  accurancy_1 <- 
    data_1 %>%
    tidy() %>%
    filter(term == "accuracy") %>%
    mutate(result = paste(round(estimate, 3), " 95% CI: ", round(conf.low, 3), "-", round(conf.high, 3))) %>%
    select(result) %>%
    pull()
  accurancy_2 <- 
    data_1 %>%
    tidy() %>%
    filter(term == "accuracy") %>%
    mutate(result = paste(round(estimate, 3), " 95% CI: ", round(conf.low, 3), "-", round(conf.high, 3))) %>%
    select(result) %>%
    pull()
  table_1 <-
    data_1 %>%
    tidy() %>%
    pivot_wider(values_from = estimate, 
                names_from = class) %>%
    filter(term == "sensitivity" | term == "specificity") %>%
    select(c(term, '1')) %>%
    mutate_at(c('1'), round, 3) %>%
    mutate(term = c("Sensitivity", "Specificity")) 

  data_2 %>%
    tidy() %>%
    pivot_wider(values_from = estimate, names_from = class) %>%
    filter(term == "sensitivity" | term == "specificity") %>%
    select(c(term, '2')) %>%
    mutate_at(c('2'), round, 3) %>%
    mutate(term = c("Sensitivity", "Specificity")) %>%
    full_join(table_1) %>%
    knitr::kable(align = c(rep('c', times = 3)), 
                 col.names = NULL) %>%
    #add_header_above(c("Benign vs. Malignant" = 1, " " = 1)) %>% 
    add_header_above(c("a\na" = 1, setNames(1, "Control vs. Malignant"), setNames(1, "Benign vs. Malignant")), line = TRUE, bold = FALSE) %>%
    add_header_above(c("Accurancy\na" = 1, setNames(1, accurancy_1), setNames(1, accurancy_2)), line = TRUE, bold = FALSE) %>%
    add_header_above(c(setNames(3, paste("Performance - ", title)))) %>% 
    kable_styling(position = "center") %>%
    column_spec (1:3, border_left = "1px solid #ddd;", border_right = "1px solid #ddd;") %>%
    row_spec(row = 1:2, font_size = 24, hline_after = "1px solid #ddd;") %>%
    column_spec(column = 1, width = "6em")
}

plot_roc <- function(roc_render_1, roc_render_2, title = "ROC") {
  plotx_1 <- rev(roc_render_1$specificities)
  ploty_1 <- rev(roc_render_1$sensitivities)
  auc_1 <- roc_render_1$auc
  auc_ci_1 <- roc_render_1 %>%
    ci.auc() %>%
    round(3)
  plotx_2 <- rev(roc_render_2$specificities)
  ploty_2 <- rev(roc_render_2$sensitivities)
  auc_2 <- roc_render_2$auc
  auc_ci_2 <- roc_render_2 %>%
    ci.auc() %>%
    round(3)
  ggplot(NULL) +
    geom_segment(aes(x = 0, y = 1, xend = 1, yend = 0), alpha = 0.5) + 
    geom_line(aes(x = plotx_1, y = ploty_1, color = "low")) +
    geom_line(aes(x = plotx_2, y = ploty_2, color = "high")) +
    scale_color_discrete(name="Contrast",
                         labels=c(paste0("Control vs. Benign\n",
                                         round(auc_1, 3), ", 95% CI: ", auc_ci_1[1], "-", auc_ci_1[3]), 
                                  paste0("Benign vs. Malignant\n",
                                         round(auc_2, 3), ", 95% CI: ", auc_ci_2[1], "-", auc_ci_2[3])),
                         breaks=c("low", "high")) +
    scale_x_reverse(name = "Specificity",
                    limits = c(1, 0), 
                    breaks = seq(0, 1, 0.2), 
                    expand = c(0.001,0.001)) + 
    scale_y_continuous(name = "Sensitivity", 
                       limits = c(0, 1), 
                       breaks = seq(0, 1, 0.2), 
                       expand = c(0.001, 0.001)) +
    theme_minimal() + 
    theme(plot.margin = unit(c(10, 10, 10, 10), "pt")) +
    coord_equal() +
    ggtitle(str_c(title))
} 

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

## Plasma


## PancRISK
form_risk = formula(diagnosis ~  REG1B + LYVE1 + TFF1 + cutoff_plasma + age + creatinine)
pancrisk_model_plasma <- 
  my_data_train %>%
  multinom(formula = diagnosis ~ REG1B + LYVE1 + TFF1 + cutoff_plasma + age + creatinine,
           data = .,
           model = TRUE,
           na.action = na.omit)

pancrisk_model_1 <- fit(my_data_train, 
                        formula(diagnosis ~ REG1B + LYVE1 + TFF1 + age + creatinine), 
                        contrast = c("control", "benign"))

pancrisk_model_2 <- fit(my_data_train, 
                        formula(diagnosis ~ REG1B + LYVE1 + TFF1 + age + creatinine), 
                        contrast = c("control", "malignant"))

pancrisk_simple_model_1 <- fit(my_data_train, 
                               formula(diagnosis ~ REG1A + LYVE1 + TFF1 + age + creatinine), 
                               contrast = c("control", "benign"))

pancrisk_simple_model_2 <- fit(my_data_train, 
                               formula(diagnosis ~ REG1A + LYVE1 + TFF1 + age + creatinine), 
                               contrast = c("control", "malignant"))

pancrisk_just_plasma_model_1 <- fit(my_data_plasma, 
                                    formula(diagnosis ~ cutoff_plasma + age), 
                                    contrast = c("control", "benign"))

pancrisk_just_plasma_model_2 <- fit(my_data_plasma, 
                                    formula(diagnosis ~ cutoff_plasma + age), 
                                    contrast = c("control", "malignant"))

pancrisk_logistic_1 <- fit(my_data_plasma, 
                           formula(diagnosis ~ REG1B + LYVE1 + TFF1 + cutoff_plasma + age + creatinine), 
                           contrast = c("control", "benign"))

pancrisk_logistic_2 <- fit(my_data_plasma, 
                           formula(diagnosis ~ REG1B + LYVE1 + TFF1 + cutoff_plasma + age + creatinine), 
                           contrast = c("control", "malignant"))

# pancrisk_model %>%
#   tidy() %>%
#   select(y.level, term, p.value) %>%
#   group_by(y.level) %>%
#   mutate(significant = case_when(p.value < 0.05 ~ TRUE,
#                                  p.value >= 0.05 ~ FALSE))

## Only plasma

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

roc1 <- pROC::roc(my_data_val_just_plasma_1$diagnosis, my_data_val_just_plasma_1$pred,
            plot = FALSE, conf.level = 0.95)

roc2 <- pROC::roc(my_data_val_just_plasma_2$diagnosis, my_data_val_just_plasma_2$pred,
            plot = FALSE, conf.level = 0.95) 

roc3 <- pROC::roc(my_data_val_simple_1$diagnosis, my_data_val_simple_1$pred,
            plot = FALSE, conf.level = 0.95)

roc4 <- pROC::roc(my_data_val_simple_2$diagnosis, my_data_val_simple_2$pred,
            plot = FALSE, conf.level = 0.95) 

roc5 <- pROC::roc(my_data_val_1$diagnosis, my_data_val_1$pred,
            plot = FALSE, conf.level = 0.95)

roc6 <- pROC::roc(my_data_val_2$diagnosis, my_data_val_2$pred,
            plot = FALSE, conf.level = 0.95) 

roc7 <- pROC::roc(my_data_val_plasma$diagnosis, my_data_val_plasma$pred,
            plot = FALSE, conf.level = 0.95, levels = c(1, 2))

roc8 <- pROC::roc(my_data_val_plasma$diagnosis, my_data_val_plasma$pred,
            plot = FALSE, conf.level = 0.95, levels = c(2, 3)) 

roc9 <- pROC::roc(my_data_val_log_1$diagnosis, my_data_val_log_1$pred,
            plot = FALSE, conf.level = 0.95)

roc10 <- pROC::roc(my_data_val_log_2$diagnosis, my_data_val_log_2$pred,
            plot = FALSE, conf.level = 0.95) 

my_data_val_just_plasma_1 <- find_cutoff(my_data_val_just_plasma_1, 1, 3)

my_data_val_just_plasma_2 <- find_cutoff(my_data_val_just_plasma_2, 2, 3)

my_data_val_simple_1 <- find_cutoff(my_data_val_simple_1, 1, 3)

my_data_val_simple_2 <- find_cutoff(my_data_val_simple_2, 2, 3)

my_data_val_1 <- find_cutoff(my_data_val_1, 1, 3)

my_data_val_2 <- find_cutoff(my_data_val_2, 2, 3)

my_data_val_log_1 <- find_cutoff(my_data_val_log_1, 1, 3)

my_data_val_log_2 <- find_cutoff(my_data_val_log_2, 2, 3)

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

# acc_5 <- 
#   table(my_data_val_log_2$pred, my_data_val_log_2$diagnosis) %>%
#   confusionMatrix() %>%
#   rename_at(c(2, 3), ~ (c(1, 2))) %>% 
#   print_acc(title = "PancRISK - seperate logistic regression")