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

print_acc <- function(data, title = "Group") {
  accurancy <- 
    data %>%
    tidy() %>%
    filter(term == "accuracy") %>%
    mutate(result = paste(round(estimate, 3), " 95% CI: ", round(conf.low, 3), "-", round(conf.high, 3))) %>%
    select(result) %>%
    pull()
  data %>%
    tidy() %>%
    pivot_wider(values_from = estimate, names_from = class) %>%
    filter(term == "sensitivity" | term == "specificity") %>%
    select(c(term, '1', '2', '3')) %>%
    mutate_at(c('1', '2', '3'), round, 3) %>%
    mutate(term = c("Sensitivity", "Specificity")) %>%
    rename_at(c("term", "1", "2", "3"), ~ (c(" ", "Control", "Benign", "Malignant"))) %>% 
    knitr::kable() %>%
    add_header_above(c("Accurancy" = 1, setNames(3, accurancy)), line = FALSE, bold = FALSE) %>%
    add_header_above(c(setNames(4, paste("Performance - ", title)))) %>% 
    kable_styling()
}

print_p <- function(data) {
  accurancy <- 
    data %>%
    tidy() %>%
    filter(term == "accuracy") %>%
    mutate(result = paste("AUC: ", round(estimate, 3), " 95% CI: ", round(conf.low, 3), "-", round(conf.high, 3))) %>%
    select(result) %>%
    pull()
  data %>%
    tidy() %>%
    pivot_wider(values_from = estimate, names_from = class) %>%
    filter(term == "sensitivity" | term == "specificity") %>%
    select(c(term, '1', '2', '3')) %>%
    mutate_at(c('1', '2', '3'), round, 3) %>%
    mutate(term = c("Sensitivity", "Specificity")) %>%
    rename_at(c("term", "1", "2", "3"), ~ (c("Group", "Control", "Benign", "Malignant"))) %>% 
    knitr::kable() %>%
    add_header_above(c(setNames(4, accurancy)), line = FALSE, bold = FALSE) %>%
    add_header_above(c("Performance" = 4)) %>%
    kable_styling()
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
                                  paste0("Control vs. Malignant\n",
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

pancrisk_just_plasma_model <- 
  my_data_train %>%
  multinom(formula = diagnosis ~ cutoff_plasma,
           data = .)

# pancrisk_model %>%
#   tidy() %>%
#   select(y.level, term, p.value) %>%
#   group_by(y.level) %>%
#   mutate(significant = case_when(p.value < 0.05 ~ TRUE,
#                                  p.value >= 0.05 ~ FALSE))

## Only plasma

my_data_val_just_plasma <-
  my_data_test %>%
  ungroup() %>%
  select(diagnosis, cutoff_plasma) %>%
  drop_na() %>%
  mutate(pred = as.numeric(predict(pancrisk_just_plasma_model, .))) %>%
  mutate(diagnosis = as.numeric(diagnosis)) %>%
  select(diagnosis, pred)

my_data_val_simple <-
  my_data_test %>%
  ungroup() %>%
  select(diagnosis, REG1A, LYVE1, TFF1, age, creatinine) %>%
  drop_na() %>%
  mutate(pred = as.numeric(predict(pancrisk_simple_model, .))) %>%
  mutate(diagnosis = as.numeric(diagnosis)) %>%
  select(diagnosis, pred)

my_data_val <-
  my_data_test %>%
  ungroup() %>%
  select(diagnosis, REG1B, LYVE1, TFF1, age, creatinine) %>%
  drop_na() %>%
  mutate(pred = as.numeric(predict(pancrisk_model, .))) %>%
  mutate(diagnosis = as.numeric(diagnosis)) %>%
  select(diagnosis, pred)

my_data_val_plasma <-
  my_data_test %>%
  ungroup() %>%
  select(diagnosis, REG1B, LYVE1, TFF1, cutoff_plasma, age, creatinine) %>%
  mutate(pred = as.numeric(predict(pancrisk_model_plasma, .))) %>%
  mutate(diagnosis = as.numeric(diagnosis)) %>%
  select(diagnosis, pred)

roc1 <- roc(my_data_val_just_plasma$diagnosis, my_data_val_just_plasma$pred,
            plot = FALSE, conf.level = 0.95, levels = levels(as.factor(c(1, 2))))

roc2 <- roc(my_data_val_just_plasma$diagnosis, my_data_val_just_plasma$pred,
            plot = FALSE, conf.level = 0.95, levels = levels(as.factor(c(1, 3)))) 

roc3 <- roc(my_data_val_simple$diagnosis, my_data_val_simple$pred,
            plot = FALSE, conf.level = 0.95, levels = levels(as.factor(c(1, 2))))

roc4 <- roc(my_data_val_simple$diagnosis, my_data_val_simple$pred,
            plot = FALSE, conf.level = 0.95, levels = levels(as.factor(c(1, 3)))) 

roc5 <- roc(my_data_val$diagnosis, my_data_val$pred,
            plot = FALSE, conf.level = 0.95, levels = levels(as.factor(c(1, 2))))

roc6 <- roc(my_data_val$diagnosis, my_data_val$pred,
            plot = FALSE, conf.level = 0.95, levels = levels(as.factor(c(1, 3)))) 

roc7 <- roc(my_data_val_plasma$diagnosis, my_data_val_plasma$pred,
            plot = FALSE, conf.level = 0.95, levels = levels(as.factor(c(1, 2))))

roc8 <- roc(my_data_val_plasma$diagnosis, my_data_val_plasma$pred,
            plot = FALSE, conf.level = 0.95, levels = levels(as.factor(c(1, 3)))) 

# acc_1 <- 
#   table(my_data_val_just_plasma$pred, my_data_val_just_plasma$diagnosis) %>%
#   confusionMatrix() %>%
#   print_acc()

acc_2 <-   
  table(my_data_val_simple$pred, my_data_val_simple$diagnosis) %>%
  confusionMatrix() %>%
  print_acc(title = "Simple")
  
acc_3 <- 
  table(my_data_val$pred, my_data_val$diagnosis) %>%
  confusionMatrix() %>%
  print_acc(title = "PancRISK without Plasma")
  
acc_4 <- 
  table(my_data_val_plasma$pred, my_data_val_plasma$diagnosis) %>%
  confusionMatrix() %>%
  print_acc(title = "PancRISK complete")