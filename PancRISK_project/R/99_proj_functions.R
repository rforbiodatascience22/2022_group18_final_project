# Validate ---------------------------------------------------------------------
validate <- function(data, 
                     model, 
                     diagnosis_1, 
                     diagnosis_2) {
  
}
# Cutoff -----------------------------------------------------------------------
find_cutoff <- function(data, 
                        const_1, 
                        const_2) {
  cutoff_val <- cutoff::roc(data$pred, 
                            data$diagnosis)$cutoff
  
  data %>%
    mutate(pred = ifelse(pred < cutoff_val, 
                         const_1, 
                         const_2))
}
# Fit --------------------------------------------------------------------------
fit <- function(data, 
                form, 
                contrast_1, 
                contrast_2) {
  data %>%
    filter(diagnosis == contrast_1 | diagnosis == contrast_2) %>%
    glm(formula = form,
        data = .,
        na.action = na.omit,
        family = binomial(link = "logit"))
}
# Accuracy ---------------------------------------------------------------------
print_acc <- function(data_1, 
                      data_2, 
                      title = "Group") {
  accurancy_1 <- 
    data_1 %>%
    tidy() %>%
    filter(term == "accuracy") %>%
    mutate(result = paste(round(estimate, 
                                digits = 3), 
                          " 95% CI: ", 
                          round(conf.low, 
                                digits = 3), 
                          "-", 
                          round(conf.high,
                                digits = 3))) %>%
    select(result) %>%
    pull()
  
  accurancy_2 <- 
    data_2 %>%
    tidy() %>%
    filter(term == "accuracy") %>%
    mutate(result = paste(round(estimate, 
                                digits = 3), 
                          " 95% CI: ", 
                          round(conf.low, 
                                digits = 3), 
                          "-", 
                          round(conf.high, 
                                digits = 3))) %>%
    select(result) %>%
    pull()
  
  table_1 <-
    data_1 %>%
    tidy() %>%
    pivot_wider(values_from = estimate, 
                names_from = class) %>%
    filter(term == "sensitivity" | term == "specificity") %>%
    select(c(term, 
             '1')) %>%
    mutate_at(c('1'), round, 3) %>%
    mutate(term = c("Sensitivity", 
                    "Specificity")) 
  
  data_2 %>%
    tidy() %>%
    pivot_wider(values_from = estimate, 
                names_from = class) %>%
    filter(term == "sensitivity" | term == "specificity") %>%
    select(c(term, 
             '2')) %>%
    mutate_at(c('2'), 
              round, 
              3) %>%
    mutate(term = c("Sensitivity", 
                    "Specificity")) %>%
    full_join(table_1) %>%
    knitr::kable(align = c(rep('c', 
                               times = 3)), 
                 col.names = NULL,) %>%
    #add_header_above(c("Benign vs. Malignant" = 1, " " = 1)) %>% 
    add_header_above(c("\ \n\n \ " = 1, 
                       setNames(object = 1, 
                                nm = "Control vs. Malignant"), 
                       setNames(object = 1, 
                                nm = "Benign vs. Malignant")), 
                     line = TRUE, 
                     bold = FALSE) %>%
    add_header_above(c("Accurancy\n\n \ " = 1, 
                       setNames(object = 1, 
                                nm = accurancy_1), 
                       setNames(object = 1, 
                                nm = accurancy_2)), 
                     line = TRUE, 
                     bold = FALSE) %>%
    add_header_above(c(setNames(object = 3, 
                                nm = paste("Performance - ", 
                                           title)))) %>% 
    kable_styling(position = "center") %>%
    column_spec (1:3, 
                 border_left = "1px solid #ddd;", 
                 border_right = "1px solid #ddd;") %>%
    row_spec(row = 1:2, 
             font_size = 20, 
             hline_after = "1px solid #ddd;") %>%
    column_spec(column = 1, 
                width = "6em")
}

# ROC --------------------------------------------------------------------------
plot_roc <- function(roc_render_1, 
                     roc_render_2, 
                     title = "ROC", 
                     caption) {
  plotx_1 <- rev(roc_render_1$specificities)
  ploty_1 <- rev(roc_render_1$sensitivities)
  auc_1 <- roc_render_1$auc
  auc_ci_1 <- roc_render_1 %>%
    ci.auc() %>%
    round(digits = 3)
  plotx_2 <- rev(roc_render_2$specificities)
  ploty_2 <- rev(roc_render_2$sensitivities)
  auc_2 <- roc_render_2$auc
  auc_ci_2 <- roc_render_2 %>%
    ci.auc() %>%
    round(digits = 3)
  ggplot(NULL) +
    geom_segment(aes(x = 0, 
                     y = 1, 
                     xend = 1, 
                     yend = 0), alpha = 0.5) + 
    geom_line(aes(x = plotx_1, 
                  y = ploty_1, 
                  color = "low")) +
    geom_line(aes(x = plotx_2, 
                  y = ploty_2, 
                  color = "high")) +
    scale_color_discrete(name="Contrast",
                         labels=c(paste0("Control vs. Benign\nAUC: ",
                                         round(auc_1, 
                                               digits = 3), 
                                         ", 95% CI: ", 
                                         auc_ci_1[1], 
                                         "-", 
                                         auc_ci_1[3]), 
                                  paste0("Benign vs. Malignant\nAUC: ",
                                         round(auc_2, 
                                               digits = 3), 
                                         ", 95% CI: ", 
                                         auc_ci_2[1], 
                                         "-", 
                                         auc_ci_2[3])),
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
    theme(plot.margin = unit(c(10, 80, 10, 10), "pt"),
          plot.caption.position = "plot") +
    coord_equal() +
    ggtitle(str_c(title)) +
    labs(caption = caption)
} 
