library(tidyverse)
library(ggplot2)
library(ggpubr)
library(yaml)
library(readxl)
library(pROC)
library(yardstick)
library(writexl)
library(caret)

yaml <- yaml.load_file(file.path(getwd(), "config.yaml"))
output_dir <- yaml$demographics$output_dir_figures
# invasive <- read_excel("C:/WorkingData/Documents/3_Research/IVUS_data.xlsx")

# Load the data
baseline <- readRDS(paste0(yaml$demographics$output_dir_data,"/baseline.rds"))

baseline <- baseline %>%
  filter(record_id != 109) %>%
  filter(!is.na(inv_ivusdobu_mla))

baseline <- baseline %>%
  mutate(percent_stenosis = cut(inv_ivusdobu_mla_ln_any, 
                        breaks = c(-Inf, 50, 70, 90, Inf), 
                        labels = c("<50", "50-70", "70-90", ">90")))
# Map the bins to specific sizes
size_values <- c("<50" = 2, "50-70" = 5, "70-90" = 10, ">90" = 20)

# ostium response variable
baseline <- baseline %>% mutate(
  risk_ostium = case_when(
    is.na(inv_ffrdobu) | is.na(inv_ivusdobu_ostial_pn) ~ NA_real_,
    inv_ffrdobu > 0.8 & inv_ivusdobu_ostial_pn < 70 ~ 0,
    inv_ffrdobu > 0.8 & inv_ivusdobu_ostial_pn >= 70 ~ 1,
    inv_ffrdobu <= 0.8 ~ 1
  ),
  risk_imla = case_when(
    is.na(inv_ffrdobu) | is.na(inv_ivusdobu_imla_pn) ~ NA_real_,
    inv_ffrdobu > 0.8 & inv_ivusdobu_imla_pn < 70 ~ 0,
    inv_ffrdobu > 0.8 & inv_ivusdobu_imla_pn >= 70 ~ 1,
    inv_ffrdobu <= 0.8 ~ 1
  ),
  risk_ostium_imla = case_when(
    is.na(inv_ffrdobu) | is.na(inv_ivusdobu_ostial_pn) | is.na(inv_ivusdobu_imla_pn) ~ NA_real_,
    inv_ffrdobu <= 0.8  ~ 1,
    inv_ffrdobu > 0.8 & (inv_ivusdobu_ostial_pn >= 70 | inv_ivusdobu_imla_pn >= 70) ~ 1,
    inv_ffrdobu > 0.8 & (inv_ivusdobu_ostial_pn < 70 & inv_ivusdobu_imla_pn < 70) ~ 0
  ),
  risk_mla_any = case_when(
    is.na(inv_ffrdobu) | is.na(inv_ivusdobu_mla_ln_any) ~ NA_real_,
    inv_ffrdobu > 0.8 & inv_ivusdobu_mla_ln_any < 70 ~ 0,
    inv_ffrdobu > 0.8 & inv_ivusdobu_mla_ln_any >= 70 ~ 1,
    inv_ffrdobu <= 0.8 ~ 1
  ),
  risk_mla_loc = case_when(
    is.na(inv_ffrdobu) | is.na(inv_ivusdobu_mla_ln_loc) ~ NA_real_,
    inv_ffrdobu > 0.8 & inv_ivusdobu_mla_ln_loc < 70 ~ 0,
    inv_ffrdobu > 0.8 & inv_ivusdobu_mla_ln_loc >= 70 ~ 1,
    inv_ffrdobu <= 0.8 ~ 1
  ),
)

# Functions
plot_ffr_data <- function(baseline, var1, var2, method1_name = "Method 1", method2_name = "Method 2", scale_breaks = seq(0.4, 1, 0.1)) {
  n <- nrow(baseline)
  inv_var1 <- pull(baseline, {{var1}})
  inv_var2 <- pull(baseline, {{var2}})
  method_1 <- rep(method1_name, n)
  method_2 <- rep(method2_name, n)
  ffr_df <- data.frame(method = c(method_1, method_2), 
                       value = c(inv_var1, inv_var2))
  colnames(ffr_df) <- c("method", "value")
  
  plot_ffr_change <- ggpaired(ffr_df, x = "method", y = "value", color = "method", 
                              line.color = "#bebebec7", line.size = 0.4, palette = c("#193bac", "#ebc90b")) + 
    stat_compare_means(paired = TRUE) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    scale_y_continuous(breaks = scale_breaks)
  
  scatter_ffr <- ggplot(ffr_df, aes(x = method, y = value)) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
    geom_boxplot(color = c("#193bac", "#ebc90b")) +
    geom_jitter(width = 0.05, color = "#bebebec7") +
    theme_classic() +
    scale_y_continuous(breaks = scale_breaks) + 
    ylab("Value")
  
  return(list(plot_ffr_change, scatter_ffr))
}

plot_with_patho_ffr <- function(baseline, var1, var2, method1_name = "Method 1", method2_name = "Method 2", scale_breaks = seq(0.4, 1, 0.1)) {
  n <- nrow(baseline)
  inv_var1 <- pull(baseline, {{var1}})
  inv_var2 <- pull(baseline, {{var2}})
  inv_ffrado <- pull(baseline, inv_ffrado)
  inv_ffrado <- ifelse(inv_ffrado <= 0.8, 1, 0)
  inv_ffrdobu <- pull(baseline, inv_ffrdobu)
  inv_ffrdobu <- ifelse(inv_ffrdobu <= 0.8, 1, 0)
  method_1 <- rep(method1_name, n)
  method_2 <- rep(method2_name, n)
  ivus_df <- data.frame(method = c(method_1, method_2), 
                       value = c(inv_var1, inv_var2),
                       ffr_value = c(inv_ffrado, inv_ffrdobu))
  ivus_df$ffr_value <- factor(ivus_df$ffr_value, levels = c(0, 1), labels = c("FFR > 0.8", "FFR <= 0.8"))
  ivus_df$method <- factor(ivus_df$method, levels = c(method1_name, method2_name), ordered = TRUE)
  colnames(ivus_df) <- c("method", "value", "ffr_value")
  
  plot_ffr_change <- ggpaired(ivus_df, ivus_df, x = "method", y = "value", color = "method", 
                        line.color = "#bebebec7", line.size = 0.4, palette = c("#ebc90b", "#193bac", "#bebebec7", "red")) +
  geom_point(aes(color = ffr_value)) +
  stat_compare_means(paired = TRUE) +
  scale_y_continuous(breaks = scale_breaks)

  ivus_df$ivus_method <- factor(ivus_df$method, levels = c(method1_name, method2_name), ordered = TRUE)
  scatter_imla <- ggplot(ivus_df, aes(x = method, y = value)) +
    geom_boxplot(aes(color = method)) +
    geom_jitter(width = 0.05, color = "#bebebec7") +
    theme_classic() +
    scale_y_continuous(breaks = scale_breaks) + 
    scale_color_manual(values = c(var1 = "#193bac", var2 = "#ebc90b"))
  
  return(list(plot_ffr_change, scatter_imla))
}

create_boxplot <- function(data, cols, y_breaks) {
  # Reshape the data to long format
  df_long <- pivot_longer(data, cols = all_of(cols), 
                          names_to = "variable", values_to = "value")
  
  # Specify the order of variable levels
  df_long$variable <- factor(df_long$variable, levels = cols)
  
  # Plotting
  plot <- ggplot(df_long, aes(x = record_id, y = value, fill = variable)) +
    geom_boxplot() +
    labs(title = paste("Change from", cols[1], "to", cols[2], "to", cols[3]),
         x = "Record ID",
         y = "Value",
         fill = "Variable") +
    scale_fill_manual(values = c("#193bac", "#ebc90b", "red")) +
    scale_y_continuous(breaks = y_breaks) +
    theme_classic()
  
  return(plot)
}

perform_logistic_regression <- function(response_var, data, ivus_vars) {
  logistic_results <- list()
  
  for (i in 1:length(ivus_vars)) {
    ivus_var <- ivus_vars[i]
    
    p_value <- NA
    odds_ratio <- NA
    ci_lower <- NA
    ci_upper <- NA
    
    tryCatch({
      ivus_glm <- glm(as.formula(paste(response_var, "~", ivus_var)), family = "binomial", data = data)
      
      p_value <- round(summary(ivus_glm)$coefficients[2, 4], 3)
      odds_ratio <- round(exp(coef(ivus_glm)[2]), 2)
      ci_lower <- round(exp(confint(ivus_glm)[2, 1]), 2)
      ci_upper <- round(exp(confint(ivus_glm)[2, 2]), 2)
    }, error = function(e) {
      # If error occurs, leave p_value and odds_ratio as NA
    })
    
    logistic_results[[i]] <- data.frame(
      Variable = ivus_var,
      P_Value = p_value,
      Odds_Ratio = odds_ratio,
      CI_Lower = ci_lower,
      CI_Upper = ci_upper
    )
  }
  
  logistic_results_df <- do.call(rbind, logistic_results)
  
  return(logistic_results_df)
}

perform_roc_analysis <- function(response_var, data, ivus_vars) {
  roc_results <- list()
  
  for (i in 1:length(ivus_vars)) {
    ivus_var <- ivus_vars[i]
    
    auc <- NA
    sensitivity <- NA
    specificity <- NA
    threshold <- NA
    
    tryCatch({
      roc_model <- glm(as.formula(paste(response_var, "~", ivus_var)), data = data, family = binomial)
      
      roc_pred <- predict(roc_model, data, type = "response")
      
      roc_curve <- roc(data[[response_var]], roc_pred)
      
      auc <- round(auc(roc_curve), 2)
      
      roc_threshold <- roc(data[[response_var]], data[[ivus_var]])
      opt_threshold <- coords(roc_threshold, "best", ret = "threshold")
      threshold <- round(opt_threshold, 2)
      
      sens_spec <- coords(roc_curve, "best", ret = c("sensitivity", "specificity"))
      sensitivity <- round(sens_spec["sensitivity"], 2)
      specificity <- round(sens_spec["specificity"], 2)
    }, error = function(e) {
      # If error occurs, leave variables as NA
    })    

    roc_results[[i]] <- data.frame(
      Variable = ivus_var,
      AUC = auc,
      Sensitivity = sensitivity,
      Specificity = specificity,
      Threshold = threshold
    )
  }

  roc_results_df <- do.call(rbind, roc_results)
  
  return(roc_results_df)
}

simple_logistic_regression <- function(baseline_data, response, explanatory, family = "binomial") {
    data <- baseline_data %>% 
        select(!!sym(response), !!sym(explanatory)) %>% 
        drop_na()
    formula <- as.formula(paste(response, "~", explanatory))
    mdl <- glm(formula, data = data, family = family)

    max <- max(data[[explanatory]], na.rm = TRUE)
    min <- min(data[[explanatory]], na.rm = TRUE)

    if (min<0) {
        min <- 0
    }

    explanatory_data <- tibble(
        !!sym(explanatory) := seq(min, max, length.out = 1000)
    )
    
    prediction_data <- explanatory_data %>% 
        mutate(
            response = predict(mdl, explanatory_data, type = "response"),
            most_likely_outcome = round(response),
            odds_ratio = response / (1 - response),
            log_odds_ratio = log(odds_ratio),
            log_odds_ratio2 = predict(mdl, explanatory_data)
        )

    # because of bug rename response here
    names(prediction_data)[2] <- response
    
    plot_log <- ggplot(data, aes(x = !!sym(explanatory), y = !!sym(response))) +
        geom_point() + 
        geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) +
        scale_x_continuous(breaks = seq(round(min, 0), round(max), round(max / 10, 1)))
    
    plot_odds <- ggplot(prediction_data, aes(x = !!sym(explanatory), y = odds_ratio)) + 
        geom_point() + 
        geom_line() +
        geom_hline(yintercept = 1, linetype = "dashed") +
        scale_x_continuous(breaks = seq(round(min, 0), round(max), round(max / 10, 1))) +
        scale_y_continuous(breaks = seq(0, round(max(prediction_data$odds_ratio)), 1))

    confusion <- NULL
    plot_acc <- NULL
    tryCatch(
        {
            actual_response <- data[[response]]
            predicted_response <- round(fitted(mdl))
            outcome <- table(predicted_response, actual_response)
            confusion <- conf_mat(outcome)
            plot_acc <- autoplot(confusion)
        },
        error = function(e) {
            print("Check model! No binary response")
        }
    )
    return(list(mdl, prediction_data, plot_log, plot_odds, confusion, plot_acc))
}

# FFR data
n <- nrow(baseline)
inv_var1 <- pull(baseline, inv_ffrado)
inv_var2 <- pull(baseline, inv_ffrdobu)
inv_ffrado <- pull(baseline, inv_ffrado)
inv_ffrado <- ifelse(inv_ffrado <= 0.8, 1, 0)
inv_ffrdobu <- pull(baseline, inv_ffrdobu)
inv_ffrdobu <- ifelse(inv_ffrdobu <= 0.8, 1, 0)
method_1 <- rep("FFRadenosine", n)
method_2 <- rep("FFRdobutamine", n)
ivus_df <- data.frame(method = c(method_1, method_2), 
                      value = c(inv_var1, inv_var2),
                      ffr_value = c(inv_ffrado, inv_ffrdobu))
ivus_df$ffr_value <- factor(ivus_df$ffr_value, levels = c(0, 1), labels = c("FFR > 0.8", "FFR <= 0.8"))
ivus_df$method <- factor(ivus_df$method, levels = c("FFRadenosine", "FFRdobutamine"), ordered = TRUE)
colnames(ivus_df) <- c("method", "value", "ffr_value")

plot_ffr_line <- ggpaired(ivus_df, ivus_df, x = "method", y = "value", color = "method", 
                        line.color = "#bebebec7", line.size = 0.4, palette = c("#193bac", "#ebc90b", "#bebebec7", "red")) +
  geom_point(aes(color = ffr_value)) +
  stat_compare_means(paired = TRUE) +
  scale_y_continuous(breaks = seq(0, 1, 0.05))  + 
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") + # Plot 1
  ylab("FFR")

pressure_rest_ado_dobu <- create_boxplot(baseline, c("inv_rfr", "inv_ffrado", "inv_ffrdobu"), seq(0, 1, 0.1))
create_boxplot(baseline, c("inv_ivusrest_mla", "inv_ivusado_mla", "inv_ivusdobu_mla"), seq(0, 15, 1))
ln_rest_ado_dobu <- create_boxplot(baseline, c("inv_ivusrest_mla_ln", "inv_ivusado_mla_ln", "inv_ivusdobu_mla_ln_any"), seq(0, 100, 10))

ggsave(paste0(output_dir,"/pressure_rest_ado_dobu.png"), pressure_rest_ado_dobu, width = 6, height = 5)
ggsave(paste0(output_dir,"/ln_rest_ado_dobu.png"), ln_rest_ado_dobu, width = 6, height = 5)

# MLA
plot_mla <- plot_with_patho_ffr(baseline, var1 = inv_ivusrest_mla, var2 = inv_ivusdobu_mla,
                       method1_name = "IVUSrest", method2_name = "IVUSdobutamine", scale_breaks = seq(0, 15, 1))[[1]] +
                       ylab("Minimal lumen area [mm²]")

plot_mla_ln_any <- plot_with_patho_ffr(baseline, var1 = inv_ivusrest_mla_ln, var2 = inv_ivusdobu_mla_ln_any,
                       method1_name = "IVUSrest", method2_name = "IVUSdobutamine", scale_breaks = seq(0, 100, 10))[[1]] +
  geom_hline(yintercept = 70, linetype = "dashed", color = "red") +
  ylab("Minimal lumen narrowing [%]")

plot_mla_ln_loc <- plot_with_patho_ffr(baseline, var1 = inv_ivusrest_mla_ln, var2 = inv_ivusdobu_mla_ln_loc,
                       method1_name = "IVUSrest", method2_name = "IVUSdobutamine", scale_breaks = seq(0, 100, 10))[[1]] +
  geom_hline(yintercept = 70, linetype = "dashed", color = "red") +
  ylab("Minimal lumen narrowing [%]")

# Ostial data
plot_ostial_a <- plot_with_patho_ffr(baseline, var1 = inv_ivusrest_ostial_a, var2 = inv_ivusdobu_ostial_a,
                       method1_name = "IVUSrest", method2_name = "IVUSdobutamine", scale_breaks = seq(0, 15, 1))[[1]] +
                       ylab("Ostial lumen area [mm²]")

plot_ostial_pn <- plot_with_patho_ffr(baseline, var1 = inv_ivusrest_ostial_pn, var2 = inv_ivusdobu_ostial_pn,
                       method1_name = "IVUSrest", method2_name = "IVUSdobutamine", scale_breaks = seq(0, 100, 10))[[1]] +
  geom_hline(yintercept = 70, linetype = "dashed", color = "red") +
  ylab("Ostial lumen narrowing [%]")

# MLA data
plot_imla <- plot_with_patho_ffr(baseline, var1 = inv_ivusrest_imla, var2 = inv_ivusdobu_imla,
                       method1_name = "IVUSrest", method2_name = "IVUSdobutamine", scale_breaks = seq(0, 15, 1))[[1]] +
                       ylab("Intramural minimal lumen area [mm²]")

plot_imla_pn <- plot_with_patho_ffr(baseline, var1 = inv_ivusrest_imla_pn, var2 = inv_ivusdobu_imla_pn,
                       method1_name = "IVUSrest", method2_name = "IVUSdobutamine", scale_breaks = seq(0, 100, 10))[[1]] +
  geom_hline(yintercept = 70, linetype = "dashed", color = "red") +
  ylab("Intramural lumen narrowing [%]")

ggsave(paste0(output_dir,"/plot_ffr_line.png"), plot_ffr_line, width = 4, height = 5)
ggsave(paste0(output_dir,"/plot_mla.png"), plot_mla, width = 4, height = 5)
ggsave(paste0(output_dir,"/plot_mla_ln_any.png"), plot_mla_ln_any, width = 4, height = 5)
ggsave(paste0(output_dir,"/plot_mla_ln_loc.png"), plot_mla_ln_loc, width = 4, height = 5)
ggsave(paste0(output_dir,"/plot_ostial_a.png"), plot_ostial_a, width = 4, height = 5)
ggsave(paste0(output_dir,"/plot_ostial_pn.png"), plot_ostial_pn, width = 4, height = 5)
ggsave(paste0(output_dir,"/plot_imla.png"), plot_imla, width = 4, height = 5)
ggsave(paste0(output_dir,"/plot_imla_pn.png"), plot_imla_pn, width = 4, height = 5)

# linear relationships
mla_ffr <- ggplot(baseline, aes(x = inv_ivusdobu_mla_bsa, y = inv_ffrdobu)) +
  geom_smooth(se = FALSE, color = "grey", linetype = "dashed") +  # Single smooth line for all points
  geom_point(aes(size = percent_stenosis, alpha = 0.7)) +  # Map color to the points
  theme_classic() + 
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  scale_y_continuous(breaks = seq(0, 1, 0.1)) +
  scale_color_identity() + 
  scale_size_manual(values = size_values) +
  labs(size = "Percent stenosis [%]") +
  xlab("Minimal lumen area (BSA adjusted) [mm²]") +
  ylab("FFR")

ggsave(paste0(output_dir,"/mla_ffr.png"), mla_ffr, width = 6, height = 5)

# risk groups
# Create a new column 'color' in the baseline dataframe based on the given conditions
baseline$color <- with(baseline, ifelse(inv_ffrdobu <= 0.8 & inv_ivusdobu_mla_ln_any >= 70, 'red',
                               ifelse(inv_ffrdobu <= 0.8 | inv_ivusdobu_mla_ln_any >= 70, '#ff7504', '#094609')))

baseline <- baseline %>% 
  mutate(synp_revasc = case_when(
    synp_surgery_yn == "yes" | synp_stent_yn == "yes" ~ "yes",
    is.na(synp_surgery_yn) & is.na(synp_stent_yn) ~ "no",
    TRUE ~ "no"
  ))

# replace value of patient_id 264 in synp_revasc with yes
baseline$synp_revasc[baseline$patient_id == "NARCO_264"] <- "yes"
baseline %>% filter(inv_ffrdobu <0.8) %>% select(patient_id, inv_ffrdobu, synp_revasc)

# calculate if the stenosis was present in rest already or developed afterwards for >70% stenosis
baseline <- baseline %>% mutate(
    dynamic_stenosis = case_when(
        inv_ivusrest_ostial_pn < 70 & inv_ivusdobu_ostial_pn >= 70 ~ "yes",
        inv_ivusrest_ostial_pn >= 70 & inv_ivusdobu_ostial_pn >= 70 ~ "no",
        TRUE ~ "no"
    )
)

# Plot using the new 'color' column
mla_ln_ffr <- ggplot(baseline, aes(x = inv_ivusdobu_mla_ln_any, y = inv_ffrdobu)) +
  # geom_point(aes(size = inv_ivusdobu_mla_a_bsa, color = color)) +  # Only map color to the points
  geom_smooth(se = FALSE, color = "grey", linetype = "dashed", alpha = 0.5) +  # Single smooth line for all points
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 70, linetype = "dashed", color = "black") +
  geom_point(aes(color = color)) +
#   geom_point(aes(color = color, shape = factor(synp_revasc))) + 
  theme_classic() + 
  scale_y_continuous(breaks = seq(0, 1, 0.05)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_color_identity() + 
  scale_shape_manual(values = c("yes" = 13, "no" = 16)) +
  xlab("Minimal lumen narrowing [%]") +
  ylab("FFR")

mla_ln_ffr_dynamic <- ggplot(baseline, aes(x = inv_ivusdobu_mla_ln_any, y = inv_ffrdobu)) +
  geom_smooth(se = FALSE, color = "grey", linetype = "dashed", alpha = 0.5) +  # Single smooth line for all points
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 70, linetype = "dashed", color = "black") +
  geom_point(aes(color = color, shape = factor(dynamic_stenosis))) +
  theme_classic() + 
  scale_y_continuous(breaks = seq(0, 1, 0.05)) +
  scale_x_continuous(breaks = seq(0, 100, 10)) +
  scale_color_identity() + 
  scale_shape_manual(values = c("yes" = 13, "no" = 16)) +
  labs(shape = "Dynamic stenosis") +
  xlab("Minimal lumen narrowing [%]") +
  ylab("FFR")

ggsave(paste0(output_dir,"/mla_ln_ffr.png"), mla_ln_ffr, width = 5, height = 4)
ggsave(paste0(output_dir,"/mla_ln_ffr_dynamic.png"), mla_ln_ffr_dynamic, width = 5, height = 3.5)

# Define the variables for comparison
rest_vars <- c("inv_ffrado", 
                "inv_rest_hr", 
                "inv_rest_aosys", 
                "inv_rest_aodia",
               "inv_rest_aomean", 
               "inv_ivusrest_ostial_a",
                "inv_ivusrest_imla", 
                "inv_ivusrest_pn",
               "inv_ivusrest_imla_bsa", 
               "inv_ivusrest_ostial_a_bsa",
               "inv_ivusrest_ostial_pn", 
               "inv_ivusrest_imla_pn", 
               "inv_ivusrest_mla",
               "inv_ivusrest_mla_bsa", 
               "inv_ivusrest_mla_ln")

dobu_vars <- c("inv_ffrdobu", 
                "inv_dobu_hr",
               "inv_dobu_aosys", 
               "inv_dobu_aodia", 
               "inv_dobu_aomean", 
               "inv_ivusdobu_ostial_a", 
               "inv_ivusdobu_imla", 
               "inv_ivusdobu_pn", 
               "inv_ivusdobu_imla_bsa", 
               "inv_ivusdobu_ostial_a_bsa", 
               "inv_ivusdobu_ostial_pn", 
               "inv_ivusdobu_imla_pn", 
               "inv_ivusdobu_mla",
               "inv_ivusdobu_mla_bsa", 
               "inv_ivusdobu_mla_ln_any")

# Initialize a list to store results
results <- list()

# Loop over variables to perform the tests and compute statistics
for (i in 1:length(rest_vars)) {
  rest_var <- rest_vars[i]
  dobu_var <- dobu_vars[i]
  
  # Initialize variables to store test results
  rest_summary <- ""
  dobu_summary <- ""
  test_method <- ""
  p_value <- NA
  
  # Try performing Shapiro-Wilk test
  tryCatch({
    rest_shapiro <- shapiro.test(as.numeric(baseline[[rest_var]]))
    dobu_shapiro <- shapiro.test(as.numeric(baseline[[dobu_var]]))
    
    # Check for normality
    if (rest_shapiro$p.value > 0.05 && dobu_shapiro$p.value > 0.05) {
      # Both variables are normally distributed, use paired t-test
      test_result <- t.test(baseline[[rest_var]], baseline[[dobu_var]], paired = TRUE)
      # Compute mean ± SD
      rest_mean <- mean(baseline[[rest_var]], na.rm = TRUE)
      rest_sd <- sd(baseline[[rest_var]], na.rm = TRUE)
      dobu_mean <- mean(baseline[[dobu_var]], na.rm = TRUE)
      dobu_sd <- sd(baseline[[dobu_var]], na.rm = TRUE)
      rest_summary <- sprintf("%.2f ± %.2f", rest_mean, rest_sd)
      dobu_summary <- sprintf("%.2f ± %.2f", dobu_mean, dobu_sd)
      test_method <- "t-test"
    } else {
      # Use Wilcoxon signed-rank test
      test_result <- wilcox.test(baseline[[rest_var]], baseline[[dobu_var]], paired = TRUE)
      # Compute median (Q1 - Q3)
      rest_median <- median(baseline[[rest_var]], na.rm = TRUE)
      rest_q1 <- quantile(baseline[[rest_var]], 0.25, na.rm = TRUE)
      rest_q3 <- quantile(baseline[[rest_var]], 0.75, na.rm = TRUE)
      dobu_median <- median(baseline[[dobu_var]], na.rm = TRUE)
      dobu_q1 <- quantile(baseline[[dobu_var]], 0.25, na.rm = TRUE)
      dobu_q3 <- quantile(baseline[[dobu_var]], 0.75, na.rm = TRUE)
      rest_summary <- sprintf("%.2f (%.2f-%.2f)", rest_median, rest_q1, rest_q3)
      dobu_summary <- sprintf("%.2f (%.2f-%.2f)", dobu_median, dobu_q1, dobu_q3)
      test_method <- "Wilcoxon"
    }
    p_value <- test_result$p.value
  }, error = function(e) {
    # Fill in blank line if error occurs
    rest_var <- ""
    dobu_var <- ""
  })
  
  # Store results
  results[[i]] <- data.frame(
    Variable = rest_var,
    Rest_Summary = rest_summary,
    Dobu_Summary = dobu_summary,
    Test = test_method,
    P_Value = round(p_value,3)
  )
}

# Combine results into a single data frame
results_df <- do.call(rbind, results)

# Print the results
print(results_df)

# Optionally, save the results to a CSV file
write.csv(results_df, "statistical_analysis/data/invasive_change.csv", row.names = FALSE)

# baseline IVUS variables to predict FFR change
ivus_vars <- c("inv_ivusrest_ostial_a",
               "inv_ivusrest_ostial_a_bsa",
               "inv_ivusrest_ostial_pn", 
               "inv_ivusrest_imla",
               "inv_ivusrest_imla_bsa", 
               "inv_ivusrest_imla_pn",
               "inv_ivusrest_mla", 
               "inv_ivusrest_mla_bsa", 
               "inv_ivusrest_mla_ln", "inv_ffrado")

logistic_results_ffr_0.8 <- perform_logistic_regression("ffr_0.8", baseline, ivus_vars)
logistic_results_risk_ostium <- perform_logistic_regression("risk_ostium", baseline, ivus_vars)
logistic_results_risk_imla <- perform_logistic_regression("risk_imla", baseline, ivus_vars)
logistic_results_risk_ostium_imla <- perform_logistic_regression("risk_ostium_imla", baseline, ivus_vars)
logistic_results_risk_mla_any <- perform_logistic_regression("risk_mla_any", baseline, ivus_vars)
logistic_results_risk_mla_loc <- perform_logistic_regression("risk_mla_loc", baseline, ivus_vars)

results_list <- list(
  "ffr_0.8" = logistic_results_ffr_0.8,
  "risk_ostium" = logistic_results_risk_ostium,
  "risk_imla" = logistic_results_risk_imla,
  "risk_ostium_imla" = logistic_results_risk_ostium_imla,
  "risk_mla_any" = logistic_results_risk_mla_any,
  "risk_mla_loc" = logistic_results_risk_mla_loc
)

# Write the list of data frames to an Excel file
write_xlsx(results_list, "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/ivus_logistic_regression_results.xlsx")

######### ROC analysis for IVUS variables ###############################################################
# Perform ROC analysis for all the specified response variables
roc_results_ffr_0.8 <- perform_roc_analysis("ffr_0.8", baseline, ivus_vars)
roc_results_risk_ostium <- perform_roc_analysis("risk_ostium", baseline, ivus_vars)
roc_results_risk_imla <- perform_roc_analysis("risk_imla", baseline, ivus_vars)
roc_results_risk_ostium_imla <- perform_roc_analysis("risk_ostium_imla", baseline, ivus_vars)
roc_results_risk_imla_50 <- perform_roc_analysis("risk_mla_any", baseline, ivus_vars)
roc_results_risk_ostium_50 <- perform_roc_analysis("risk_mla_loc", baseline, ivus_vars)

# Create a list of data frames to be written to the Excel file
roc_results_list <- list(
  "ffr_0.8" = roc_results_ffr_0.8,
  "risk_ostium" = roc_results_risk_ostium,
  "risk_imla" = roc_results_risk_imla,
  "risk_ostium_imla" = roc_results_risk_ostium_imla,
  "risk_mla_any" = roc_results_risk_imla_50,
  "risk_mla_loc" = roc_results_risk_ostium_50
)

# Write the list of data frames to an Excel file
write_xlsx(roc_results_list, "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/ivus_roc_analysis_results.xlsx")

############################################################################################################
baseline <- baseline %>% mutate(
  risk_group = case_when(
    inv_ivusdobu_mla_ln_any < 70 & inv_ffrdobu > 0.8  ~ 0,
    inv_ivusdobu_mla_ln_any >= 70 & inv_ffrdobu > 0.8 ~ 1,
    inv_ivusdobu_mla_ln_any < 70 & inv_ffrdobu <= 0.8 ~ 2,
    inv_ivusdobu_mla_ln_any >= 70 & inv_ffrdobu <= 0.8 ~ 3
  ) 
)

baseline$risk_group <- factor(baseline$risk_group, 
                              levels = c(0, 1, 2, 3), 
                              labels = c("FFR>0.8 & stenosis<70%", 
                                         "FFR>0.8 & stenosis>=70%", 
                                         "FFR<=0.8 & stenosis<70%", 
                                         "FFR<=0.8 & stenosis>=70%"), 
                              ordered = TRUE)

# bar graph with risk groups with n per group
risk_groups_plot <- baseline %>%
  filter(!is.na(risk_group)) %>%
  ggplot(aes(x = risk_group, fill = risk_group)) +
  geom_bar() +
  geom_text(stat = "count", aes(label = ..count..), vjust = -0.5) +
  theme_classic() +
  scale_fill_manual(values = c("FFR>0.8 & stenosis<70%" = "#094609", 
                               "FFR>0.8 & stenosis>=70%" = "#ff7504", 
                               "FFR<=0.8 & stenosis<70%" = "#ff7504", 
                               "FFR<=0.8 & stenosis>=70%" = "red")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(paste0(output_dir,"/risk_groups_plot.png"), risk_groups_plot, width = 7, height = 5)

############################################################################################################
# best cutoffs
# Fit models and get prediction data
mdl_ost_a_bsa_ffr <- simple_logistic_regression(baseline, "risk_mla_any", "inv_ivusrest_ostial_a_bsa")[[1]]
mdl_ost_pn_ffr <- simple_logistic_regression(baseline, "risk_mla_any", "inv_ivusrest_ostial_pn")[[1]]
mdl_imla_bsa_ffr <- simple_logistic_regression(baseline, "risk_mla_any", "inv_ivusrest_imla_bsa")[[1]]
mdl_imla_pn_ffr <- simple_logistic_regression(baseline, "risk_mla_any", "inv_ivusrest_imla_pn")[[1]]
mdl_mla_bsa_ffr <- simple_logistic_regression(baseline, "risk_mla_any", "inv_ivusrest_mla_bsa")[[1]]
mdl_mla_ln_ffr <- simple_logistic_regression(baseline, "risk_mla_any", "inv_ivusrest_mla_ln")[[1]]
mdl_ffrado <- simple_logistic_regression(baseline, "risk_mla_any", "inv_ffrado")[[1]]

prediction_data_ost_a <- simple_logistic_regression(baseline, "risk_mla_any", "inv_ivusrest_ostial_a_bsa")[[2]]
prediction_data_ost_pn <- simple_logistic_regression(baseline, "risk_mla_any", "inv_ivusrest_ostial_pn")[[2]]
prediction_data_imla <- simple_logistic_regression(baseline, "risk_mla_any", "inv_ivusrest_imla_bsa")[[2]]
prediction_data_imla_pn <- simple_logistic_regression(baseline, "risk_mla_any", "inv_ivusrest_imla_pn")[[2]]
prediction_data_mla <- simple_logistic_regression(baseline, "risk_mla_any", "inv_ivusrest_mla_bsa")[[2]]
prediction_data_mla_ln <- simple_logistic_regression(baseline, "risk_mla_any", "inv_ivusrest_mla_ln")[[2]]
prediction_data_ffrado <- simple_logistic_regression(baseline, "risk_mla_any", "inv_ffrado")[[2]]


data_ivus <- baseline %>% select(risk_mla_any, inv_ivusrest_ostial_a) %>% drop_na()
roc_ost_a <- roc(data_ivus$risk_mla_any, mdl_ost_a_bsa_ffr$fitted.values)
data_ivus <- baseline %>% select(risk_mla_any, inv_ivusrest_ostial_pn) %>% drop_na()
roc_ost_pn <- roc(data_ivus$risk_mla_any, mdl_ost_pn_ffr$fitted.values)
data_ivus <- baseline %>% select(risk_mla_any, inv_ivusrest_imla_bsa) %>% drop_na()
roc_imla <- roc(data_ivus$risk_mla_any, mdl_imla_bsa_ffr$fitted.values)
data_ivus <- baseline %>% select(risk_mla_any, inv_ivusrest_imla_pn) %>% drop_na()
roc_imla_pn <- roc(data_ivus$risk_mla_any, mdl_imla_pn_ffr$fitted.values)
data_ivus <- baseline %>% select(risk_mla_any, inv_ivusrest_mla_bsa) %>% drop_na()
roc_mla <- roc(data_ivus$risk_mla_any, mdl_mla_bsa_ffr$fitted.values)
data_ivus <- baseline %>% select(risk_mla_any, inv_ivusrest_mla_ln) %>% drop_na()
roc_mla_ln <- roc(data_ivus$risk_mla_any, mdl_mla_ln_ffr$fitted.values)
data_ivus <- baseline %>% select(risk_mla_any, inv_ffrado) %>% drop_na()
roc_ffrado <- roc(data_ivus$risk_mla_any, mdl_ffrado$fitted.values)

# Define test sets
test_set <- list(
  coords(roc_ost_a), coords(roc_ost_pn), 
  coords(roc_imla), coords(roc_imla_pn),
  coords(roc_mla), coords(roc_mla_ln),
  coords(roc_ffrado)
)

prediction_data <- list(
  prediction_data_ost_a, prediction_data_ost_pn, 
  prediction_data_imla, prediction_data_imla_pn,
  prediction_data_mla, prediction_data_mla_ln,
  prediction_data_ffrado
)

variables <- c(
  "IVUS rest ostial area (BSA)", "IVUS rest ostial luminal narrowing", 
  "IVUS rest IMLA (BSA)", "IVUS rest IMLA luminal narrowing",
  "IVUS rest MLA (BSA)", "IVUS rest MLA luminal narrowing",
  "FFR adenosine"
)

# Extract thresholds
thresholds <- list()

for (i in 1:length(test_set)) {
  coordinates <- test_set[[i]]
  coordinates <- coordinates[-c(1, nrow(coordinates)), ]
  coordinates <- as_tibble(coordinates)
  
  result <- coordinates %>%
    group_by(sensitivity) %>%
    filter(specificity == max(specificity)) %>%
    select(threshold) %>%
    distinct()
  
  threshold_values <- result$threshold
  thresholds[[variables[i]]] <- threshold_values
}

# Generate confusion matrices and calculate metrics
cutoff_tables <- list()

thresholds_testing <- thresholds[[1]]
prediction <- prediction_data[[1]]
name <- names(prediction)[1]
threshold <- thresholds_testing[1]
for (i in 1:length(thresholds_testing)) {
  threshold <- thresholds_testing[i]
  closest_index <- which.min(abs(prediction$risk_mla_any - threshold))
  value <- prediction[[name]][closest_index]
}

for (i in 1:length(thresholds)) {
  thresholds_testing <- thresholds[[i]]
  prediction <- prediction_data[[i]]
  name <- names(prediction)[1]
  
  table_thresholds <- tibble()
  for (j in 1:length(thresholds_testing)) {
    threshold <- thresholds_testing[j]
    closest_index <- which.min(abs(prediction$risk_mla_any - threshold))
    value <- prediction[[name]][closest_index]
    
    if (name %in% c("inv_ivusrest_ostial_a_bsa", "inv_ivusrest_imla_bsa", "inv_ivusrest_mla_bsa", "inv_ffrado")) {
      table <- baseline %>%
        mutate(pred_label = ifelse(!!sym(name) > value, 0, 1)) %>%
        select(pred_label, risk_mla_any) %>%
        drop_na() %>%
        table()
    } else {
      table <- baseline %>%
        mutate(pred_label = ifelse(!!sym(name) < value, 0, 1)) %>%
        select(pred_label, risk_mla_any) %>%
        drop_na() %>%
        table()
    }
    
    confusion <- conf_mat(table)
    result <- summary(confusion, event_level = "second")
    
    sens <- result %>% filter(.metric == "sens") %>% select(.estimate) %>% pull()
    spec <- result %>% filter(.metric == "spec") %>% select(.estimate) %>% pull()
    ppv <- result %>% filter(.metric == "ppv") %>% select(.estimate) %>% pull()
    npv <- result %>% filter(.metric == "npv") %>% select(.estimate) %>% pull()
    acc <- result %>% filter(.metric == "accuracy") %>% select(.estimate) %>% pull()
    
    df <- tibble(
      value = as.numeric(round(value, 2)),
      sensitivity = as.numeric(round(sens * 100, 0)),
      specificity = as.numeric(round(spec * 100, 0)),
      ppv = as.numeric(round(ppv * 100, 0)),
      npv = as.numeric(round(npv * 100, 0)),
      accuracy = as.numeric(round(acc * 100, 0)),
      number_of_true_negatives = as.numeric(table[1, 1]),
      number_of_false_negatives = as.numeric(table[1, 2]),
      number_of_false_positives = as.numeric(table[2, 1]),
      number_of_true_positives = as.numeric(table[2, 2])
    )
    names(df)[1] <- name
    table_thresholds <- bind_rows(table_thresholds, df)
  }
  cutoff_tables[[variables[i]]] <- table_thresholds
}

# Save all as CSV
path <- "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/"
for (i in 1:length(cutoff_tables)) {
  name <- variables[i]
  table <- cutoff_tables[[i]]
  write_csv(table, paste0(path, name, ".csv"))
}

roc_list_ost <- list(roc_ost_a, roc_ost_pn)
ggroc_ost <- ggroc(roc_list_ost, legacy.axes = TRUE) +
    scale_color_manual(values = c("darkred", "darkblue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Ostial lumen") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC OLA =", round(roc_ost_a$auc, 2)), color = "darkred") + 
    annotate("text", x = 0.7, y = 0.2, label = paste("AUC OLN =", round(roc_ost_pn$auc, 2)), color = "darkblue")

roc_list_imla <- list(roc_imla, roc_imla_pn)
ggroc_imla <- ggroc(roc_list_imla, legacy.axes = TRUE) +
    scale_color_manual(values = c("darkred", "darkblue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Intramural lumen") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC IMLA =", round(roc_imla$auc, 2)), color = "darkred") + 
    annotate("text", x = 0.7, y = 0.2, label = paste("AUC IMLN =", round(roc_imla_pn$auc, 2)), color = "darkblue")

roc_list_mla <- list(roc_mla, roc_mla_ln)
ggroc_mla <- ggroc(roc_list_mla, legacy.axes = TRUE) +
    scale_color_manual(values = c("darkred", "darkblue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Minimal lumen") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC MLA =", round(roc_mla$auc, 2)), color = "darkred") + 
    annotate("text", x = 0.7, y = 0.2, label = paste("AUC MLN =", round(roc_mla_ln$auc, 2)), color = "darkblue")

roc_list <- list(roc_ost_a, roc_imla, roc_mla)
ggroc_all_a <- ggroc(roc_list, legacy.axes = TRUE) +
    scale_color_manual(values = c("darkred", "darkblue", "darkgreen")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Lumen area comparison") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC OLA =", round(roc_ost_a$auc, 2)), color = "darkred") + 
    annotate("text", x = 0.7, y = 0.2, label = paste("AUC IMLA =", round(roc_imla$auc, 2)), color = "darkblue") +
    annotate("text", x = 0.7, y = 0.1, label = paste("AUC MLA =", round(roc_mla$auc, 2)), color = "darkgreen") +
    annotate("text", x = 0.1, y = 1, label = paste("p-value (OLA versus IMLA) =", round(roc.test(roc_ost_a, roc_imla, method = "delong")$p.value, 3))) +
    annotate("text", x = 0.1, y = 0.9, label = paste("p-value (IMLA versus MLA) =", round(roc.test(roc_imla, roc_mla, method = "delong")$p.value, 3))) +
    annotate("text", x = 0.1, y = 0.8, label = paste("p-value (OLA versus MLA) =", round(roc.test(roc_ost_a, roc_mla, method = "delong")$p.value, 3)))

roc_list <- list(roc_ost_pn, roc_imla_pn, roc_mla_ln)
ggroc_all_ln <- ggroc(roc_list, legacy.axes = TRUE) +
    scale_color_manual(values = c("darkred", "darkblue", "darkgreen")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Lumen narrowing comparison") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC OLN =", round(roc_ost_pn$auc, 2)), color = "darkred") + 
    annotate("text", x = 0.7, y = 0.2, label = paste("AUC IMLN =", round(roc_imla_pn$auc, 2)), color = "darkblue") +
    annotate("text", x = 0.7, y = 0.1, label = paste("AUC MLN =", round(roc_mla_ln$auc, 2)), color = "darkgreen") +
    annotate("text", x = 0.1, y = 1, label = paste("p-value (OLA versus IMLA) =", round(roc.test(roc_ost_pn, roc_imla_pn, method = "delong")$p.value, 3))) +
    annotate("text", x = 0.1, y = 0.9, label = paste("p-value (IMLA versus MLA) =", round(roc.test(roc_imla_pn, roc_mla_ln, method = "delong")$p.value, 3))) +
    annotate("text", x = 0.1, y = 0.8, label = paste("p-value (OLA versus MLA) =", round(roc.test(roc_ost_pn, roc_mla_ln, method = "delong")$p.value, 3)))

ggsave(paste0(output_dir,"/roc_ost.png"), ggroc_ost, width = 6, height = 5)
ggsave(paste0(output_dir,"/roc_imla.png"), ggroc_imla, width = 6, height = 5)
ggsave(paste0(output_dir,"/roc_mla.png"), ggroc_mla, width = 6, height = 5)
ggsave(paste0(output_dir,"/roc_all_area.png"), ggroc_all_a, width = 6, height = 5)
ggsave(paste0(output_dir,"/roc_all_ln.png"), ggroc_all_ln, width = 6, height = 5)

# Ensure dynamic_stenosis is numeric
baseline <- baseline %>%
  mutate(dynamic_stenosis = ifelse(dynamic_stenosis == "yes", 1, 0))

# Calculate the total counts and dynamic_stenosis proportions per risk_group
total_count <- baseline %>% filter(!is.na(inv_ivusrest_ostial_a)) %>% nrow()

plot_data <- baseline %>%
  filter(!is.na(risk_group)) %>%
  group_by(risk_group) %>%
  summarise(count = n(),
            dynamic_stenosis_sum = sum(dynamic_stenosis, na.rm = TRUE),
            dynamic_stenosis_percent = mean(dynamic_stenosis, na.rm = TRUE)) 

# Create the plot with the percentage bars
risk_groups_plot_dynamic <- plot_data %>%
  ggplot(aes(x = risk_group, fill = risk_group)) +
  geom_bar(aes(y = count), stat = "identity") +
  geom_bar(aes(y = count * dynamic_stenosis_percent), fill = "white", alpha = 0.5, stat = "identity") +
  geom_text(aes(y = count, label = count), vjust = -0.5) +
  geom_text(aes(y = count * dynamic_stenosis_percent, label = paste0(round(dynamic_stenosis_percent * 100, 0), "%")), vjust = -0.5, color = "grey") +
  theme_classic() +
  scale_fill_manual(values = c("FFR>0.8 & stenosis<70%" = "#094609", 
                               "FFR>0.8 & stenosis>=70%" = "#ff7504", 
                               "FFR<=0.8 & stenosis<70%" = "#ff7504", 
                               "FFR<=0.8 & stenosis>=70%" = "red",
                               "Dummy" = "grey"),
                    labels = c("FFR>0.8 & stenosis<70%", 
                               "FFR>0.8 & stenosis>=70%", 
                               "FFR<=0.8 & stenosis<70%", 
                               "FFR<=0.8 & stenosis>=70%",
                               "Dynamic Stenosis Percentage")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(fill = "Risk Group") +
  scale_y_continuous(
    name = "Count",
    sec.axis = sec_axis(~ . / total_count * 100, name = "Percentage"),
    labels = scales::number_format(accuracy = 1)
  )

ggsave(paste0(output_dir,"/risk_groups_plot_dynamic.png"), risk_groups_plot_dynamic, width = 7, height = 5)

sum(baseline$inv_ivusrest_mla_ln >= 70, na.rm = T)

# proportion test between 9/28 and 3/26
prop.test(x = c(9, 3), n = c(28, 26), alternative = "two.sided", correct = FALSE)

baseline <- baseline %>% mutate(
  percent_compression = (inv_ivusrest_mla - inv_ivusdobu_mla) / inv_ivusrest_mla * 100
)

# get number of inv_ivusdobu_mla_ln_any >=70% that were inv_ivusrest_mla <70%
sum(baseline$inv_ivusdobu_mla_ln_any >= 70 & baseline$inv_ivusrest_mla_ln < 70, na.rm = T)

shapiro.test(baseline$inv_rfr)
median(baseline$inv_rfr, na.rm = T)
quantile(baseline$inv_rfr, c(0.25, 0.75), na.rm = T)

# ANOVA paired for inv_rfr, inv_ffrado, inv_ffrdobu
anova_test <- aov(inv_rfr ~ inv_ffrado + inv_ffrdobu, data = baseline)
summary(anova_test)
# post hoc
wilcox.test(baseline$inv_rfr, baseline$inv_ffrado, paired = T)
wilcox.test(baseline$inv_rfr, baseline$inv_ffrdobu, paired = T)
wilcox.test(baseline$inv_ffrado, baseline$inv_ffrdobu, paired = T)

# same with inv_ivusrest_mla_ln, inv_ivusado_mla_ln, inv_ivusdobu_mla_ln_any
anova_test <- aov(inv_ivusrest_mla_ln ~ inv_ivusado_mla_ln + inv_ivusdobu_mla_ln_any, data = baseline)
summary(anova_test)
# post hoc
wilcox.test(baseline$inv_ivusrest_mla_ln, baseline$inv_ivusado_mla_ln, paired = T)
wilcox.test(baseline$inv_ivusrest_mla_ln, baseline$inv_ivusdobu_mla_ln_any, paired = T)
wilcox.test(baseline$inv_ivusado_mla_ln, baseline$inv_ivusdobu_mla_ln_any, paired = T)

median(baseline$inv_ivusrest_mla_ln, na.rm = T)
quantile(baseline$inv_ivusrest_mla_ln, c(0.25, 0.75), na.rm = T)
median(baseline$inv_ivusado_mla_ln, na.rm = T)
quantile(baseline$inv_ivusado_mla_ln, c(0.25, 0.75), na.rm = T)
median(baseline$inv_ivusdobu_mla_ln_any, na.rm = T)
quantile(baseline$inv_ivusdobu_mla_ln_any, c(0.25, 0.75), na.rm = T)


############################################################################################################
# roc cross validation
cross <- baseline %>% drop_na(risk_mla_any, all_of(ivus_vars))
cross$risk_mla_any <- factor(cross$risk_mla_any, levels = c(0, 1), labels = c("Class0", "Class1"))

# Create a list of 50 seeds
set.seed(69)
seeds <- sample(1:1000, 100)

results_list <- list()

# Function to train and evaluate models
train_and_evaluate <- function(seed) {
  set.seed(seed)
  
  # Split the data into training (80%) and testing (20%) sets
  train_indices <- createDataPartition(cross$risk_mla_any, p = 0.8, list = FALSE)
  train_data <- cross[train_indices, ]
  test_data <- cross[-train_indices, ]
  
  # Initialize a list to store results for this seed
  seed_results <- list()
  
  # Train logistic regression models and evaluate ROC for each IVUS variable
  for (ivus_var in ivus_vars) {
    formula <- as.formula(paste("risk_mla_any ~", ivus_var))
    model <- glm(formula, data = train_data, family = "binomial")
    
    predictions <- predict(model, test_data, type = "response")

    prediction_long <- simple_logistic_regression(train_data, "risk_mla_any", ivus_var)[[2]]
    roc_curve <- roc(test_data$risk_mla_any, predictions)
    
    best_coords <- tryCatch({
      coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"))
    }, error = function(e) {
      message(sprintf("Error in coords for variable %s: %s", ivus_var, e))
      return(data.frame(threshold = NA, sensitivity = NA, specificity = NA))
    })
    
    if (nrow(best_coords) > 1) {
      best_coords <- best_coords[1,]
    }

    # All other cut-offs and then highest sensitivity with highest spec, and highest spec with highest sens
    coords_list <- coords(roc_curve, "all", ret = c("threshold", "sensitivity", "specificity"))
    coords_list <- coords_list[-c(1, nrow(coords_list)), ]
    best_coords_sens <- coords_list[order(coords_list$sensitivity, decreasing = TRUE),][1,]
    best_coords_spec <- coords_list[order(coords_list$specificity, decreasing = TRUE),][1,]

    thresholds <- c(as.numeric(best_coords$threshold), as.numeric(best_coords_sens$threshold), as.numeric(best_coords_spec$threshold))
    sensitivities <- c(as.numeric(best_coords$sensitivity), as.numeric(best_coords_sens$sensitivity), as.numeric(best_coords_spec$sensitivity))
    specificities <- c(as.numeric(best_coords$specificity), as.numeric(best_coords_sens$specificity), as.numeric(best_coords_spec$specificity))
    
    for (i in 1:length(thresholds)) {
      threshold <- thresholds[i]
      sens_other <- sensitivities[i]
      spec_other <- specificities[i]
      closest_index <- which.min(abs(prediction_long[["risk_mla_any"]] - threshold))
      value <- prediction_long[[ivus_var]][closest_index]

      if (name %in% c("inv_ivusrest_ostial_a_bsa", "inv_ivusrest_imla_bsa", "inv_ivusrest_mla_bsa", "inv_ffrado")) {
        table <- test_data %>%
          mutate(pred_label = ifelse(!!sym(ivus_var) > value, 0, 1)) %>%
          mutate(risk_mla_any = ifelse(risk_mla_any == "Class0", 0, 1)) %>%
          select(pred_label, risk_mla_any) %>%
          drop_na() %>%
          table()
      } else {
        table <- test_data %>%
          mutate(pred_label = ifelse(!!sym(ivus_var) < value, 0, 1)) %>%
          mutate(risk_mla_any = ifelse(risk_mla_any == "Class0", 0, 1)) %>%
          select(pred_label, risk_mla_any) %>%
          drop_na() %>%
          table()
      }

      confusion <- conf_mat(table)
      result <- summary(confusion, event_level = "second")

      sens <- result %>% filter(.metric == "sens") %>% select(.estimate) %>% pull()
      spec <- result %>% filter(.metric == "spec") %>% select(.estimate) %>% pull()
      ppv <- result %>% filter(.metric == "ppv") %>% select(.estimate) %>% pull()
      npv <- result %>% filter(.metric == "npv") %>% select(.estimate) %>% pull()
      acc <- result %>% filter(.metric == "accuracy") %>% select(.estimate) %>% pull()


      # Naming convention ivus_var_1, ivus_var_2, ivus_var_3
      seed_results[[paste(ivus_var, i, sep = "_")]] <- list(
        AUC = as.numeric(auc(roc_curve)),
        Threshold = as.numeric(value),
        Sensitivity = as.numeric(sens),
        Specificity = as.numeric(spec),
        PPV = as.numeric(ppv),
        NPV = as.numeric(npv),
        Accuracy = as.numeric(acc),
        TrueNegatives = as.numeric(table[1, 1]),
        FalseNegatives = as.numeric(table[1, 2]),
        FalsePositives = as.numeric(table[2, 1]),
        TruePositives = as.numeric(table[2, 2]),
        Sensitivity_other = as.numeric(sens_other),
        Specificity_other = as.numeric(spec_other)
      )
    }
  }
  
  return(seed_results)
}

# Run the training and evaluation for each seed and store results
for (seed in seeds) {
  seed_results <- train_and_evaluate(seed)
  results_list[[as.character(seed)]] <- seed_results
}

# Aggregate the results across different seeds
aggregate_results <- function(results_list) {
  aggregate_list <- list()
  
  for (ivus_var in ivus_vars) {
    for (i in 1:3) {
      var_name <- paste(ivus_var, i, sep = "_")
      aucs <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$AUC else NA)
      thresholds <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$Threshold else NA)
      sensitivities <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$Sensitivity else NA)
      specificities <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$Specificity else NA)
      ppvs <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$PPV else NA)
      npvs <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$NPV else NA)
      accuracies <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$Accuracy else NA)
      true_negatives <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$TrueNegatives else NA)
      false_negatives <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$FalseNegatives else NA)
      false_positives <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$FalsePositives else NA)
      true_positives <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$TruePositives else NA)
      sensitivity_others <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$Sensitivity_other else NA)
      specificity_others <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$Specificity_other else NA)
      
      aggregate_list[[var_name]] <- list(
        AUC = mean(aucs, na.rm = TRUE),
        Threshold = mean(thresholds, na.rm = TRUE),
        Sensitivity = mean(sensitivities, na.rm = TRUE),
        Specificity = mean(specificities, na.rm = TRUE),
        PPV = mean(ppvs, na.rm = TRUE),
        NPV = mean(npvs, na.rm = TRUE),
        Accuracy = mean(accuracies, na.rm = TRUE),
        TrueNegatives = mean(true_negatives, na.rm = TRUE),
        FalseNegatives = mean(false_negatives, na.rm = TRUE),
        FalsePositives = mean(false_positives, na.rm = TRUE),
        TruePositives = mean(true_positives, na.rm = TRUE),
        Sensitivity_other = mean(sensitivity_others, na.rm = TRUE),
        Specificity_other = mean(specificity_others, na.rm = TRUE)
      )
    }
  }
  return(aggregate_list)
}

aggregated_results <- aggregate_results(results_list)

# Convert list of lists to a data frame
roc_summary_df <- do.call(rbind, lapply(names(aggregated_results), function(var_name) {
  res <- aggregated_results[[var_name]]
  data.frame(
    Variable = var_name,
    AUC = round(res$AUC, 3),
    Threshold = round(res$Threshold, 2),
    Sensitivity = round(res$Sensitivity * 100, 0),
    Specificity = round(res$Specificity * 100, 0),
    PPV = round(res$PPV * 100, 0),
    NPV = round(res$NPV * 100, 0),
    Accuracy = round(res$Accuracy * 100, 0),
    TrueNegatives = round(res$TrueNegatives, 0),
    FalseNegatives = round(res$FalseNegatives, 0),
    FalsePositives = round(res$FalsePositives, 0),
    TruePositives = round(res$TruePositives, 0),
    Sensitivity_other = round(res$Sensitivity_other * 100, 0),
    Specificity_other = round(res$Specificity_other * 100, 0)
  )
}))

print(roc_summary_df)

write_xlsx(roc_summary_df, "roc_summary_aggregated.xlsx")




########################### FFR0.8 as prediction variable ##################################################
# best cutoffs
# Fit models and get prediction data
mdl_mla_bsa_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_bsa")[[1]]
mdl_mla_ln_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_ln")[[1]]
mdl_ffrado <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ffrado")[[1]]

prediction_data_mla_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_bsa")[[2]]
prediction_data_mla_ln_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_ln")[[2]]
prediction_data_ffrado <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ffrado")[[2]]

data_ivus <- baseline %>% select(ffr_0.8, inv_ivusrest_mla_bsa) %>% drop_na()
roc_mla_ffr <- roc(data_ivus$ffr_0.8, mdl_mla_bsa_ffr$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, inv_ivusrest_mla_ln) %>% drop_na()
roc_mln_ffr <- roc(data_ivus$ffr_0.8, mdl_mla_ln_ffr$fitted.values)
data_ivus <- baseline %>% select(ffr_0.8, inv_ffrado) %>% drop_na()
roc_ffrado <- roc(data_ivus$ffr_0.8, mdl_ffrado$fitted.values)

# Define test sets
test_set <- list(
    coords(roc_mla_ffr), coords(roc_mln_ffr), coords(roc_ffrado)
)

prediction_data <- list(
    prediction_data_mla_ffr, prediction_data_mla_ln_ffr, prediction_data_ffrado
)

variables <- c(
    "Supplemental IVUS rest MLA (BSA)", "Supplemental IVUS rest MLA luminal narrowing", "Supplemental FFR adenosine"
)

baseline <- baseline %>% mutate(
    ffr_0.8 = ifelse(ffr_0.8 == "yes", 1, 0)
)
# Extract thresholds
thresholds <- list()

for (i in 1:length(test_set)) {
  coordinates <- test_set[[i]]
  coordinates <- coordinates[-c(1, nrow(coordinates)), ]
  coordinates <- as_tibble(coordinates)
  
  result <- coordinates %>%
    group_by(sensitivity) %>%
    filter(specificity == max(specificity)) %>%
    select(threshold) %>%
    distinct()
  
  threshold_values <- result$threshold
  thresholds[[variables[i]]] <- threshold_values
}

# Generate confusion matrices and calculate metrics
cutoff_tables <- list()

for (i in 1:length(thresholds)) {
  thresholds_testing <- thresholds[[i]]
  prediction <- prediction_data[[i]]
  name <- names(prediction)[1]
  
  table_thresholds <- tibble()
  for (j in 1:length(thresholds_testing)) {
    threshold <- thresholds_testing[j]
    closest_index <- which.min(abs(prediction$ffr_0.8 - threshold))
    value <- prediction[[name]][closest_index]
    
    if (name %in% c("inv_ivusrest_ostial_a_bsa", "inv_ivusrest_imla_bsa", "inv_ivusrest_mla_bsa", "inv_ffrado")) {
      table <- baseline %>%
        mutate(pred_label = ifelse(!!sym(name) > value, 0, 1)) %>%
        select(pred_label, ffr_0.8) %>%
        drop_na() %>%
        table()
    } else {
      table <- baseline %>%
        mutate(pred_label = ifelse(!!sym(name) < value, 0, 1)) %>%
        select(pred_label, ffr_0.8) %>%
        drop_na() %>%
        table()
    }
    
    confusion <- conf_mat(table)
    result <- summary(confusion, event_level = "second")
    
    sens <- result %>% filter(.metric == "sens") %>% select(.estimate) %>% pull()
    spec <- result %>% filter(.metric == "spec") %>% select(.estimate) %>% pull()
    ppv <- result %>% filter(.metric == "ppv") %>% select(.estimate) %>% pull()
    npv <- result %>% filter(.metric == "npv") %>% select(.estimate) %>% pull()
    acc <- result %>% filter(.metric == "accuracy") %>% select(.estimate) %>% pull()
    
    df <- tibble(
      value = as.numeric(round(value, 2)),
      sensitivity = as.numeric(round(sens * 100, 0)),
      specificity = as.numeric(round(spec * 100, 0)),
      ppv = as.numeric(round(ppv * 100, 0)),
      npv = as.numeric(round(npv * 100, 0)),
      accuracy = as.numeric(round(acc * 100, 0)),
      number_of_true_negatives = as.numeric(table[1, 1]),
      number_of_false_negatives = as.numeric(table[1, 2]),
      number_of_false_positives = as.numeric(table[2, 1]),
      number_of_true_positives = as.numeric(table[2, 2])
    )
    names(df)[1] <- name
    table_thresholds <- bind_rows(table_thresholds, df)
  }
  cutoff_tables[[variables[i]]] <- table_thresholds
}

# Save all as CSV
path <- "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/"
for (i in 1:length(cutoff_tables)) {
  name <- variables[i]
  table <- cutoff_tables[[i]]
  write_csv(table, paste0(path, name, ".csv"))
}

roc_list_mla_ffr <- list(roc_mla_ffr, roc_mln_ffr)
ggroc_mla_ffr <- ggroc(roc_list_mla_ffr, legacy.axes = TRUE) +
    scale_color_manual(values = c("darkred", "darkblue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Minimal lumen") +
    theme_classic() + 
    annotate("text", x = 0.7, y = 0.3, label = paste("AUC MLA =", round(roc_mla_ffr$auc, 2)), color = "darkred") + 
    annotate("text", x = 0.7, y = 0.2, label = paste("AUC MLN =", round(roc_mln_ffr$auc, 2)), color = "darkblue")

ggsave(paste0(output_dir,"/supplemental_roc_mla_ffr.png"), ggroc_mla_ffr, width = 6, height = 5)

############################################################################################################
# roc cross validation for ffr0.8
cross <- baseline %>% drop_na(ffr_0.8, all_of(ivus_vars))
cross$ffr_0.8 <- factor(cross$ffr_0.8, levels = c(0, 1), labels = c("Class0", "Class1"))

# Create a list of 50 seeds
set.seed(69)
seeds <- sample(1:1000, 100)

results_list <- list()

# Function to train and evaluate models
train_and_evaluate <- function(seed) {
  set.seed(seed)
  
  # Split the data into training (80%) and testing (20%) sets
  train_indices <- createDataPartition(cross$ffr_0.8, p = 0.8, list = FALSE)
  train_data <- cross[train_indices, ]
  test_data <- cross[-train_indices, ]
  
  # Initialize a list to store results for this seed
  seed_results <- list()
  
  # Train logistic regression models and evaluate ROC for each IVUS variable
  for (ivus_var in ivus_vars) {
    formula <- as.formula(paste("ffr_0.8 ~", ivus_var))
    model <- glm(formula, data = train_data, family = "binomial")
    
    predictions <- predict(model, test_data, type = "response")

    prediction_long <- simple_logistic_regression(train_data, "ffr_0.8", ivus_var)[[2]]
    roc_curve <- roc(test_data$ffr_0.8, predictions)
    
    best_coords <- tryCatch({
      coords(roc_curve, "best", ret = c("threshold", "sensitivity", "specificity"))
    }, error = function(e) {
      message(sprintf("Error in coords for variable %s: %s", ivus_var, e))
      return(data.frame(threshold = NA, sensitivity = NA, specificity = NA))
    })
    
    if (nrow(best_coords) > 1) {
      best_coords <- best_coords[1,]
    }

    # All other cut-offs and then highest sensitivity with highest spec, and highest spec with highest sens
    coords_list <- coords(roc_curve, "all", ret = c("threshold", "sensitivity", "specificity"))
    coords_list <- coords_list[-c(1, nrow(coords_list)), ]
    best_coords_sens <- coords_list[order(coords_list$sensitivity, decreasing = TRUE),][1,]
    best_coords_spec <- coords_list[order(coords_list$specificity, decreasing = TRUE),][1,]

    thresholds <- c(as.numeric(best_coords$threshold), as.numeric(best_coords_sens$threshold), as.numeric(best_coords_spec$threshold))
    sensitivities <- c(as.numeric(best_coords$sensitivity), as.numeric(best_coords_sens$sensitivity), as.numeric(best_coords_spec$sensitivity))
    specificities <- c(as.numeric(best_coords$specificity), as.numeric(best_coords_sens$specificity), as.numeric(best_coords_spec$specificity))
    
    for (i in 1:length(thresholds)) {
      threshold <- thresholds[i]
      sens_other <- sensitivities[i]
      spec_other <- specificities[i]
      closest_index <- which.min(abs(prediction_long[["ffr_0.8"]] - threshold))
      value <- prediction_long[[ivus_var]][closest_index]

      if (name %in% c("inv_ivusrest_ostial_a_bsa", "inv_ivusrest_imla_bsa", "inv_ivusrest_mla_bsa", "inv_ffrado")) {
        table <- test_data %>%
          mutate(pred_label = ifelse(!!sym(ivus_var) > value, 0, 1)) %>%
          mutate(ffr_0.8 = ifelse(ffr_0.8 == "Class0", 0, 1)) %>%
          select(pred_label, ffr_0.8) %>%
          drop_na() %>%
          table()
      } else {
        table <- test_data %>%
          mutate(pred_label = ifelse(!!sym(ivus_var) < value, 0, 1)) %>%
          mutate(ffr_0.8 = ifelse(ffr_0.8 == "Class0", 0, 1)) %>%
          select(pred_label, ffr_0.8) %>%
          drop_na() %>%
          table()
      }

      confusion <- conf_mat(table)
      result <- summary(confusion, event_level = "second")

      sens <- result %>% filter(.metric == "sens") %>% select(.estimate) %>% pull()
      spec <- result %>% filter(.metric == "spec") %>% select(.estimate) %>% pull()
      ppv <- result %>% filter(.metric == "ppv") %>% select(.estimate) %>% pull()
      npv <- result %>% filter(.metric == "npv") %>% select(.estimate) %>% pull()
      acc <- result %>% filter(.metric == "accuracy") %>% select(.estimate) %>% pull()


      # Naming convention ivus_var_1, ivus_var_2, ivus_var_3
      seed_results[[paste(ivus_var, i, sep = "_")]] <- list(
        AUC = as.numeric(auc(roc_curve)),
        Threshold = as.numeric(value),
        Sensitivity = as.numeric(sens),
        Specificity = as.numeric(spec),
        PPV = as.numeric(ppv),
        NPV = as.numeric(npv),
        Accuracy = as.numeric(acc),
        TrueNegatives = as.numeric(table[1, 1]),
        FalseNegatives = as.numeric(table[1, 2]),
        FalsePositives = as.numeric(table[2, 1]),
        TruePositives = as.numeric(table[2, 2]),
        Sensitivity_other = as.numeric(sens_other),
        Specificity_other = as.numeric(spec_other)
      )
    }
  }
  
  return(seed_results)
}

# Run the training and evaluation for each seed and store results
for (seed in seeds) {
  seed_results <- train_and_evaluate(seed)
  results_list[[as.character(seed)]] <- seed_results
}

# Aggregate the results across different seeds
aggregate_results <- function(results_list) {
  aggregate_list <- list()
  
  for (ivus_var in ivus_vars) {
    for (i in 1:3) {
      var_name <- paste(ivus_var, i, sep = "_")
      aucs <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$AUC else NA)
      thresholds <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$Threshold else NA)
      sensitivities <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$Sensitivity else NA)
      specificities <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$Specificity else NA)
      ppvs <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$PPV else NA)
      npvs <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$NPV else NA)
      accuracies <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$Accuracy else NA)
      true_negatives <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$TrueNegatives else NA)
      false_negatives <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$FalseNegatives else NA)
      false_positives <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$FalsePositives else NA)
      true_positives <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$TruePositives else NA)
      sensitivity_others <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$Sensitivity_other else NA)
      specificity_others <- sapply(results_list, function(x) if(!is.null(x[[var_name]])) x[[var_name]]$Specificity_other else NA)
      
      aggregate_list[[var_name]] <- list(
        AUC = mean(aucs, na.rm = TRUE),
        Threshold = mean(thresholds, na.rm = TRUE),
        Sensitivity = mean(sensitivities, na.rm = TRUE),
        Specificity = mean(specificities, na.rm = TRUE),
        PPV = mean(ppvs, na.rm = TRUE),
        NPV = mean(npvs, na.rm = TRUE),
        Accuracy = mean(accuracies, na.rm = TRUE),
        TrueNegatives = mean(true_negatives, na.rm = TRUE),
        FalseNegatives = mean(false_negatives, na.rm = TRUE),
        FalsePositives = mean(false_positives, na.rm = TRUE),
        TruePositives = mean(true_positives, na.rm = TRUE),
        Sensitivity_other = mean(sensitivity_others, na.rm = TRUE),
        Specificity_other = mean(specificity_others, na.rm = TRUE)
      )
    }
  }
  return(aggregate_list)
}

aggregated_results <- aggregate_results(results_list)

# Convert list of lists to a data frame
roc_summary_df_ffr0.8 <- do.call(rbind, lapply(names(aggregated_results), function(var_name) {
  res <- aggregated_results[[var_name]]
  data.frame(
    Variable = var_name,
    AUC = round(res$AUC, 3),
    Threshold = round(res$Threshold, 2),
    Sensitivity = round(res$Sensitivity * 100, 0),
    Specificity = round(res$Specificity * 100, 0),
    PPV = round(res$PPV * 100, 0),
    NPV = round(res$NPV * 100, 0),
    Accuracy = round(res$Accuracy * 100, 0),
    TrueNegatives = round(res$TrueNegatives, 0),
    FalseNegatives = round(res$FalseNegatives, 0),
    FalsePositives = round(res$FalsePositives, 0),
    TruePositives = round(res$TruePositives, 0),
    Sensitivity_other = round(res$Sensitivity_other * 100, 0),
    Specificity_other = round(res$Specificity_other * 100, 0)
  )
}))

print(roc_summary_df_ffr0.8)

write_xlsx(roc_summary_df_ffr0.8, "roc_summary_ffr0.8_aggregated.xlsx")


############################################################################################################
# cut-off 0.8 for ffr_ado
# Fit models and get prediction data
baseline <- baseline %>% mutate(
  ffr_ado0.8 = ifelse(inv_ffrado <= 0.8, 1, 0)
)

table <- baseline %>%
  select(ffr_ado0.8, ffr_0.8) %>%
  drop_na() %>%
  table()

confusion <- conf_mat(table)
result <- summary(confusion, event_level = "second")

############################################################################################################
table <- baseline %>%
  select(ffr_ado0.8, risk_mla_any) %>%
  drop_na() %>%
  table()

confusion <- conf_mat(table)
result <- summary(confusion, event_level = "second")

############################################################################################################
baseline <- baseline %>% mutate(
  mla_ln50 = ifelse(inv_ivusrest_mla_ln >= 50, 1, 0)
)

table <- baseline %>%
  select(mla_ln50, risk_mla_any) %>%
  drop_na() %>%
  table()

confusion <- conf_mat(table)
result <- summary(confusion, event_level = "second")

############################################################################################################
table <- baseline %>%
  select(mla_ln50, ffr_0.8) %>%
  drop_na() %>%
  table()

confusion <- conf_mat(table)
result <- summary(confusion, event_level = "second")

# roc curves
roc_invffrado <- roc(baseline$ffr_0.8, baseline$inv_ffrado)
ggroc(roc_invffrado, legacy.axes = TRUE, color = "darkblue") +
  theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ggtitle("FFR adenosine") +
  theme_classic() + 
  annotate("text", x = 0.7, y = 0.3, label = paste("AUC FFR adenosine =", round(roc_invffrado$auc, 2)))

roc_invffrado_risk_mla <- roc(baseline$risk_mla_any, baseline$inv_ffrado)
ggroc(roc_invffrado_risk_mla, legacy.axes = TRUE, color = "darkblue") +
  theme(legend.position = "none") +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  ggtitle("FFR adenosine") +
  theme_classic() + 
  annotate("text", x = 0.7, y = 0.3, label = paste("AUC FFR adenosine =", round(roc_invffrado_risk_mla$auc, 2)))


# forest plot
# for hemodynamic and anatomic relevance as predictor
forest_df_bsa <- data.frame(
    var = c(
            "Minimal lumen area (BSA) [mm2]", 
            "Minimal lumen narrowing [%]", "FFR adenosine"),
    odds_ratio = c(
            exp(mdl_mla_bsa_ffr$coefficients[2]), 
            exp(mdl_mla_ln_ffr$coefficients[2]),
            exp(mdl_ffrado$coefficients[2])),
    lower_ci = c(
            exp(confint(mdl_mla_bsa_ffr))[2], 
            exp(confint(mdl_mla_ln_ffr))[2], 
            exp(confint(mdl_ffrado))[2]), 
    upper_ci = c(
            exp(confint(mdl_mla_bsa_ffr))[4], 
            exp(confint(mdl_mla_ln_ffr))[4],  
            exp(confint(mdl_ffrado))[4])
)

forest_df_bsa$var <- factor(forest_df_bsa$var, levels = rev(forest_df_bsa$var))

forest_cut <- forest_df_bsa %>%
    mutate(
        upper_ci = ifelse(upper_ci > 6, 6, upper_ci)
    )

forest_risk_mla_any <- ggplot(data = forest_cut, aes(y = var, x = odds_ratio, xmin = lower_ci, xmax = upper_ci)) +
    geom_point() +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, 2, 0.5), labels = parse(text = seq(0, 2, 0.5))) +
    geom_vline(xintercept = 2, linetype = "solid", color = "white", lwd = 2) +
    theme_classic()

ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/forest_risk_mla_any.png", forest_risk_mla_any, width = 5, height = 2)

mdl_mla_bsa_supp <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_bsa")[[1]]
mdl_ln_supp <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_ln")[[1]]
mdl_ffrado_supp <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ffrado")[[1]]

forest_df_supp <- data.frame(
    var = c(
            "Minimal lumen area (BSA) [mm2]", 
            "Minimal lumen narrowing [%]", "FFR adenosine"),
    odds_ratio = c(
            exp(mdl_mla_bsa_supp$coefficients[2]), 
            exp(mdl_ln_supp$coefficients[2]),
            exp(mdl_ffrado_supp$coefficients[2])),
    lower_ci = c(
            exp(confint(mdl_mla_bsa_supp))[2], 
            exp(confint(mdl_ln_supp))[2], 
            exp(confint(mdl_ffrado_supp))[2]), 
    upper_ci = c(
            exp(confint(mdl_mla_bsa_supp))[4], 
            exp(confint(mdl_ln_supp))[4],  
            exp(confint(mdl_ffrado_supp))[4])
)

forest_df_supp$var <- factor(forest_df_supp$var, levels = rev(forest_df_supp$var))

forest_cut_supp <- forest_df_supp %>%
    mutate(
        upper_ci = ifelse(upper_ci > 6, 6, upper_ci)
    )

forest_risk_ffr_supp <- ggplot(data = forest_cut_supp, aes(y = var, x = odds_ratio, xmin = lower_ci, xmax = upper_ci)) +
    geom_point() +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, 2, 0.5), labels = parse(text = seq(0, 2, 0.5))) +
    geom_vline(xintercept = 2, linetype = "solid", color = "white", lwd = 2) +
    theme_classic()
  
ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/forest_risk_ffr_supp.png", forest_risk_ffr_supp, width = 5, height = 2)


# # create one dataframe by combining the rows with inv_ivusrest_mla_ln3, inv_ivusrest_mla_bsa_3 from all roc_summary_df
# roc_summary_df_10 <- roc_summary_df_10 %>% filter(Variable %in% c("inv_ivusrest_mla_ln_3", "inv_ivusrest_mla_bsa_3"))
# roc_summary_df_20 <- roc_summary_df_20 %>% filter(Variable %in% c("inv_ivusrest_mla_ln_3", "inv_ivusrest_mla_bsa_3"))
# roc_summary_df_50 <- roc_summary_df_50 %>% filter(Variable %in% c("inv_ivusrest_mla_ln_3", "inv_ivusrest_mla_bsa_3"))
# roc_summary_df_75 <- roc_summary_df_75 %>% filter(Variable %in% c("inv_ivusrest_mla_ln_3", "inv_ivusrest_mla_bsa_3"))
# roc_summary_df_100 <- roc_summary_df_100 %>% filter(Variable %in% c("inv_ivusrest_mla_ln_3", "inv_ivusrest_mla_bsa_3"))
# roc_summary_df_150 <- roc_summary_df_150 %>% filter(Variable %in% c("inv_ivusrest_mla_ln_3", "inv_ivusrest_mla_bsa_3"))
# roc_summary_df_175 <- roc_summary_df_175 %>% filter(Variable %in% c("inv_ivusrest_mla_ln_3", "inv_ivusrest_mla_bsa_3"))

# df <- tibble(
#   iteration = c(10, 20, 50, 75, 100, 150, 175),
#   AUC_mla_ln = c(roc_summary_df_10[1, "AUC"], roc_summary_df_20[1, "AUC"], roc_summary_df_50[1, "AUC"], roc_summary_df_75[1, "AUC"], roc_summary_df_100[1, "AUC"], roc_summary_df_150[1, "AUC"], roc_summary_df_175[1, "AUC"]),
#   AUC_mla_bsa = c(roc_summary_df_10[2, "AUC"], roc_summary_df_20[2, "AUC"], roc_summary_df_50[2, "AUC"], roc_summary_df_75[2, "AUC"],roc_summary_df_100[2, "AUC"], roc_summary_df_150[2, "AUC"], roc_summary_df_175[2, "AUC"]),
# )

# ggplot(df, aes(x = iteration)) +
#   geom_line(aes(y = AUC_mla_ln, color = "MLA Luminal Narrowing")) +
#   geom_line(aes(y = AUC_mla_bsa, color = "MLA BSA")) +
#   scale_color_manual(values = c("MLA Luminal Narrowing" = "darkred", "MLA BSA" = "darkblue")) +
#   theme_classic() +
#   labs(title = "AUC comparison for different iterations",
#        x = "Iteration",
#        y = "AUC") +
#   theme(legend.position = "top")

# ggsave(paste0(output_dir,"/roc_summary_comparison.png"), last_plot(), width = 6, height = 5)

# baseline %>% filter(inv_ivusrest_mla_ln <70 & inv_ivusdobu_mla_ln_any >=70) %>% select(patient_id, inv_ivusrest_mla_ln, inv_ivusdobu_mla_ln_any, inv_ffrdobu) %>% print(n = 50)

# ggroc_mla <- ggroc(roc_mla, legacy.axes = TRUE, color = "darkblue") +
#   theme(legend.position = "none") +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#   ggtitle("Minimal lumen area (BSA adjusted)") +
#   theme_classic() + 
#   annotate("text", x = 0.7, y = 0.3, label = paste("AUC MLA (BSA adjusted) =", round(roc_mla$auc, 2))) 

# ggroc_mln <- ggroc(roc_mla_ln, legacy.axes = TRUE, color = "darkblue") +
#   theme(legend.position = "none") +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
#   ggtitle("Minimal lumen narrowing") +
#   theme_classic() + 
#   annotate("text", x = 0.7, y = 0.3, label = paste("AUC MLN =", round(roc_mla_ln$auc, 2)))

# ggsave(paste0(output_dir,"/roc_mla.png"), ggroc_mla, width = 6, height = 5)
# ggsave(paste0(output_dir,"/roc_mla_ln.png"), ggroc_mln, width = 6, height = 5)


# baseline %>%
#   filter(is.na(inv_ivusrest_mla)) %>%
#   select(caa_origin___0, caa_origin___1, caa_origin___2, caa_origin___3, caa_origin___4, caa_origin___5, caa_origin___6) %>%
#   mutate(across(everything(), ~ as.numeric(.) - 1)) %>%
#   summarise_all(sum, na.rm = TRUE) %>%
#   View()