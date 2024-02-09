library(tidyverse)
library(ggplot2)
library(readxl)
library(yaml)

library(DescTools)
library(infer)

## DATA LOADING AND PREPPING ################################################################################
yaml <- yaml.load_file(file.path(getwd(), "config.yaml")) # only works after activating venv

baseline <- read_excel(paste0(yaml$first_analysis$output_dir, "/baseline.xlsx"))
# remove the first row, which is test data
baseline <- baseline[-1, ]

patterns <- list(
  hx = "^hx_",
  sports = "^sports_",
  ccta = "^ccta_",
  inv = "^inv_",
  caa = "^caa_"
)

filtered_columns <- list()
for (pattern_name in names(patterns)) {
  filtered_columns[[pattern_name]] <- grep(patterns[[pattern_name]], names(baseline), value = TRUE)
  if (pattern_name == "ccta") {
    filtered_columns[[pattern_name]] <- filtered_columns[[pattern_name]][!grepl("_postop", filtered_columns[[pattern_name]])]
  }
}

baseline <- baseline %>% select(all_of(yaml$demographics$dict_variables), all_of(unlist(filtered_columns, use.names = FALSE)))

for (i in seq_along(yaml$demographics$dict_rename)) {
    new_col_name <- names(yaml$demographics$dict_rename)[i]
    old_col_name <- yaml$demographics$dict_rename[[i]]

    baseline <- baseline %>% rename(!!new_col_name := !!old_col_name)
}

## MUTATE ##########################################################################################################
baseline <- baseline %>% 
  mutate(
    sports_classification = case_when(
      sports_static_component == 0 & sports_dynamic_component == 0 ~ 1,
      (sports_static_component == 1 & sports_dynamic_component == 0) | (sports_static_component == 0 & sports_dynamic_component ==1) ~ 2,
      (sports_static_component == 2 & sports_dynamic_component == 0) | (sports_static_component == 1 & sports_dynamic_component == 1) | (sports_static_component == 0 & sports_dynamic_component == 2) ~ 3,
      (sports_static_component == 2 & sports_dynamic_component == 1) | (sports_static_component == 1 & sports_dynamic_component ==2) ~ 4,
      sports_static_component == 2 & sports_dynamic_component == 2 ~ 5,
    ),
    sports_classification_ch = case_when(
      sports_static_component_ch == 0 & sports_dynamic_comp_ch == 0 ~ 1,
      (sports_static_component_ch == 1 & sports_dynamic_comp_ch == 0) | (sports_static_component_ch == 0 & sports_dynamic_comp_ch ==1) ~ 2,
      (sports_static_component_ch == 2 & sports_dynamic_comp_ch == 0) | (sports_static_component_ch == 1 & sports_dynamic_comp_ch == 1) | (sports_static_component_ch == 0 & sports_dynamic_comp_ch == 2) ~ 3,
      (sports_static_component_ch == 2 & sports_dynamic_comp_ch == 1) | (sports_static_component_ch == 1 & sports_dynamic_comp_ch ==2) ~ 4,
      sports_static_component_ch == 2 & sports_dynamic_comp_ch == 2 ~ 5,
    ),
    ffr_0.8 = ifelse(inv_ffrdobu <= 0.8, 1, 0),
    ffr_0.81 = ifelse(inv_ffrdobu <= 0.81, 1, 0)
  )

## FACTORIZE ##########################################################################################################
baseline$sex_0_male_1_female <- factor(baseline$sex_0_male_1_female, levels = c(0, 1), labels = c("male", "female"))
baseline$caa_modality <- factor(baseline$caa_modality, levels = c(0, 1, 2, 3, 4, 5), labels = c("ccta", "angiography", "echocardiography", "cmr", "surgery", "other"))
baseline$hx_sym_ap_ccs <- factor(baseline$hx_sym_ap_ccs, levels = 1:4, ordered = TRUE, labels = c("I", "II", "III", "IV"))
baseline$hx_sym_dysp_nyha <- factor(baseline$hx_sym_dysp_nyha, levels = 1:4, ordered = TRUE, labels = c("I", "II", "III", "IV"))
baseline$sports_level <- factor(baseline$sports_level, levels = 1:3, ordered = TRUE, labels = c("recreational", "competitive", "elite"))
baseline$sports_static_component <- factor(baseline$sports_static_component, levels = 0:2, ordered = TRUE, labels = c("I", "II", "III"))
baseline$sports_dynamic_component <- factor(baseline$sports_dynamic_component, levels = 0:2, ordered = TRUE, labels = c("A", "B", "C"))
baseline$sports_classification <- factor(baseline$sports_classification, levels = 1:5, ordered = TRUE, labels = c("very low", "low", "moderate", "high", "very high"))
baseline$caa_malignancy <- factor(baseline$caa_malignancy, levels = c(0, 1), labels = c("benign", "malign"))
binary_vars <- names(baseline)[
  sapply(names(baseline), function(var_name) {
    all(baseline[[var_name]] %in% c(0, 1, NA))
  })
]
for (var in binary_vars) {
  baseline[[var]] <- factor(baseline[[var]], levels = c(0, 1), labels = c("no", "yes"))
}

saveRDS(baseline, file = paste0(yaml$demographics$output_dir_data, "/baseline.rds"))
############################################################################################################
# Two splits: one for benign and malign and the other for ffr above and below equal 0.8
# check if contains any NA --> shouldn't be the case, if yes check database
baseline %>% select(record_id, caa_malignancy) %>% filter(is.na(caa_malignancy))
baseline %>% select(record_id, inv_ffrdobu, caa_malignancy) %>% filter(is.na(inv_ffrdobu) & caa_malignancy == 1)

## CREATE RESULTS DATAFRAME ##########################################################################################################
# initialize the dataframe
dataframe_calculations <- data.frame(variable = character(),
                                      type_numeric = numeric(),
                                      n_all = numeric(),
                                      per_all = numeric(),
                                      mean_all = numeric(), 
                                      sd_all = numeric(),
                                      median_all = numeric(), 
                                      iqr_all = numeric(),
                                      n_group1 = numeric(),
                                      per_group1 = numeric(),
                                      mean_group1 = numeric(),
                                      sd_group1 = numeric(),
                                      median_group1 = numeric(),
                                      iqr_group1 = numeric(),
                                      n_group2 = numeric(),
                                      per_group2 = numeric(),
                                      mean_group2 = numeric(), 
                                      sd_group2 = numeric(),
                                      median_group2 = numeric(), 
                                      iqr_group2 = numeric(),
                                      stringsAsFactors = FALSE,
                                      row.names = TRUE)

process_factors <- function(data, var) {
  if (length(levels(data[[var]])) == 2) {
    counts <- as.integer(table(data[[var]])[2])
    prop <- counts / sum(table(data[[var]]))
    return(c(counts, prop))
  }
  else {
    desc_df <- Desc(data[[var]], plot = FALSE)
    counts <- desc_df[1][[1]]$freq["freq"]
    prop <- desc_df[1][[1]]$freq["perc"]
    names <- desc_df[1][[1]]$freq["level"]
    result_df <- data.frame(levels = names, counts = counts, prop = prop)
    return(result_df)
  }
}

process_numerical <- function(data, var) {
  n <- sum(!is.na(data[[var]]))
  per <- n / nrow(data)
  mean <- mean(data[[var]], na.rm = TRUE)
  sd <- sd(data[[var]], na.rm = TRUE)
  median <- median(data[[var]], na.rm = TRUE)
  q25 <- as.numeric(quantile(data[[var]], 0.25, na.rm = TRUE))
  q75 <- as.numeric(quantile(data[[var]], 0.75, na.rm = TRUE))
  iqr <- q75 - q25

  return(c(n, per, mean, sd, median, iqr))
}

fill_dataframe <- function(data, var_to_group_by, df) {
  group1 <- data %>% filter(!!sym(var_to_group_by) == levels(!!sym(var_to_group_by))[1])
  group2 <- data %>% filter(!!sym(var_to_group_by) == levels(!!sym(var_to_group_by))[2])
  
  for (var in names(data)){
    if (is.numeric(data[[var]])) {
      num_data <- process_numerical(data, var)
      num_group1 <- process_numerical(group1, var)
      num_group2 <- process_numerical(group2, var)
      df <- rbind(df, c(variable = var, type_numeric = 1, n_all = num_data[1], per_all = num_data[2], 
      mean_all = num_data[3], sd_all = num_data[4], median_all = num_data[5], 
      iqr_all = num_data[6], n_group1 = num_group1[1], 
      per_group1 = num_group1[2], mean_group1 = num_group1[3], sd_group1 = num_group1[4], 
      median_group1 = num_group1[5], iqr_group1 = num_group1[6]
      , n_group2 = num_group2[1], per_group2 = num_group2[2],
      mean_group2 = num_group2[3], sd_group2 = num_group2[4], median_group2 = num_group2[5],
      iqr_group2 = num_group2[6]))
    }
    else if (is.factor(data[[var]])) {
      factor_all <- process_factors(data, var)
      factor_group1 <- process_factors(group1, var)
      factor_group2 <- process_factors(group2, var)
      if (length(factor_all) == 2) {
        df <- rbind(df, c(variable = var, type_numeric = 0, n_all = factor_all[1], per_all = factor_all[2], 
        mean_all = NA, sd_all = NA, median_all = NA, 
        iqr_all = NA, n_group1 = factor_group1[1], 
        per_group1 = factor_group1[2], mean_group1 = NA, sd_group1 = NA, 
        median_group1 = NA, iqr_group1 = NA,
        n_group2 = factor_group2[1], per_group2 = factor_group2[2],
        mean_group2 = NA, sd_group2 = NA, median_group2 = NA,
        iqr_group2 = NA))
      }
      else {
        for (i in seq_len(nrow(factor_all))) {
          df <- rbind(df, c(variable = paste0(var, "_",factor_all[i, 1]), type_numeric = 0, 
          n_all = factor_all[i, 2], per_all = factor_all[i, 3], 
          mean_all = NA, sd_all = NA, median_all = NA, 
          iqr_all = NA, n_group1 = factor_group1[i, 2], 
          per_group1 = factor_group1[i, 3], mean_group1 = NA, sd_group1 = NA, 
          median_group1 = NA, iqr_group1 = NA,
          n_group2 = factor_group1[i, 2], per_group2 = factor_group1[i, 3],
          mean_group2 = NA, sd_group2 = NA, median_group2 = NA,
          iqr_group2 = NA))
        }
      }
    }  
  }
  suffix_group1 <- paste0("_", levels(data[[var_to_group_by]])[1])
  suffix_group2 <- paste0("_", levels(data[[var_to_group_by]])[2])
  
  colnames(df) <- c("Variable", "Type_Numeric", "N_All", "Per_All", 
                    "Mean_All", "SD_All", "Median_All", "IQR_All", 
                    paste0("N", suffix_group1), paste0("Per", suffix_group1),
                    paste0("Mean", suffix_group1), paste0("SD", suffix_group1),
                    paste0("Median", suffix_group1), paste0("IQR", suffix_group1),
                    paste0("N", suffix_group2), paste0("Per", suffix_group2),
                    paste0("Mean", suffix_group2), paste0("SD", suffix_group2),
                    paste0("Median", suffix_group2), paste0("IQR", suffix_group2))
  return(df)
}

dataframe_malignancy <- fill_dataframe(baseline, "caa_malignancy", dataframe_calculations)
dataframe_ffr0.8 <- fill_dataframe(baseline, "ffr_0.8", dataframe_calculations)
dataframe_ffr0.81 <- fill_dataframe(baseline, "ffr_0.81", dataframe_calculations)

## ADD STATISTICAL TESTS ##########################################################################################
diff_means_simulation <- function(data, response_var, predictor_var, reps = 5000, direction = "two sided") {
  data <- data %>% filter(!is.na(!!sym(predictor_var)))
  formula <- as.formula(paste(response_var, "~", predictor_var))
  null_distn <- data %>%
    specify(formula) %>%
    hypothesize(null = "independence") %>%
    generate(reps = reps, type = "permute") %>%
    calculate(
      stat = "diff in means", 
      order = c(levels(data[[predictor_var]])[1], levels(data[[predictor_var]])[2]))
  obs_stat <- data %>%
    specify(formula) %>%
    calculate(
      stat = "diff in means", 
      order = c(levels(data[[predictor_var]])[1], levels(data[[predictor_var]])[2]))
  p_value <- get_p_value(null_distn, obs_stat, direction = direction)
  return(p_value)
}

diff_props_simulation <- function(data, response_var, predictor_var, reps = 5000, direction = "two sided") {
  data <- data %>% filter(!is.na(!!sym(predictor_var)))
  formula <- as.formula(paste(response_var, "~", predictor_var))

  response_levels <- levels(data[[response_var]])
  null_distn <- data %>%
    specify(formula, success = response_levels[2]) %>%
    hypothesize(null = "independence") %>%
    generate(reps = reps, type = "permute") %>%
    calculate(
      stat = "diff in props",
      order = c(levels(data[[predictor_var]])[2], levels(data[[predictor_var]])[1]))
  obs_stat <- data %>%
    specify(formula, success = response_levels[2]) %>%
    calculate(stat = "diff in props",
      order = c(levels(data[[predictor_var]])[2], levels(data[[predictor_var]])[1]))
  p_value <- get_p_value(null_distn, obs_stat, direction = direction)
  return(p_value)
}

chisq_simulation <- function(data, response_var, predictor_var, reps = 5000, direction = "greater") {
  data <- data %>% filter(!is.na(!!sym(predictor_var)))
  formula <- as.formula(paste(response_var, "~", predictor_var))
  null_distn <- data %>%
    specify(formula) %>%
    hypothesize(null = "independence") %>%
    generate(reps = reps, type = "permute") %>%
    calculate(stat = "Chisq")
  obs_stat <- data %>%
    specify(formula) %>%
    calculate(stat = "Chisq")
  p_value <- get_p_value(null_distn, obs_stat, direction = direction)
  return(p_value)
}

statistics_dataframe <- function(data, df, predictor_var, reps = 5000, direction = "two sided") {
  group1 <- data %>% filter(!!sym(predictor_var) == levels(!!sym(predictor_var))[1])
  group2 <- data %>% filter(!!sym(predictor_var) == levels(!!sym(predictor_var))[2])
  df_stats <- data.frame(variable = character(), p_norm_all = numeric(), p_norm_group1 = numeric(), 
                        p_norm_group2 = numeric(), p_value_t = numeric(), p_value_wilcox = numeric(), 
                        p_value_prop = numeric(), p_value_chi = numeric(), p_value_sim_mean = numeric(), 
                        p_value_sim_prop = numeric(), p_value_sim_chi = numeric(), stringsAsFactors = FALSE)
  for (var in names(data)){
    tryCatch({
      if (is.numeric(data[[var]])) {
        p_sim <- diff_means_simulation(data, var, predictor_var, reps, direction)
        p_all <- shapiro.test(data[[var]])$p.value
        p_group1 <- shapiro.test(group1[[var]])$p.value
        p_group2 <- shapiro.test(group2[[var]])$p.value
        p_t <- t.test(data[[var]] ~ data[[predictor_var]], alternative = "two.sided")$p.value
        p_wilcox <- wilcox.test(data[[var]] ~ data[[predictor_var]], alternative = "two.sided")$p.value
        df_stats <- rbind(df_stats, c(variable = var, p_norm_all = p_all, p_norm_group1 = p_group1, 
                                      p_norm_group2 = p_group2, p_value_t = p_t, p_value_wilcox = p_wilcox, 
                                      p_value_prop = NA, p_value_chi = NA, p_value_sim_mean = as.numeric(p_sim), 
                                      p_value_sim_prop = NA, p_value_sim_chi = NA))
      }
      else if (is.factor(data[[var]])) {
        p_values_prop <- c()
        p_values_chi <- c()
        if (length(levels(data[[var]])) == 2) {
          p_sim_prop <- diff_props_simulation(data, var, predictor_var, reps, direction)
          p_prop <- prop.test(table(data[[var]], data[[predictor_var]]))$p.value
          p_values_prop <- c(as.numeric(p_sim_prop), p_prop)
        } else if (length(levels(data[[var]])) > 2) {  # Perform chi-square test
          p_sim_chi <- chisq_simulation(data, var, predictor_var, reps, direction)
          p_chi <- chisq.test(table(data[[var]], data[[predictor_var]]))$p.value
          p_values_chi <- c(as.numeric(p_sim_chi), p_chi)
        }
        df_stats <- rbind(df_stats, c(variable = var, p_norm_all = NA, p_norm_group1 = NA, 
                                      p_norm_group2 = NA, p_value_t = NA, p_value_wilcox = NA, 
                                      p_value_prop = p_values_prop[2], p_value_chi = p_values_chi[2], 
                                      p_value_sim_mean = NA, p_value_sim_prop = p_values_prop[1], 
                                      p_value_sim_chi = p_values_chi[1]))
      }
      else {
      print(paste("Skipping", var, "as it has only one level."))
      df_stats <- rbind(df_stats, c(variable = var, p_norm_all = NA, p_norm_group1 = NA, 
                                  p_norm_group2 = NA, p_value_t = NA, p_value_wilcox = NA, 
                                  p_value_prop = NA, p_value_chi = NA, p_value_sim_mean = NA, 
                                  p_value_sim_prop = NA, p_value_sim_chi = NA))
    }
    }, error = function(e) {
      cat("Error occurred for variable:", var, "\n")
      message("Error message:", e$message, "\n")
      df_stats <- rbind(df_stats, c(variable = var, p_norm_all = NA, p_norm_group1 = NA, 
                                  p_norm_group2 = NA, p_value_t = NA, p_value_wilcox = NA, 
                                  p_value_prop = NA, p_value_chi = NA, p_value_sim_mean = NA, 
                                  p_value_sim_prop = NA, p_value_sim_chi = NA))
    })
  }
  colnames(df_stats) <- c("Variable", "P_Norm_All", "P_Norm_Group1", "P_Norm_Group2", 
                          "P_Value_T", "P_Value_Wilcox", "P_Value_Prop", "P_Value_Chi", 
                          "P_Value_Sim_Mean", "P_Value_Sim_Prop", "P_Value_Sim_Chi")
  return(df_stats)
}

test2 <- statistics_dataframe(baseline, dataframe_malignancy, "caa_malignancy")