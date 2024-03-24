library(tidyverse)
library(ggplot2)
library(readxl)
library(yaml)
library(writexl)
library(DescTools)
library(infer)

## DATA LOADING AND PREPPING ################################################################################
yaml <- yaml.load_file(file.path(getwd(), "config.yaml")) # only works after activating venv

baseline <- read_excel(paste0(yaml$first_analysis$output_dir, "/baseline.xlsx"))
# remove the first row, which is test data
baseline <- baseline[-1, ]

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
    ffr_0.81 = ifelse(inv_ffrdobu <= 0.81, 1, 0),
    slo_cheezum = ccta_ostial_w / (2 * sqrt(ccta_dist_a / pi))
  )

## FACTORIZE ##########################################################################################################
baseline$sex_0_male_1_female <- factor(baseline$sex_0_male_1_female, levels = c(0, 1), labels = c("male", "female"))
baseline$caa_modality <- factor(baseline$caa_modality, levels = c(0, 1, 2, 3, 4, 5), labels = c("ccta", "angiography", "echocardiography", "cmr", "surgery", "other"))
baseline$hx_sym_ap_ccs <- factor(baseline$hx_sym_ap_ccs, levels = 1:4, ordered = TRUE, labels = c("I", "II", "III", "IV"))
baseline$hx_sym_dysp_nyha <- factor(baseline$hx_sym_dysp_nyha, levels = 1:4, ordered = TRUE, labels = c("I", "II", "III", "IV"))
baseline$sports_level <- factor(baseline$sports_level, levels = 1:3, ordered = TRUE, labels = c("recreational", "competitive", "elite"))
baseline$sports_level_ch <- factor(baseline$sports_level_ch, levels = 1:3, ordered = TRUE, labels = c("recreational", "competitive", "elite"))
baseline$sports_static_component <- factor(baseline$sports_static_component, levels = 0:2, ordered = TRUE, labels = c("I", "II", "III"))
baseline$sports_static_component_ch <- factor(baseline$sports_static_component_ch, levels = 0:2, ordered = TRUE, labels = c("I", "II", "III"))
baseline$sports_dynamic_component <- factor(baseline$sports_dynamic_component, levels = 0:2, ordered = TRUE, labels = c("A", "B", "C"))
baseline$sports_dynamic_comp_ch <- factor(baseline$sports_dynamic_comp_ch, levels = 0:2, ordered = TRUE, labels = c("A", "B", "C"))
baseline$sports_classification <- factor(baseline$sports_classification, levels = 1:5, ordered = TRUE, labels = c("very low", "low", "moderate", "high", "very high"))
baseline$caa_malignancy <- factor(baseline$caa_malignancy, levels = c(0, 1), labels = c("benign", "malign"))
baseline$caa_ostia <- factor(baseline$caa_ostia, levels = c(0, 1), labels = c("one", "two"))
baseline$dgn_echo_diafun <- factor(baseline$dgn_echo_diafun, levels = c(1, 2, 3, 4, 5, 6), labels = c("normal", "grade 1", "grade 2", "grade 3", "undefined arrythmia", "undefined"))
baseline$dgn_ecg_classification <- factor(baseline$dgn_ecg_classification, levels = c(0, 1), labels = c("normal", "abnormal"))
baseline$ccta_quali <- factor(baseline$ccta_quali, levels = c(0, 1, 2, 3), ordered = TRUE, labels = c("1", "2", "3", "4"))
baseline$ccta_dominance <- factor(baseline$ccta_dominance, levels = c(0, 1, 2), labels = c("right", "left", "balanced"))
baseline$funct_ergo_findings <- factor(baseline$funct_ergo_findings, levels = c(1, 2, 3, 4), labels = c("clinical and electrical negative", "clinical negative, electrical positive", "clinical positive, electrical negative", "clinical and electrical positive"))
baseline$inv_dominance <- factor(baseline$inv_dominance, levels = c(0, 1, 2), labels = c("right", "left", "balanced"))
baseline$synp_symp <- factor(baseline$synp_symp, levels = c(0, 1, 2, 3), labels = c("incidental finding", "symptoms linked to caa", "symptoms linked to cad", "symptoms linked to other"))
baseline$hx_dm <- factor(baseline$hx_dm, levels = c(2, 1, 3, 0), labels = c("DM Typ1", "DM Typ2", "Prediabetes", "No DM"))
baseline$hx_hxtobacco <- factor(baseline$hx_hxtobacco, levels = c(1, 2), labels = c("Hx Tobacco", "No Hx Tobacco"))
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
          n_group2 = factor_group2[i, 2], per_group2 = factor_group2[i, 3],
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

table_cleaner <- function(df_all, df_stats) {
  new_df <- data.frame(Variable = character(), N = character(), mean_all = character(), mean_group1 = character(), mean_group2 = character(), 
  p_value_classic = character(), p_value_sim = character(), p_value_classic_corrected = character(), p_value_sim_corrected = character(), stringsAsFactors = FALSE)
  
  for (i in 1:nrow(df_stats)) {
    var <- df_stats[i, 1]
    n <- df_stats[i, 3]
    sum_all <- nrow(df_all)
    # sum_group1 <- max(df_stats[, 9], na.rm = TRUE)
    # sum_group2 <- max(df_stats[, 15], na.rm = TRUE)
    sum_group1 <- max(df_stats[, 9], na.rm = TRUE)
    sum_group2 <- max(df_stats[, 15], na.rm = TRUE)
    all <- c()
    group1 <- c()
    group2 <- c()
    p_value_classic <- c()
    p_value_sim <- c()
    if (df_stats[i, 2] == 1 && all(!is.na(df_stats[i, 21:23])) && all(df_stats[i, 21:23] > 0.05)){
      mean_all <- round(df_stats[i, 5], 1)
      sd_all <- round(df_stats[i, 6], 1)
      mean_group1 <- round(df_stats[i, 11], 1)
      sd_group1 <- round(df_stats[i, 12], 1)
      mean_group2 <- round(df_stats[i, 17], 1)
      sd_group2 <- round(df_stats[i, 18], 1)
      all <- c(all, paste(mean_all,"±", sd_all))
      group1 <- c(group1, paste(mean_group1, "±", sd_group1))
      group2 <- c(group2, paste(mean_group2, "±", sd_group2))
      p_value_classic <- c(p_value_classic, df_stats[i, 24])
      p_value_sim <- c(p_value_sim, df_stats[i, 28])
    }
    else if (df_stats[i, 2] == 1 && all(!is.na(df_stats[i, 21:23])) && any(df_stats[i, 21:23] < 0.05)){
      median_all <- round(df_stats[i, 7], 1)
      iqr_all <- round(df_stats[i, 8], 1)
      median_group1 <- round(df_stats[i, 13], 1)
      iqr_group1 <- round(df_stats[i, 14], 1)
      median_group2 <- round(df_stats[i, 19], 1)
      iqr_group2 <- round(df_stats[i, 20], 1)
      all <- c(all, paste(median_all,"(", iqr_all, ")"))
      group1 <- c(group1, paste(median_group1,"(", iqr_group1, ")"))
      group2 <- c(group2, paste(median_group2,"(", iqr_group2, ")"))
      p_value_classic <- c(p_value_classic, df_stats[i, 25])
      p_value_sim <- c(p_value_sim, df_stats[i, 28])
    }
    else {
      n_all <- df_stats[i, 3]
      per_all <- df_stats[i, 4] * 100
      n_group1 <- df_stats[i, 9]
      per_group1 <- (n_group1 / sum_group1) * 100
      n_group2 <- df_stats[i, 15]
      per_group2 <- (n_group2 / sum_group2) * 100
      all <- c(all, paste(n_all, "(", round(per_all), "%)"))
      group1 <- c(group1, paste(n_group1, "(", round(per_group1), "%)"))
      group2 <- c(group2, paste(n_group2, "(", round(per_group2), "%)"))
      p_value_classic <- c(p_value_classic, df_stats[i, 26])
      p_value_sim <- c(p_value_sim, df_stats[i, 28])
    }
    p_value_classic <- ifelse(is.na(p_value_classic), "NA", ifelse(p_value_classic < 0.05, round(p_value_classic, 3), round(p_value_classic, 2)))
    p_value_sim <- ifelse(is.na(p_value_sim), "NA", ifelse(p_value_sim < 0.05, round(p_value_sim, 3), round(p_value_sim, 2)))

    new_df <- rbind(new_df, c(Variable = var, N = paste(n, "(", round((n / sum_all) * 100), "%)"), mean_all = all[1],
    mean_group1 = group1[1], mean_group2 = group2[1], p_value_classic = p_value_classic, 
    p_value_sim = p_value_sim, p_value_classic_corrected = NA, p_value_sim_corrected = NA))
  }
  new_df <- new_df[, -c(8:9)]
  colnames(new_df) <- c("Variable", "N", "All", "Group 1", 
                        "Group 2", "P value classic", "P value simulated")
  return(new_df)
}

## ACTUAL ANALYSIS ########################################################################################################## 
vars <- c()
for (i in seq_along(yaml$demographics$dict_rename)) {
    vars <- c(vars, yaml$demographics$dict_rename[[i]])
}
data <- baseline %>% 
  select("record_id", all_of(vars), "ffr_0.8") %>%
  filter(!is.na(inv_ivusrest_mla))
# dataframe_malignancy <- fill_dataframe(data, "caa_malignancy", dataframe_calculations)
dataframe_ffr0.8 <- fill_dataframe(data, "ffr_0.8", dataframe_calculations)
# p_values_malignancy <- statistics_dataframe(data, dataframe_malignancy, "caa_malignancy")
p_values_ffr0.8 <- statistics_dataframe(data, dataframe_ffr0.8, "ffr_0.8")
# malignancy_df <- merge(dataframe_malignancy, p_values_malignancy, by = "Variable", all = TRUE)
ffr0.8_df <- merge(dataframe_ffr0.8, p_values_ffr0.8, by = "Variable", all = TRUE)

for (i in seq_along(yaml$demographics$dict_rename)) { 
    new_names <- names(yaml$demographics$dict_rename)[i] 
    old_names <- yaml$demographics$dict_rename[[i]] 
    ffr0.8_df$Variable <- gsub(old_names, new_names, ffr0.8_df$Variable)
}

for(i in seq_along(yaml$demographics$dict_rename)){
  print(names(yaml$demographics$dict_rename)[i])
  print(yaml$demographics$dict_rename[[i]])
}

# malignancy_df <- malignancy_df %>% select(-c("P_Value_Sim_Prop", "P_Value_Sim_Chi")) %>% rename(P_Value_Simulation = P_Value_Sim_Mean)
ffr0.8_df <- ffr0.8_df %>% select(-c("P_Value_Sim_Prop", "P_Value_Sim_Chi")) %>% rename(P_Value_Simulation = P_Value_Sim_Mean)
# ffr0.81_df <- ffr0.81_df %>% select(-c("P_Value_Sim_Prop", "P_Value_Sim_Chi")) %>% rename(P_Value_Simulation = P_Value_Sim_Mean)
# malignancy_df[, c(2:28)] <- sapply(malignancy_df[, c(2:28)], as.numeric)
ffr0.8_df[, c(2:28)] <- sapply(ffr0.8_df[, c(2:28)], as.numeric)
# ffr0.81_df[, c(2:28)] <- sapply(ffr0.81_df[, c(2:28)], as.numeric)

# malignancy_df_clean <- table_cleaner(baseline, malignancy_df)
# malignancy_df_clean <- malignancy_df_clean %>% mutate(p_value_classic_fdr = p.adjust(`P value classic`, method = "fdr"),
#                                 p_value_sim_fdr = p.adjust(`P value simulated`, method = "fdr"))
ffr0.8_df_clean <- table_cleaner(data, ffr0.8_df)
ffr0.8_df_clean <- ffr0.8_df_clean %>% mutate(p_value_classic_fdr = p.adjust(`P value classic`, method = "fdr"),
                                p_value_sim_fdr = p.adjust(`P value simulated`, method = "fdr"),
                                p_value_classic_holm = p.adjust(`P value classic`, method = "holm"),
                                p_value_sim_holm = p.adjust(`P value simulated`, method = "holm"))

write.csv(ffr0.8_df_clean, "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/ffr0.8_df.csv", row.names = FALSE)

saveRDS(ffr0.8_df_clean, file = paste0(yaml$demographics$output_dir_data, "/ffr0.8_df.rds"))

ffr_0.8 <- readRDS("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/ffr0.8_df.rds")
ffr_0.8 <- ffr_0.8 %>% mutate(p_value_classic_holm = p.adjust(`P value classic`, method = "holm"),
                                p_value_sim_holm = p.adjust(`P value simulated`, method = "holm"))
R.version.string
write.csv(ffr_0.8, "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/ffr0.8_df.csv", row.names = FALSE)