library(tidyverse)
library(ggplot2)
library(readxl)
library(yaml)

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
    )
  )

# baseline <- baseline %>% 
#   mutate(
#     sports_classification_ch = case_when(
#       sports_static_component_ch == 0 & sports_dynamic_component_ch == 0 ~ 1,
#       (sports_static_component_ch == 1 & sports_dynamic_component_ch == 0) | (sports_static_component_ch == 0 & sports_dynamic_component_ch ==1) ~ 2,
#       (sports_static_component_ch == 2 & sports_dynamic_component_ch == 0) | (sports_static_component_ch == 1 & sports_dynamic_component_ch == 1) | (sports_static_component_ch == 0 & sports_dynamic_component_ch == 2) ~ 3,
#       (sports_static_component_ch == 2 & sports_dynamic_component_ch == 1) | (sports_static_component_ch == 1 & sports_dynamic_component_ch ==2) ~ 4,
#       sports_static_component_ch == 2 & sports_dynamic_component_ch == 2 ~ 5,
#     )
#   )

## FACTORIZE ##########################################################################################################
binary_vars <- names(baseline)[
  sapply(names(baseline), function(var_name) {
    !grepl("sex_0_male_1_female", var_name) &&
    all(baseline[[var_name]] %in% c(0, 1, NA))
  })
]
for (var in binary_vars) {
  baseline[[var]] <- factor(baseline[[var]], levels = c(0, 1), labels = c("no", "yes"))
}
baseline$sex_0_male_1_female <- factor(baseline$sex_0_male_1_female, levels = c(0, 1), labels = c("male", "female"))
baseline$caa_modality <- factor(baseline$caa_modality, levels = c(0, 1, 2, 3, 4, 5), labels = c("ccta", "angiography", "echocardiography", "cmr", "surgery", "other"))
baseline$hx_sym_ap_ccs <- factor(baseline$hx_sym_ap_ccs, levels = 1:4, ordered = TRUE, labels = c("I", "II", "III", "IV"))
baseline$hx_sym_dysp_nyha <- factor(baseline$hx_sym_dysp_nyha, levels = 1:4, ordered = TRUE, labels = c("I", "II", "III", "IV"))
baseline$sports_level <- factor(baseline$sports_level, levels = 1:3, ordered = TRUE, labels = c("recreational", "competitive", "elite"))
baseline$sports_static_component <- factor(baseline$sports_static_component, levels = 0:2, ordered = TRUE, labels = c("I", "II", "III"))
baseline$sports_dynamic_component <- factor(baseline$sports_dynamic_component, levels = 0:2, ordered = TRUE, labels = c("A", "B", "C"))
baseline$sports_classification <- factor(baseline$sports_classification, levels = 1:5, ordered = TRUE, labels = c("very low", "low", "moderate", "high", "very high"))
baseline$caa_malignancy <- factor(baseline$caa_malignancy, levels = c(0, 1), labels = c("benign", "malign"))
############################################################################################################
# Two splits: one for benign and malign and the other for ffr above and below equal 0.8
# check if contains any NA --> shouldn't be the case, if yes check database
baseline %>% select(record_id, caa_malignancy) %>% filter(is.na(caa_malignancy))
baseline %>% select(record_id, inv_ffrdobu, caa_malignancy) %>% filter(is.na(inv_ffrdobu) & caa_malignancy == 1)

############################################################################################################
# simulation based approach
numerical_simulation_based <- function(data, numerical_variable, categorical_variable, reps = 1000) {
    data[[categorical_variable]] <- factor(data[[categorical_variable]], levels = c(0, 1), labels = c("success", "failure"))
    null_distn <- data %>%
        specify(response = !!sym(numerical_variable), explanatory = !!sym(categorical_variable)) %>%
        hypothesize(null = "independence") %>%
        generate(reps = reps, type = "permute") %>%
        calculate(stat = "diff in means", order = c("success", "failure"))
    obs_stat <- data %>%
        specify(response = !!sym(numerical_variable), explanatory = !!sym(categorical_variable)) %>%
        calculate(stat = "diff in means", order = c("success", "failure"))
    p_value <- get_p_value(null_distn, obs_stat, direction = "two sided")
    return(p_value)
}

baseline$caa_malignancy <- factor(baseline$caa_malignancy, levels = c(0, 1), labels = c("benign", "malign"))
test <- baseline %>% filter(!is.na(baseline$caa_malignancy))
null_distn <- test %>%
    specify(age ~ caa_malignancy) %>%
    hypothesize(null = "independence") %>%
    generate(reps = 5000, type = "permute") %>%
    calculate(stat = "diff in means", order = c("malign", "benign"))
obs_stat <- test %>%
    specify(age ~ caa_malignancy) %>%
    calculate(stat = "diff in means", order = c("malign", "benign"))
p_value <- get_p_value(null_distn, obs_stat, direction = "two sided")
p_value

visualize(null_distn) +
    geom_vline(aes(xintercept = stat), 
    data = obs_stat, 
    color = "red")
