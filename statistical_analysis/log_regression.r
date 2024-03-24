library(tidyverse)
library(broom)
library(ggplot2)
library(yardstick)
library(pROC)
library(randomForest)
library(car)
library(readxl)

baseline <- readRDS("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/baseline.rds")
invasive <- read_excel("C:/WorkingData/Documents/3_Research/IVUS_data.xlsx")

invasive <- invasive %>% rename(
    patient_id = `NARCO ID`,
    inv_ivusrest_mla = MLA_rest,
    inv_ivusrest_mla_c = MLA_circ_rest,
    inv_ivusrest_mla_h = MLA_h_rest,
    inv_ivusrest_mla_w = MLA_w_rest,
    inv_ivusrest_mla_ellip = MLA_elllip_rest,
    inv_ivusrest_ostial_a = ostial_a_rest,
    inv_ivusrest_ostial_c = ostial_circ_rest,
    inv_ivusrest_ostial_h = ostial_h_rest,
    inv_ivusrest_ostial_w = ostial_w_rest,
    inv_ivusrest_ostial_ellip = ostial_ellip_rest,
    inv_ivusrest_dist_a = reference_a_rest,
    inv_ivusrest_dist_c = reference_circ_rest,
    inv_ivusdobu_mla = MLA_stress,
    inv_ivusdobu_mla_c = MLA_circ_stress,
    inv_ivusdobu_mla_h = MLA_h_stress,
    inv_ivusdobu_mla_w = MLA_w_stress,
    inv_ivusdobu_mla_ellip = MLA_elllip_stress,
    inv_ivusdobu_ostial_a = ostial_a_stress,
    inv_ivusdobu_ostial_c = ostial_circ_stress,
    inv_ivusdobu_ostial_h = ostial_h_stress,
    inv_ivusdobu_ostial_w = ostial_w_stress,
    inv_ivusdobu_ostial_ellip = ostial_ellip_stress,
    inv_ivusdobu_dist_a = reference_a_stress,
    inv_ivusdobu_dist_c = reference_circ_stress 
) %>% mutate(inv_ivusdobu_mla = as.numeric(inv_ivusdobu_mla))
invasive <- invasive %>% mutate(
    patient_id = paste0("NARCO_", patient_id)
)

# drop last 3 columns
invasive <- invasive %>% select(-c(26:28))

# merge with baseline by record_id
baseline <- left_join(as_tibble(baseline), as_tibble(invasive), by = "patient_id")
baseline <- baseline %>% mutate(
    inv_ivusrest_mla = ifelse(is.na(inv_ivusrest_mla.y), inv_ivusrest_mla.x, inv_ivusrest_mla.y),
    inv_ivusrest_mla_ellip = ifelse(is.na(inv_ivusrest_mla_ellip.y), inv_ivusrest_mla_ellip.x, inv_ivusrest_mla_ellip.y),
    inv_ivusrest_dist_a = ifelse(is.na(inv_ivusrest_dist_a.y), inv_ivusrest_dist_a.x, inv_ivusrest_dist_a.y),
    inv_ivusdobu_mla = ifelse(is.na(inv_ivusdobu_mla.y), inv_ivusdobu_mla.x, inv_ivusdobu_mla.y),
    inv_ivusdobu_mla_ellip = ifelse(is.na(inv_ivusdobu_mla_ellip.y), inv_ivusdobu_mla_ellip.x, inv_ivusdobu_mla_ellip.y),
    inv_ivusrest_dist_a = ifelse(is.na(inv_ivusrest_dist_a.y), inv_ivusrest_dist_a.x, inv_ivusrest_dist_a.y),
    inv_ivusdobu_dist_a = ifelse(is.na(inv_ivusdobu_dist_a.y), inv_ivusdobu_dist_a.x, inv_ivusdobu_dist_a.y)
) %>% 
select(-c("inv_ivusrest_mla.y", "inv_ivusrest_mla.x", "inv_ivusrest_dist_a.y", "inv_ivusrest_dist_a.x", 
          "inv_ivusdobu_mla.y", "inv_ivusdobu_mla.x", "inv_ivusdobu_mla_ellip.y", "inv_ivusdobu_mla_ellip.x", 
          "inv_ivusdobu_dist_a.y", "inv_ivusdobu_dist_a.x"))

# quick data prepping
baseline <- baseline %>% 
    mutate(
        ccta_ostial_pn = ccta_ostial_a / ccta_dist_a,
        ccta_quali = as.numeric(ccta_quali),
        ccta_quali = ifelse(ccta_photon == "yes", ccta_quali + 1, ifelse(is.na(ccta_quali), 2, ccta_quali)),
        ffr_0.8 = ifelse(ffr_0.8 == "yes", 1, 0),
        ffr_0.81 = ifelse(ffr_0.81 == "yes", 1, 0),
        ccta_ostial_a_bsa = ccta_ostial_a / bsa,
        ccta_mla_a_bsa = ccta_mla_a / bsa,
        ccta_ostial_h_bsa = ccta_ostial_h / bsa,
        ccta_mla_h_bsa = ccta_mla_h / bsa,
        ccta_ostial_w_bsa = ccta_ostial_w / bsa,
        ccta_mla_w_bsa = ccta_mla_w / bsa,
        ccta_imc_length_bsa = ccta_imc_length / bsa,
        ccta_aa_degree_bsa = ccta_aa_degree / bsa,
        inv_ivus_dist = pmax(inv_ivusrest_dist_a, inv_ivusdobu_dist_a),
        inv_ivus_dist_c = pmax(inv_ivusrest_dist_c, inv_ivusdobu_dist_c),
        inv_ivusrest_pn = inv_ivusrest_mla / inv_ivus_dist,
        inv_ivusdobu_pn = inv_ivusdobu_mla / inv_ivus_dist,
        inv_ivusrest_mla_bsa = inv_ivusrest_mla / bsa,
        inv_ivusdobu_mla_bsa = inv_ivusdobu_mla / bsa,
        inv_ivusrest_mla_h_bsa = inv_ivusrest_mla_h / bsa,
        inv_ivusdobu_mla_h_bsa = inv_ivusdobu_mla_h / bsa,
        inv_ivusrest_mla_w_bsa = inv_ivusrest_mla_w / bsa,
        inv_ivusdobu_mla_w_bsa = inv_ivusdobu_mla_w / bsa,
        inv_ivusrest_hypo = inv_ivusrest_mla_c / inv_ivus_dist_c,
        inv_ivusdobu_hypo = inv_ivusdobu_mla_c / inv_ivus_dist_c,
        inv_ivus_delta = inv_ivusdobu_mla - inv_ivusrest_mla,
        inv_ivusrest_ostial_a_bsa = inv_ivusrest_ostial_a / bsa,
        inv_ivusdobu_ostial_a_bsa = inv_ivusdobu_ostial_a / bsa,
        inv_ivusrest_ostial_h_bsa = inv_ivusrest_ostial_h / bsa,
        inv_ivusdobu_ostial_h_bsa = inv_ivusdobu_ostial_h / bsa,
        inv_ivusrest_ostial_w_bsa = inv_ivusrest_ostial_w / bsa,
        inv_ivusdobu_ostial_w_bsa = inv_ivusdobu_ostial_w / bsa,
        inv_ivusrest_ostial_pn = inv_ivusrest_ostial_a / inv_ivus_dist,
        inv_ivusdobu_ostial_pn = inv_ivusdobu_ostial_a / inv_ivus_dist,
    )
saveRDS(baseline, "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/baseline.rds")

# functions needed
simple_logistic_regression <- function(baseline_data, response, explanatory, family = "binomial") {
    data <- baseline_data %>% 
        select(!!sym(response), !!sym(explanatory)) %>% 
        drop_na()
    formula <- as.formula(paste(response, "~", explanatory))
    mdl <- glm(formula, data = data, family = family)

    max <- max(data[[explanatory]])
    min <- min(data[[explanatory]])
    buffer <- (max - min) / 50

    explanatory_data <- tibble(
        !!sym(explanatory) := seq(min - 2 * buffer, max + 2 * buffer, length.out = 50)
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
        scale_x_continuous(breaks = seq(round(min - buffer, 0), round(max + buffer), round(max / 10, 1)))
    
    plot_odds <- ggplot(prediction_data, aes(x = !!sym(explanatory), y = odds_ratio)) + 
        geom_point() + 
        geom_line() +
        geom_hline(yintercept = 1, linetype = "dashed") +
        scale_x_continuous(breaks = seq(round(min - buffer, 0), round(max + buffer), round(max / 10, 1))) +
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

logistic_regression_categorical <- function(baseline_data, response, explanatory, family = "binomial") {
    data <- baseline_data %>% 
        select(!!sym(response), !!sym(explanatory)) %>% 
        drop_na()

    formula <- as.formula(paste(response, "~", explanatory))
    mdl <- glm(formula, data = data, family = family)
    odds_ratio <- exp(coef(mdl))[2]
    lower_ci <- exp(confint(mdl))[2]
    upper_ci <- exp(confint(mdl))[4]
    print(formula)
    tryCatch(
        {
            outcome <- table(baseline_data[[explanatory]], baseline_data[[response]])
            confusion <- conf_mat(outcome)
            plot_acc <- autoplot(confusion)
        },
        error = function(e) {
            print("Unfactor the explanatory variable first for confusion matrix and accuracy plot!")
        }
    )
    return(list(mdl, odds_ratio, lower_ci, upper_ci, confusion, plot_acc))
}

# all models
# ostium
all_ost_a_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_a_bsa")
mdl_ost_a_bsa_ffr <- all_ost_a_bsa[[1]] # significant
all_ost_ellip <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_elliptic")
mdl_ost_ellip_ffr <- all_ost_ellip[[1]] # significant
all_ost_h_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_h_bsa")
mdl_ost_h_bsa_ffr <- all_ost_h_bsa[[1]] # not significant!!!!
all_ost_w_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_w_bsa")
mdl_ost_w_bsa_ffr <- all_ost_w_bsa[[1]] # significant
all_ost_pn <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_pn")
mdl_ost_pn_ffr <- all_ost_pn[[1]] # significant

# mla
all_mla_a_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_a_bsa")
mdl_mla_a_bsa_ffr <- all_mla_a_bsa[[1]] # significant
all_mla_ellip <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_elliptic")
mdl_mla_ellip_ffr <- all_mla_ellip[[1]] # not significant!!!
all_mla_h_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_h_bsa")
mdl_mla_h_bsa_ffr <- all_mla_h_bsa[[1]] # not significant!!!
all_mla_w_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_w_bsa")
mdl_mla_w_bsa_ffr <- all_mla_w_bsa[[1]] # significant
all_mla_pn <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_pn_dist")
mdl_mla_pn_ffr <- all_mla_pn[[1]] # significant
all_ostial_pn <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_pn")
mdl_ostial_pn_ffr <- all_ostial_pn[[1]] # significant

# rest of the variables
all_imc_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_imc_length_bsa")
mdl_imc_bsa_ffr <- all_imc_bsa[[1]] # not significant!!!
all_aa_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_aa_degree_bsa")
mdl_aa_bsa_ffr <- all_aa_bsa[[1]] # not significant!!!

# invasive variables
mdl_ivusdobu_ffr <- all_ivusdobu[[1]] # significant
all_ivusrest_pn <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_pn")
mdl_ivusrest_pn_ffr <- all_ivusrest_pn[[1]] # significant
all_ivusrest_ostial_pn <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_ostial_pn")
mdl_ivusrest_ostial_pn_ffr <- all_ivusrest_ostial_pn[[1]] # significant
all_ivusdobu_ostial_pn <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusdobu_ostial_pn")
mdl_ivusdobu_ostial_pn_ffr <- all_ivusdobu_ostial_pn[[1]] # significant
all_ivusdobu_pn <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusdobu_pn")
mdl_ivusdobu_pn_ffr <- all_ivusdobu_pn[[1]] # significant
all_ivusrest_ellip <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_ellip")
mdl_ivusrest_ellip_ffr <- all_ivusrest_ellip[[1]] # significant
all_ivusdobu_ellip <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusdobu_mla_ellip")
mdl_ivusdobu_ellip_ffr <- all_ivusdobu_ellip[[1]] # significant 

all_ivusrest_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_bsa")
mdl_ivusrest_bsa_ffr <- all_ivusrest_bsa[[1]] # significant
all_ivusdobu_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusdobu_mla_bsa")
mdl_ivusdobu_bsa_ffr <- all_ivusdobu_bsa[[1]] # significant
all_ivusrest_h_bsa_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_h_bsa")
mdl_ivusrest_h_bsa_ffr <- all_ivusrest_h_bsa_ffr[[1]] # not significant
all_ivusdobu_h_bsa_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusdobu_mla_h_bsa")
mdl_ivusdobu_h_bsa_ffr <- all_ivusdobu_h_bsa_ffr[[1]] # not significant
all_ivusrest_w_bsa_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_mla_w_bsa")
mdl_ivusrest_w_bsa_ffr <- all_ivusrest_w_bsa_ffr[[1]] # significant
all_ivusdobu_w_bsa_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusdobu_mla_w_bsa")
mdl_ivusdobu_w_bsa_ffr <- all_ivusdobu_w_bsa_ffr[[1]] # not significant
all_ivusrest_hypo <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_hypo")
mdl_ivusrest_hypo_ffr <- all_ivusrest_hypo[[1]] # not significant
all_ivusdobu_hypo <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusdobu_hypo")
mdl_ivusdobu_hypo_ffr <- all_ivusdobu_hypo[[1]] # not significant

# invasive ostial measurements
all_ivusrest_ostial_a <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_ostial_a_bsa")
mdl_ivusrest_ostial_a_ffr <- all_ivusrest_ostial_a[[1]] # significant
all_ivusdobu_ostial_a_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusdobu_ostial_a_bsa")
mdl_ivusdobu_ostial_a_bsa_ffr <- all_ivusdobu_ostial_a_bsa[[1]] # significant
all_ivusrest_ostial_h_bsa_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_ostial_h_bsa")
mdl_ivusrest_ostial_h_bsa_ffr <- all_ivusrest_ostial_h_bsa_ffr[[1]] # not significant
all_ivusdobu_ostial_h_bsa_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusdobu_ostial_h_bsa")
mdl_ivusdobu_ostial_h_bsa_ffr <- all_ivusdobu_ostial_h_bsa_ffr[[1]] # not significant
all_ivusrest_ostial_w_bsa_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_ostial_w_bsa")
mdl_ivusrest_ostial_w_bsa_ffr <- all_ivusrest_ostial_w_bsa_ffr[[1]] # significant
all_ivusdobu_ostial_w_bsa_ffr <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusdobu_ostial_w_bsa")
mdl_ivusdobu_ostial_w_bsa_ffr <- all_ivusdobu_ostial_w_bsa_ffr[[1]] # significant
all_ivusrest_ostial_ellip <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusrest_ostial_ellip")
mdl_ivusrest_ostial_ellip_ffr <- all_ivusrest_ostial_ellip[[1]] # significant
all_ivusdobu_ostial_ellip <- simple_logistic_regression(baseline, "ffr_0.8", "inv_ivusdobu_ostial_ellip")
mdl_ivusdobu_ostial_ellip_ffr <- all_ivusdobu_ostial_ellip[[1]] # significant


# forest plot
# bsa adjusted
forest_df_bsa <- data.frame(
    var = c("Ostial area (BSA) [mm2]", "Ostial elliptic ratio", "Ostial max. diameter (BSA) [mm]", "Ostial min. diameter (BSA) [mm]", 
            "Ostial proximal narrowing [%]", "MLA (BSA) [mm2]", "ML elliptic ratio", "ML max. diameter (BSA) [mm]", "ML min. diameter (BSA) [mm]", 
            "ML proximal narrowing [%]", "Intramural length (BSA) [mm]", "Take-off angle (BSA) [°]"),
    odds_ratio = c(exp(mdl_ost_a_bsa_ffr$coefficients[2]), exp(mdl_ost_ellip_ffr$coefficients[2]), 
            exp(mdl_ost_h_bsa_ffr$coefficients[2]), exp(mdl_ost_w_bsa_ffr$coefficients[2]), 
            exp(mdl_ost_pn_ffr$coefficients[2]), exp(mdl_mla_a_bsa_ffr$coefficients[2]), 
            exp(mdl_mla_ellip_ffr$coefficients[2]), exp(mdl_mla_h_bsa_ffr$coefficients[2]), 
            exp(mdl_mla_w_bsa_ffr$coefficients[2]), exp(mdl_mla_pn_ffr$coefficients[2]),
            exp(mdl_imc_bsa_ffr$coefficients[2]), exp(mdl_aa_bsa_ffr$coefficients[2])),
    lower_ci = c(exp(confint(mdl_ost_a_bsa_ffr))[2], exp(confint(mdl_ost_ellip_ffr))[2], 
            exp(confint(mdl_ost_h_bsa_ffr))[2], exp(confint(mdl_ost_w_bsa_ffr))[2], 
            exp(confint(mdl_ost_pn_ffr))[2], exp(confint(mdl_mla_a_bsa_ffr))[2], 
            exp(confint(mdl_mla_ellip_ffr))[2], exp(confint(mdl_mla_h_bsa_ffr))[2], 
            exp(confint(mdl_mla_w_bsa_ffr))[2], exp(confint(mdl_mla_pn_ffr))[2], 
            exp(confint(mdl_imc_bsa_ffr))[2], exp(confint(mdl_aa_bsa_ffr))[2]),
    upper_ci = c(exp(confint(mdl_ost_a_bsa_ffr))[4], exp(confint(mdl_ost_ellip_ffr))[4], 
            exp(confint(mdl_ost_h_bsa_ffr))[4], exp(confint(mdl_ost_w_bsa_ffr))[4], 
            exp(confint(mdl_ost_pn_ffr))[4], exp(confint(mdl_mla_a_bsa_ffr))[4], 
            exp(confint(mdl_mla_ellip_ffr))[4], exp(confint(mdl_mla_h_bsa_ffr))[4], 
            exp(confint(mdl_mla_w_bsa_ffr))[4], exp(confint(mdl_mla_pn_ffr))[4], 
            exp(confint(mdl_imc_bsa_ffr))[4], exp(confint(mdl_aa_bsa_ffr))[4])
)
forest_continuous_ffr_bsa <- ggplot(data = forest_df_bsa, aes(y = var, x = odds_ratio, xmin = lower_ci, xmax = upper_ci)) +
    geom_point() +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, 6, 1)) +
    theme_classic()

forest_cut <- forest_df_bsa %>%
    mutate(
        upper_ci = ifelse(upper_ci > 6, 6, upper_ci)
    )

forest_cut_ffr_bsa <- ggplot(data = forest_cut, aes(y = var, x = odds_ratio, xmin = lower_ci, xmax = upper_ci)) +
    geom_point() +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, 7, 1)) +
    geom_vline(xintercept = 6, linetype = "solid", color = "white", lwd = 2) +
    geom_segment(aes(x = 6.1, xend = 7, y = 2, yend = 2), linetype = "dashed") +
    geom_segment(aes(x = 6.1, xend = 7, y = 8, yend = 8), linetype = "dashed") +
    theme_classic()

ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/forest_ccta_ffr.png", forest_cut_ffr_bsa, width = 5, height = 4)


# forest invasive bsa adjusted
forest_df_invasive <- data.frame(
    var = c("ivusrest_pn", "ivusdobu_pn", "ivusrest_mla_ellip", "ivusdobu_mla_ellip", 
            "ivusrest_mla_bsa", "ivusdobu_mla_bsa",
            "ivusrest_mla_h_bsa", "ivusdobu_mla_h_bsa", "ivusrest_mla_w_bsa", "ivusdobu_mla_w_bsa", 
            "ivusrest_ostial_a", "ivusdobu_ostial_a_bsa", 
            "ivusrest_ostial_h_bsa", "ivusdobu_ostial_h_bsa", "ivusrest_ostial_w_bsa", 
            "ivusdobu_ostial_w_bsa", "ivusrest_ostial_ellip", "ivusdobu_ostial_ellip"),
    odds_ratio = c(
        exp(coef(mdl_ivusrest_pn_ffr)[2]), exp(coef(mdl_ivusdobu_pn_ffr)[2]),
        exp(coef(mdl_ivusrest_ellip_ffr)[2]), exp(coef(mdl_ivusdobu_ellip_ffr)[2]),
        exp(coef(mdl_ivusrest_bsa_ffr)[2]), exp(coef(mdl_ivusdobu_bsa_ffr)[2]),
        exp(coef(mdl_ivusrest_h_bsa_ffr)[2]), exp(coef(mdl_ivusdobu_h_bsa_ffr)[2]),
        exp(coef(mdl_ivusrest_w_bsa_ffr)[2]), exp(coef(mdl_ivusdobu_w_bsa_ffr)[2]),
        exp(coef(mdl_ivusrest_ostial_a_ffr)[2]), exp(coef(mdl_ivusdobu_ostial_a_bsa_ffr)[2]),
        exp(coef(mdl_ivusrest_ostial_h_bsa_ffr)[2]), exp(coef(mdl_ivusdobu_ostial_h_bsa_ffr)[2]),
        exp(coef(mdl_ivusrest_ostial_w_bsa_ffr)[2]), exp(coef(mdl_ivusdobu_ostial_w_bsa_ffr)[2]),
        exp(coef(mdl_ivusrest_ostial_ellip_ffr)[2]), exp(coef(mdl_ivusdobu_ostial_ellip_ffr)[2])
    ),
    lower_ci = c(
        exp(confint(mdl_ivusrest_pn_ffr)[2, 1]), exp(confint(mdl_ivusdobu_pn_ffr)[2, 1]),
        exp(confint(mdl_ivusrest_ellip_ffr)[2, 1]), exp(confint(mdl_ivusdobu_ellip_ffr)[2, 1]),
        exp(confint(mdl_ivusrest_bsa_ffr)[2, 1]), exp(confint(mdl_ivusdobu_bsa_ffr)[2, 1]),
        exp(confint(mdl_ivusrest_h_bsa_ffr)[2, 1]), exp(confint(mdl_ivusdobu_h_bsa_ffr)[2, 1]),
        exp(confint(mdl_ivusrest_w_bsa_ffr)[2, 1]), exp(confint(mdl_ivusdobu_w_bsa_ffr)[2, 1]),
        exp(confint(mdl_ivusrest_ostial_a_ffr)[2, 1]), exp(confint(mdl_ivusdobu_ostial_a_bsa_ffr)[2, 1]),
        exp(confint(mdl_ivusrest_ostial_h_bsa_ffr)[2, 1]), exp(confint(mdl_ivusdobu_ostial_h_bsa_ffr)[2, 1]),
        exp(confint(mdl_ivusrest_ostial_w_bsa_ffr)[2, 1]), exp(confint(mdl_ivusdobu_ostial_w_bsa_ffr)[2, 1]),
        exp(confint(mdl_ivusrest_ostial_ellip_ffr)[2, 1]), exp(confint(mdl_ivusdobu_ostial_ellip_ffr)[2, 1])
    ),
    upper_ci = c(
        exp(confint(mdl_ivusrest_pn_ffr)[2, 2]), exp(confint(mdl_ivusdobu_pn_ffr)[2, 2]),
        exp(confint(mdl_ivusrest_ellip_ffr)[2, 2]), exp(confint(mdl_ivusdobu_ellip_ffr)[2, 2]),
        exp(confint(mdl_ivusrest_bsa_ffr)[2, 2]), exp(confint(mdl_ivusdobu_bsa_ffr)[2, 2]),
        exp(confint(mdl_ivusrest_h_bsa_ffr)[2, 2]), exp(confint(mdl_ivusdobu_h_bsa_ffr)[2, 2]),
        exp(confint(mdl_ivusrest_w_bsa_ffr)[2, 2]), exp(confint(mdl_ivusdobu_w_bsa_ffr)[2, 2]),
        exp(confint(mdl_ivusrest_ostial_a_ffr)[2, 2]), exp(confint(mdl_ivusdobu_ostial_a_bsa_ffr)[2, 2]),
        exp(confint(mdl_ivusrest_ostial_h_bsa_ffr)[2, 2]), exp(confint(mdl_ivusdobu_ostial_h_bsa_ffr)[2, 2]),
        exp(confint(mdl_ivusrest_ostial_w_bsa_ffr)[2, 2]), exp(confint(mdl_ivusdobu_ostial_w_bsa_ffr)[2, 2]),
        exp(confint(mdl_ivusrest_ostial_ellip_ffr)[2, 2]), exp(confint(mdl_ivusdobu_ostial_ellip_ffr)[2, 2])
    )
)
forest_invasive_ffr <- ggplot(data = forest_df_invasive, aes(y = var, x = odds_ratio, xmin = lower_ci, xmax = upper_ci)) +
    geom_point() +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 9.5, ymax = Inf), fill = "#193bac", alpha = 0.05) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 9.5), fill = "#ebc90b", alpha = 0.05) +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    geom_hline(yintercept = 9.5, linetype = "dotted") +
    scale_x_continuous(breaks = seq(0, 7, 1)) +
    theme_classic()

ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/forest_invasive_ffr.png", forest_invasive_ffr, width = 5, height = 4)
# regression models #################################################################################################
# exploration of continuous variables that were significant
## ostium
### area
summary(mdl_ost_a_bsa_ffr) # significant
prediction_data_ost_a <- all_ost_a_bsa[[2]]
print(prediction_data_ost_a, n = 50)
plot_log_ost_a <- all_ost_a_bsa[[3]]
plot_log_ost_a
plot_odds_ost_a <- all_ost_a_bsa[[4]]
plot_odds_ost_a
confusion_ost_a <- all_ost_a_bsa[[5]]
confusion_ost_a
summary(confusion_ost_a, event_level = "second")
plot_acc_ost_a <- all_ost_a_bsa[[6]]
plot_acc_ost_a

### elliptic
summary(mdl_ost_ellip_ffr) # significant
prediction_data_ost_ellip <- all_ost_ellip[[2]]
print(prediction_data_ost_ellip, n = 50)
plot_log_ost_ellip <- all_ost_ellip[[3]]
plot_log_ost_ellip
plot_odds_ost_ellip <- all_ost_ellip[[4]]
plot_odds_ost_ellip

### width
summary(mdl_ost_w_bsa_ffr) # significant
prediction_data_ost_w <- all_ost_w_bsa[[2]]
print(prediction_data_ost_w, n = 50)
plot_log_ost_w <- all_ost_w_bsa[[3]]
plot_log_ost_w
plot_odds_ost_w <- all_ost_w_bsa[[4]]
plot_odds_ost_w

### proximal narrowing
summary(mdl_ost_pn_bsa_ffr) # not significant
prediction_data_ost_pn <- all_ost_pn_bsa[[2]]
print(prediction_data_ost_pn, n = 50)
plot_log_ost_pn <- all_ost_pn_bsa[[3]]
plot_log_ost_pn
plot_odds_ost_pn <- all_ost_pn_bsa[[4]]
plot_odds_ost_pn

## mla
### area
summary(mdl_mla_a_bsa_ffr) # significant
prediction_data_mla <- all_mla_a_bsa[[2]]
print(prediction_data_mla, n = 50)
plot_log_mla <- all_mla_a_bsa[[3]]
plot_log_mla
plot_odds_mla <- all_mla_a_bsa[[4]]
plot_odds_mla
confusion_mla <- all_mla_a_bsa[[5]]
confusion_mla
summary(confusion_mla, event_level = "second")
plot_acc_mla <- all_mla_a_bsa[[6]]
plot_acc_mla

### elliptic
summary(mdl_mla_ellip_bsa_ffr) # significant
prediction_data_mla_ellip <- all_mla_ellip_bsa[[2]]
print(prediction_data_mla_ellip, n = 50)
plot_log_mla_ellip <- all_mla_ellip_bsa[[3]]
plot_log_mla_ellip
plot_odds_mla_ellip <- all_mla_ellip_bsa[[4]]
plot_odds_mla_ellip

### width
summary(mdl_mla_w_bsa_ffr) # significant
prediction_data_mla_w <- all_mla_w_bsa[[2]]
print(prediction_data_mla, n = 50)
plot_log_mla_w <- all_mla_w_bsa[[3]]
plot_log_mla_w
plot_odds_mla_w <- all_mla_w_bsa[[4]]
plot_odds_mla_w

### proximal narrowing
summary(mdl_mla_pn_bsa_ffr) # not significant
prediction_data_mla_pn <- all_mla_pn[[2]]
print(prediction_data_mla, n = 50)
plot_log_mla_pn <- all_mla_pn_bsa[[3]]
plot_log_mla_pn
plot_odds_mla_pn <- all_mla_pn_bsa[[4]]
plot_odds_mla_pn

## invasive
### ivus rest
summary(mdl_ivusrest_bsa_ffr) # significant
prediction_data_ivusrest <- all_ivusrest_bsa[[2]]
print(prediction_data_ivusrest, n = 50)
plot_log_ivusrest <- all_ivusrest_bsa[[3]]
plot_log_ivusrest
plot_odds_ivusrest <- all_ivusrest_bsa[[4]]
plot_odds_ivusrest
confusion_ivusrest <- all_ivusrest_bsa[[5]]
confusion_ivusrest
summary(confusion_ivusrest, event_level = "second")
plot_acc_ivusrest <- all_ivusrest_bsa[[6]]
plot_acc_ivusrest

### ivus dobu
summary(mdl_ivusdobu_bsa_ffr) # significant
prediction_data_ivusdobu <- all_ivusdobu_bsa[[2]]
print(prediction_data_ivusdobu, n = 50)
plot_log_ivusdobu <- all_ivusdobu_bsa[[3]]
plot_log_ivusdobu
plot_odds_ivusdobu <- all_ivusdobu_bsa[[4]]
plot_odds_ivusdobu
confusion_ivusdobu <- all_ivusdobu_bsa[[5]]
confusion_ivusdobu
summary(confusion_ivusdobu, event_level = "second")
plot_acc_ivusdobu <- all_ivusdobu_bsa[[6]]
plot_acc_ivusdobu

### ivus dobu pn
summary(mdl_ivusdobu_pn_bsa_ffr) # significant
prediction_data_ivusdobu_pn <- all_ivusdobu_pn_bsa[[2]]
print(prediction_data_ivusdobu_pn, n = 50)
plot_log_ivusdobu_pn <- all_ivusdobu_pn_bsa[[3]]
plot_log_ivusdobu_pn
plot_odds_ivusdobu_pn <- all_ivusdobu_pn_bsa[[4]]
plot_odds_ivusdobu_pn

### ivus rest ellip
summary(mdl_ivusrest_ellip_bsa_ffr) # significant
prediction_data_ivusrest_ellip <- all_ivusrest_ellip_bsa[[2]]
print(prediction_data_ivusrest_ellip, n = 50)
plot_log_ivusrest_ellip <- all_ivusrest_ellip_bsa[[3]]
plot_log_ivusrest_ellip
plot_odds_ivusrest_ellip <- all_ivusrest_ellip_bsa[[4]]
plot_odds_ivusrest_ellip

### ivus dobu ellip
summary(mdl_ivusdobu_ellip_bsa_ffr) # significant
prediction_data_ivusdobu_ellip <- all_ivusdobu_ellip_bsa[[2]]
print(prediction_data_ivusdobu_ellip, n = 50)
plot_log_ivusdobu_ellip <- all_ivusdobu_ellip_bsa[[3]]
plot_log_ivusdobu_ellip
plot_odds_ivusdobu_ellip <- all_ivusdobu_ellip_bsa[[4]]
plot_odds_ivusdobu_ellip

# roc analysis ######################################################################################################
## ostium
### area
data_ostium <- baseline %>% select(ffr_0.8, ccta_ostial_a_bsa) %>% drop_na()
roc_ost_a <- roc(data_ostium$ffr_0.8, mdl_ost_a_bsa_ffr$fitted.values)
cutoff_best <- coords(roc_ost_a, x = "best", best.method = "youden", ret = "all")
cutoff_trueneg <- coords(roc_ost_a, x = "best", best.method = "youden", best.weights = c(1, 0.24), ret = "all") # with weights for high true negative
cutoff_truepos <- coords(roc_ost_a, x = "best", best.method = "youden", best.weights = c(1, 0.76), ret = "all") # with weights for high true positive
plot_log_ost_a +
    geom_hline(yintercept = cutoff_trueneg$threshold, linetype = "dashed") + 
    geom_hline(yintercept = cutoff_best$threshold, linetype = "solid") + 
    geom_hline(yintercept = cutoff_truepos$threshold, linetype = "dashed")

### elliptic
data_ostium <- baseline %>% select(ffr_0.8, ccta_ostial_elliptic) %>% drop_na()
roc_ost_ellip <- roc(data_ostium$ffr_0.8, mdl_ost_ellip_ffr$fitted.values)
cutoff_best <- coords(roc_ost_ellip, x = "best", best.method = "youden", ret = "all")
cutoff_trueneg <- coords(roc_ost_ellip, x = "best", best.method = "youden", best.weights = c(1, 1/3), ret = "all") # with weights for high true negative
cutoff_truepos <- coords(roc_ost_ellip, x = "best", best.method = "youden", best.weights = c(1, 2/3), ret = "all") # with weights for high true positive
plot_log_ost_ellip +
    geom_hline(yintercept = cutoff_trueneg$threshold, linetype = "dashed") + 
    geom_hline(yintercept = mean(cutoff_best$threshold), linetype = "solid") + 
    geom_hline(yintercept = cutoff_truepos$threshold, linetype = "dashed")

### width
data_ostium <- baseline %>% select(ffr_0.8, ccta_ostial_w_bsa) %>% drop_na()
roc_ost_w <- roc(data_ostium$ffr_0.8, mdl_ost_w_bsa_ffr$fitted.values)
cutoff_best <- coords(roc_ost_w, x = "best", best.method = "youden", ret = "all")
cutoff_trueneg <- coords(roc_ost_w, x = "best", best.method = "youden", best.weights = c(1, 1/3), ret = "all") # with weights for high true negative
cutoff_truepos <- coords(roc_ost_w, x = "best", best.method = "youden", best.weights = c(1, 2/3), ret = "all") # with weights for high true positive
plot_log_ost_w +
    geom_hline(yintercept = cutoff_trueneg$threshold, linetype = "dashed") + 
    geom_hline(yintercept = cutoff_best$threshold, linetype = "solid") + 
    geom_hline(yintercept = cutoff_truepos$threshold, linetype = "dashed")

### proximal narrowing
data_ostium <- baseline %>% select(ffr_0.8, ccta_ostial_pn) %>% drop_na()
roc_ost_pn <- roc(data_ostium$ffr_0.8, mdl_ost_pn_ffr$fitted.values)
cutoff_best <- coords(roc_ost_pn, x = "best", best.method = "youden", ret = "all")
cutoff_trueneg <- coords(roc_ost_pn, x = "best", best.method = "youden", best.weights = c(1, 1/3), ret = "all") # with weights for high true negative
cutoff_truepos <- coords(roc_ost_pn, x = "best", best.method = "youden", best.weights = c(1, 2/3), ret = "all") # with weights for high true positive
plot_log_ost_pn +
    geom_hline(yintercept = cutoff_trueneg$threshold, linetype = "dashed") + 
    geom_hline(yintercept = cutoff_best$threshold, linetype = "solid") + 
    geom_hline(yintercept = cutoff_truepos$threshold, linetype = "dashed")

## mla
data_mla <- baseline %>% select(ffr_0.8, ccta_mla_a_bsa) %>% drop_na()
### area
roc_mla <- roc(data_mla$ffr_0.8, mdl_mla_a_bsa_ffr$fitted.values)
cutoff_best <- coords(roc_mla, x = "best", best.method = "youden", ret = "all")
cutoff_trueneg <- coords(roc_mla, x = "best", best.method = "youden", best.weights = c(1, 1/3), ret = "all") # with weights for high true negative
cutoff_truepos <- coords(roc_mla, x = "best", best.method = "youden", best.weights = c(1, 2/3), ret = "all") # with weights for high true positive
plot_log_mla +
    geom_hline(yintercept = cutoff_trueneg$threshold, linetype = "dashed") + 
    geom_hline(yintercept = cutoff_best$threshold, linetype = "solid") + 
    geom_hline(yintercept = cutoff_truepos$threshold, linetype = "dashed")

### width
data_mla <- baseline %>% select(ffr_0.8, ccta_mla_w_bsa) %>% drop_na()
roc_mla_w <- roc(data_mla$ffr_0.8, mdl_mla_w_bsa_ffr$fitted.values)
cutoff_best <- coords(roc_mla_w, x = "best", best.method = "youden", ret = "all")
cutoff_trueneg <- coords(roc_mla_w, x = "best", best.method = "youden", best.weights = c(1, 1/3), ret = "all") # with weights for high true negative
cutoff_truepos <- coords(roc_mla_w, x = "best", best.method = "youden", best.weights = c(1, 2/3), ret = "all") # with weights for high true positive
plot_log_mla_w +
    geom_hline(yintercept = cutoff_trueneg$threshold, linetype = "dashed") + 
    geom_hline(yintercept = cutoff_best$threshold, linetype = "solid") + 
    geom_hline(yintercept = cutoff_truepos$threshold, linetype = "dashed")

### proximal narrowing
data_mla <- baseline %>% select(ffr_0.8, ccta_pn_dist) %>% drop_na()
roc_mla_pn <- roc(data_mla$ffr_0.8, mdl_mla_pn_ffr$fitted.values)
cutoff_best <- coords(roc_mla_pn, x = "best", best.method = "youden", ret = "all")


## invasive
data_ivus <- baseline %>% select(ffr_0.8, inv_ivusrest_mla_bsa) %>% drop_na()
### ivus rest
roc_ivusrest <- roc(data_mla$ffr_0.8, mdl_ivusrest_bsa_ffr$fitted.values)
cutoff_best <- coords(roc_ivusrest, x = "best", best.method = "youden", ret = "all")
cutoff_trueneg <- coords(roc_ivusrest, x = "best", best.method = "youden", best.weights = c(1, 1/3), ret = "all") # with weights for high true negative
cutoff_truepos <- coords(roc_ivusrest, x = "best", best.method = "youden", best.weights = c(1, 2/3), ret = "all") # with weights for high true positive
plot_log_ivusrest +
    geom_hline(yintercept = cutoff_trueneg$threshold, linetype = "dashed") + 
    geom_hline(yintercept = cutoff_best$threshold, linetype = "solid") + 
    geom_hline(yintercept = cutoff_truepos$threshold, linetype = "dashed")

### ivus dobu
data_ivus <- baseline %>% select(ffr_0.8, inv_ivusdobu_ostial_a_bsa) %>% drop_na()
roc_ivus_ostial_a <- roc(data_ivus$ffr_0.8, mdl_ivusdobu_ostial_a_bsa_ffr$fitted.values)
cutoff_best <- coords(roc_ivus_ostial_a, x = "best", best.method = "youden", ret = "all")

### ostial width
data_ivus <- baseline %>% select(ffr_0.8, inv_ivusdobu_ostial_w_bsa) %>% drop_na()
roc_ivus_ostial_w <- roc(data_ivus$ffr_0.8, mdl_ivusdobu_ostial_w_bsa_ffr$fitted.values)
cutoff_best <- coords(roc_ivus_ostial_w, x = "best", best.method = "youden", ret = "all")

### ostial elliptic
data_ivus <- baseline %>% select(ffr_0.8, inv_ivusdobu_ostial_ellip) %>% drop_na()
roc_ivus_ostial_ellip <- roc(data_ivus$ffr_0.8, mdl_ivusdobu_ostial_ellip_ffr$fitted.values)
cutoff_best <- coords(roc_ivus_ostial_ellip, x = "best", best.method = "youden", ret = "all")

### ostial proximal narrowing
data_ivus <- baseline %>% select(ffr_0.8, inv_ivusdobu_ostial_pn) %>% drop_na()
roc_ivus__ostial_pn <- roc(data_ivus$ffr_0.8, mdl_ivusdobu_ostial_pn_ffr$fitted.values)
cutoff_best <- coords(roc_ivus_pn, x = "best", best.method = "youden", ret = "all")

### MLA
data_ivus <- baseline %>% select(ffr_0.8, inv_ivusdobu_mla_bsa) %>% drop_na()
roc_ivus_mla <- roc(data_ivus$ffr_0.8, mdl_ivusdobu_bsa_ffr$fitted.values)
cutoff_best <- coords(roc_ivus_mla, x = "best", best.method = "youden", ret = "all")

### MLA width
data_ivus <- baseline %>% select(ffr_0.8, inv_ivusdobu_mla_w_bsa) %>% drop_na()
roc_ivus_mla_w <- roc(data_ivus$ffr_0.8, mdl_ivusdobu_w_bsa_ffr$fitted.values)
cutoff_best <- coords(roc_ivus_mla_w, x = "best", best.method = "youden", ret = "all")

### MLA proximal narrowing
data_ivus <- baseline %>% select(ffr_0.8, inv_ivusdobu_pn) %>% drop_na()
roc_ivus_mla_pn <- roc(data_ivus$ffr_0.8, mdl_ivusdobu_pn_ffr$fitted.values)
cutoff_best <- coords(roc_ivus_mla_pn, x = "best", best.method = "youden", ret = "all")


ggplot(baseline, aes(inv_ivusrest_ostial_a, ccta_ostial_a)) + 
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_y_continuous(limits = c(2, 12)) +
    scale_x_continuous(limits = c(2, 12))

ggplot(baseline, aes(inv_ivusrest_mla, ccta_mla_a)) + 
    geom_point() + 
    geom_smooth(method = "lm", se = FALSE) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_y_continuous(limits = c(2, 12)) +
    scale_x_continuous(limits = c(2, 12))


# ostial area
roc_list <- list("IVUS, AUC: 0.944" = roc_ivus_ostial_a, "CCTA, AUC: 0.894" = roc_ost_a)
best_ivus <- coords(roc_ivus_ostial_a, x = "best", best.method = "youden", ret = "all")
# weighted <- coords(roc_ivus_ostial_a, x = "best", best.method = "youden", best.weights = c(1, 38/50), ret = "all")
# x_weighted <- 1 - weighted$specificity
# y_weighted <- weighted$sensitivity
x_ivus <- 1 - best_ivus$specificity
y_ivus <- best_ivus$sensitivity
best_ccta <- coords(roc_ost_a, x = "best", best.method = "youden", ret = "all")
x_ccta <- 1 - best_ccta$specificity
y_ccta <- best_ccta$sensitivity
fig_roc_ost_a <- ggroc(roc_list, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Ostial area") +
    geom_point(aes(x = x_ivus, y = y_ivus), color = "red", size = 2) +
    geom_point(aes(x = x_ccta, y = y_ccta), color = "red", size = 2) +
    # geom_point(aes(x = x_weighted, y = y_weighted), color = "red", size = 2) +
    theme_classic()
ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/roc_ostial_a.png", fig_roc_ost_a, width = 6, height = 4)


# ostial elliptic ratio
roc_list <- list("IVUS, AUC: 0.745" = roc_ivus_ostial_ellip, "CCTA, AUC: 0.742" = roc_ost_ellip)
best_ivus <- coords(roc_ivus_ostial_ellip, x = "best", best.method = "youden", ret = "all")
x_ivus <- 1 - best_ivus$specificity
y_ivus <- best_ivus$sensitivity
best_ccta <- coords(roc_ost_ellip, x = "best", best.method = "youden", ret = "all")
x_ccta <- 1 - best_ccta$specificity
y_ccta <- best_ccta$sensitivity
ggroc(roc_list, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Ostial elliptic ratio") +
    geom_point(aes(x = x_ivus, y = y_ivus), color = "red", size = 2) +
    geom_point(aes(x = x_ccta, y = y_ccta), color = "red", size = 2) +
    theme_classic()

# ostial width
roc_list <- list("IVUS, AUC: 0.836" = roc_ivus_ostial_w, "CCTA, AUC: 0.832" = roc_ost_w)
best_ivus <- coords(roc_ivus_ostial_w, x = "best", best.method = "youden", ret = "all")
x_ivus <- 1 - best_ivus$specificity
y_ivus <- best_ivus$sensitivity
best_ccta <- coords(roc_ost_w, x = "best", best.method = "youden", ret = "all")
x_ccta <- 1 - best_ccta$specificity
y_ccta <- best_ccta$sensitivity
fig_roc_ost_w <- ggroc(roc_list, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Ostial width") +
    geom_point(aes(x = x_ivus, y = y_ivus), color = "red", size = 2) +
    geom_point(aes(x = x_ccta, y = y_ccta), color = "red", size = 2) +
    theme_classic()
ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/roc_ostial_w.png", fig_roc_ost_w, width = 6, height = 4)

# ostial proximal narrowing
roc_list <- list("IVUS, AUC: 0.878" = roc_ivus_ostial_pn, "CCTA, AUC: 0.770" = roc_ost_pn)
best_ivus <- coords(roc_ivus_ostial_pn, x = "best", best.method = "youden", ret = "all")
#weighted <- coords(roc_ivus_ostial_pn, x = "best", best.method = "youden", best.weights = c(1, 0.24), ret = "all")
x_ivus <- 1 - best_ivus$specificity
y_ivus <- best_ivus$sensitivity
#x_weighted <- 1 - weighted$specificity
#y_weighted <- weighted$sensitivity
best_ccta <- coords(roc_ost_pn, x = "best", best.method = "youden", ret = "all")
x_ccta <- 1 - best_ccta$specificity
y_ccta <- best_ccta$sensitivity
ggroc(roc_list, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Ostial proximal narrowing") +
    geom_point(aes(x = x_ivus, y = y_ivus), color = "red", size = 2) +
    geom_point(aes(x = x_ccta, y = y_ccta), color = "red", size = 2) +
    #geom_point(aes(x = x_weighted, y = y_weighted), color = "red", size = 2) +
    theme_classic()

# MLA
roc_list <- list("IVUS, AUC: 0.894" = roc_ivus_mla, "CCTA, AUC: 0.790" = roc_mla)
best_ivus <- coords(roc_ivus_mla, x = "best", best.method = "youden", ret = "all")
x_ivus <- 1 - best_ivus$specificity
y_ivus <- best_ivus$sensitivity
best_ccta <- coords(roc_mla, x = "best", best.method = "youden", ret = "all")
x_ccta <- 1 - best_ccta$specificity
y_ccta <- best_ccta$sensitivity
ggroc(roc_list, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("MLA") +
    geom_point(aes(x = x_ivus, y = y_ivus), color = "red", size = 2) +
    geom_point(aes(x = x_ccta, y = y_ccta), color = "red", size = 2) +
    theme_classic()

# MLA width
roc_list <- list("IVUS, AUC: 0.797" = roc_ivus_mla_w, "CCTA, AUC: 0.768" = roc_mla_w)
best_ivus <- coords(roc_ivus_mla_w, x = "best", best.method = "youden", ret = "all")
x_ivus <- 1 - best_ivus$specificity
y_ivus <- best_ivus$sensitivity
best_ccta <- coords(roc_mla_w, x = "best", best.method = "youden", ret = "all")
x_ccta <- 1 - best_ccta$specificity
y_ccta <- best_ccta$sensitivity
ggroc(roc_list, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("MLA width") +
    geom_point(aes(x = x_ivus, y = y_ivus), color = "red", size = 2) +
    geom_point(aes(x = x_ccta, y = y_ccta), color = "red", size = 2) +
    theme_classic()

# MLA proximal narrowing
roc_list <- list("IVUS, AUC: 0.799" = roc_ivus_mla_pn, "CCTA, AUC: 0.753" = roc_mla_pn)
best_ivus <- coords(roc_ivus_mla_pn, x = "best", best.method = "youden", ret = "all")
x_ivus <- 1 - best_ivus$specificity
y_ivus <- best_ivus$sensitivity
best_ccta <- coords(roc_mla_pn, x = "best", best.method = "youden", ret = "all")
x_ccta <- 1 - best_ccta$specificity
y_ccta <- best_ccta$sensitivity
ggroc(roc_list, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("MLA proximal narrowing") +
    geom_point(aes(x = x_ivus, y = y_ivus), color = "red", size = 2) +
    geom_point(aes(x = x_ccta, y = y_ccta), color = "red", size = 2) +
    theme_classic()


# test different thresholds
list_prediction_data <- list("CCTA ostial area (BSA)" = all_ost_a_bsa[[2]], 
                            "IVUS ostial area (BSA)" = all_ivusdobu_ostial_a_bsa[[2]],
                            "CCTA ostial elliptic" = all_ost_ellip[[2]],
                            "IVUS ostial elliptic" = all_ivusrest_ostial_ellip[[2]],
                            "CCTA ostial width (BSA)" = all_ost_w_bsa[[2]],
                            "IVUS ostial width (BSA)" = all_ivusdobu_ostial_w_bsa_ffr[[2]],
                            "CCTA ostial proximal narrowing" = all_ost_pn[[2]],
                            "IVUS ostial proximal narrowing" = all_ivusdobu_ostial_pn[[2]],
                            "CCTA MLA (BSA)" = all_mla_a_bsa[[2]],
                            "IVUS MLA (BSA)" = all_ivusdobu_bsa[[2]],
                            "CCTA MLA width (BSA)" = all_mla_w_bsa[[2]],
                            "IVUS MLA width (BSA)" = all_ivusdobu_w_bsa_ffr[[2]],
                            "CCTA MLA proximal narrowing" = all_mla_pn[[2]],
                            "IVUS MLA proximal narrowing" = all_ivusdobu_pn[[2]])
list_roc <- list("CCTA ostial area (BSA)" = roc_ost_a, 
                "IVUS ostial area (BSA)" = roc_ivus_ostial_a,
                "CCTA ostial elliptic" = roc_ost_ellip,
                "IVUS ostial elliptic" = roc_ivus_ostial_ellip,
                "CCTA ostial width (BSA)" = roc_ost_w,
                "IVUS ostial width (BSA)" = roc_ivus_ostial_w,
                "CCTA ostial proximal narrowing" = roc_ost_pn,
                "IVUS ostial proximal narrowing" = roc_ivus_ostial_pn,
                "CCTA MLA (BSA)" = roc_mla,
                "IVUS MLA (BSA)" = roc_ivus_mla,
                "CCTA MLA width (BSA)" = roc_mla_w,
                "IVUS MLA width (BSA)" = roc_ivus_mla_w,
                "CCTA MLA proximal narrowing" = roc_mla_pn,
                "IVUS MLA proximal narrowing" = roc_ivus_mla_pn)

roc_table_creater <- function(list_roc, list_prediction_data, wgt) {
    roc_table <- data.frame()
    for (i in 1:length(list_roc)) {
        roc <- list_roc[[i]]
        prediction_data <- list_prediction_data[[i]]
        best_coords <- coords(roc, x = "best", best.method = "youden", ret = "all")
        wgt_coords <- coords(roc, x = "best", best.method = "youden", best.weights = wgt, ret = "all")

        value_cutoff <- prediction_data %>% filter(prediction_data[,2] > best_coords$threshold)
        value_cutoff_wgt <- prediction_data %>% filter(prediction_data[,2] > wgt_coords$threshold)
        
        roc_table <- rbind(roc_table, data.frame(
            "Predictor" = names(list_roc)[i],
            "AUC" = round(roc$auc, 3),
            "Best threshold" = round(best_coords$threshold, 3),
            "Accuracy" = round(best_coords$accuracy*100, 0),
            "Best sensitivity" = round(best_coords$sensitivity *100, 0),
            "Best specificity" = round(best_coords$specificity*100, 0),
            "Best PPV" = round(best_coords$ppv*100, 0),
            "Best NPV" = round(best_coords$npv*100, 0),
            "Cutoff value" = round(max(value_cutoff[,1]), 2),
            "Weighted threshold" = round(wgt_coords$threshold, 3),
            "Weighted accuracy" = round(wgt_coords$accuracy*100, 0),
            "Weighted sensitivity" = round(wgt_coords$sensitivity*100, 0),
            "Weighted specificity" = round(wgt_coords$specificity*100, 0),
            "Weighted PPV" = round(wgt_coords$ppv*100, 0),
            "Weighted NPV" = round(wgt_coords$npv*100, 0),
            "Weighted cutoff value" = round(max(value_cutoff_wgt[,1]),2)
        ))
    }
    return(roc_table)
}

roc_table <- roc_table_creater(list_roc, list_prediction_data, c(1, 38/50))
write_csv(roc_table, "roc_table.csv")













# ########################################################################################
# data <- baseline %>% select(ffr_0.8, ccta_ostial_a) %>% drop_na()
# roc_ost_a <- roc(data$ffr_0.8, mdl_ost_a_ffr$fitted.values)
# coords(roc_ost_a, x = "best", best.method = "youden", best.weights = c(1, 0.2), ret = "all") # with weights for high true negative
# coords(roc_ost_a, x = "best", best.method = "youden", best.weights = c(1, 0.8), ret = "all") # with weights for high true positive
# coords(roc_ost_a, x = "best", best.method = "youden", ret = "all")
# plot_log_ost_a +
#     geom_hline(yintercept = 0.52, linetype = "dashed") + 
#     geom_hline(yintercept = 0.24, linetype = "solid") + 
#     geom_hline(yintercept = 0.05, linetype = "dashed")

# # double check that elliptic is not significant
# all_ost_ellip <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_elliptic")
# mdl_ost_ellip_ffr <- all_ost_ellip[[1]]
# summary(mdl_ost_ellip_ffr) # significant
# prediction_data_ost_ellip <- all_ost_ellip[[2]]
# print(prediction_data_ost_ellip, n = 50)


# # ostium proximal narrowing
# all_ost_pn <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_pn")
# mdl_ost_pn_ffr <- all_ost_pn[[1]]
# summary(mdl_ost_pn_ffr) # significant
# prediction_data_ost_pn <- all_ost_pn[[2]]
# print(prediction_data_ost_pn, n = 50)

# # ostium bsa adjusted
# all_ost_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_a_bsa")
# mdl_ost_bsa_ffr <- all_ost_bsa[[1]]
# summary(mdl_ost_bsa_ffr) # significant
# prediction_data_ost_bsa <- all_ost_bsa[[2]]
# print(prediction_data_ost_bsa, n = 50)
# plot_log_ost_bsa <- all_ost_bsa[[3]]
# plot_log_ost_bsa
# plot_odds_ost_bsa <- all_ost_bsa[[4]]
# plot_odds_ost_bsa
# data <- baseline %>% select(ffr_0.8, ccta_ostial_a_bsa) %>% drop_na()
# roc_ost_bsa <- roc(data$ffr_0.8, mdl_ost_bsa_ffr$fitted.values)
# coords(roc_ost_bsa, x = "best", best.method = "youden", best.weights = c(1, 0.2), ret = "all") # with weights for high true negative
# coords(roc_ost_bsa, x = "best", best.method = "youden", best.weights = c(1, 0.8), ret = "all") # with weights for high true positive
# coords(roc_ost_bsa, x = "best", best.method = "youden", ret = "all")
# plot_log_ost_bsa +
#     geom_hline(yintercept = 0.52, linetype = "dashed") + 
#     geom_hline(yintercept = 0.24, linetype = "solid") + 
#     geom_hline(yintercept = 0.05, linetype = "dashed")
# ggroc(list(roc_ost_a, roc_ost_bsa), legacy.axes = TRUE, lwd = 3) +
#     scale_color_manual(values = c("#1818aa", "#18aa7e")) +
#     theme(legend.position = "none") +
#     theme_classic()

# # mla
# all_mla <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_a")
# mdl_mla_ffr <- all_mla[[1]]
# summary(mdl_mla_ffr) # significant
# prediction_data_mla <- all_mla[[2]]
# print(prediction_data_mla, n = 50)
# plot_log_mla <- all_mla[[3]]
# plot_log_mla
# plot_odds_mla <- all_mla[[4]]
# plot_odds_mla
# confusion_mla <- all_mla[[5]]
# confusion_mla
# summary(confusion_mla, event_level = "second")
# plot_acc_mla <- all_mla[[6]]
# plot_acc_mla
# data <- baseline %>% select(ffr_0.8, ccta_mla_a) %>% drop_na()
# roc_mla <- roc(data$ffr_0.8, mdl_mla_ffr$fitted.values)

# # mla bsa adjusted
# all_mla_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_a_bsa")
# mdl_mla_bsa_ffr <- all_mla_bsa[[1]]
# summary(mdl_mla_bsa_ffr) # significant
# prediction_data_mla_bsa <- all_mla_bsa[[2]]
# print(prediction_data_mla_bsa, n = 50)
# plot_log_mla_bsa <- all_mla_bsa[[3]]
# plot_log_mla_bsa
# data <- baseline %>% select(ffr_0.8, ccta_mla_a_bsa) %>% drop_na()
# roc_mla_bsa <- roc(data$ffr_0.8, mdl_mla_bsa_ffr$fitted.values)
# coords(roc_mla_bsa, x = "best", best.method = "youden", best.weights = c(1, 0.2), ret = "all") # with weights for high true negative
# coords(roc_mla_bsa, x = "best", best.method = "youden", best.weights = c(1, 0.8), ret = "all") # with weights for high true positive
# coords(roc_mla_bsa, x = "best", best.method = "youden", ret = "all")
# plot_log_mla_bsa +
#     geom_hline(yintercept = 0.52, linetype = "dashed") + 
#     geom_hline(yintercept = 0.24, linetype = "solid") + 
#     geom_hline(yintercept = 0.05, linetype = "dashed")
# ggroc(list(roc_ost_a, roc_ost_bsa, roc_mla, roc_mla_bsa), legacy.axes = TRUE, lwd = 2) +
#     scale_color_manual(values = c("#1818aa", "#18aa7e", "#aa1818", "orange")) +
#     theme(legend.position = "none") +
#     theme_classic()

# # exploration of ccta_mla_elliptic
# all_ellip <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_elliptic")
# mdl_ellip <- all_ellip[[1]] 
# summary(mdl_ellip) # significant
# prediction_data_ellip <- all_ellip[[2]]
# print(prediction_data_ellip, n = 50)
# plot_log_ellip <- all_ellip[[3]]
# plot_odds_ellip <- all_ellip[[4]]
# plot_odds_ellip
# confusion_ellip <- all_ellip[[5]]
# summary(confusion_ellip, event_level = "second")
# plot_acc <- all_ellip[[6]]


# # exploration of ccta_pn_dist
# all_pn <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_pn_dist")
# mdl_pn_ffr <- all_pn[[1]]
# summary(mdl_pn_ffr) # significant
# prediction_data_pn <- all_pn[[2]]
# print(prediction_data_pn, n = 50)
# plot_log_pn <- all_pn[[3]]
# plot_log_pn
# plot_odds_pn <- all_pn[[4]]
# plot_odds_pn
# confusion_pn <- all_pn[[5]]
# summary(confusion_pn, event_level = "second")
# plot_acc_pn <- all_pn[[6]]
# plot_acc_pn


# # exploration of ccta_imc_length
# all_imc <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_imc_length")
# mdl_imc_ffr <- all_imc[[1]]
# summary(mdl_imc_ffr) # not significant
# prediction_data_imc <- all_imc[[2]]
# print(prediction_data_imc, n = 50)
# plot_log_imc <- all_imc[[3]]
# plot_log_imc
# plot_odds_imc <- all_imc[[4]]
# plot_odds_imc
# confusion_imc <- all_imc[[5]]
# confusion_imc
# summary(confusion_imc, event_level = "second")
# plot_acc_imc <- all_imc[[6]]
# plot_acc_imc


# # exploration of ccta_aa_degree
# all_aa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_aa_degree")
# mdl_aa_ffr <- all_aa[[1]]
# summary(mdl_aa_ffr) # not significant
# prediction_data_aa <- all_aa[[2]]
# print(prediction_data_aa, n = 50)
# plot_log_aa <- all_aa[[3]]
# plot_log_aa
# plot_odds_aa <- all_aa[[4]]
# plot_odds_aa
# confusion_aa <- all_aa[[5]]
# summary(confusion_aa, event_level = "second")
# plot_acc_aa <- all_aa[[6]]
# plot_acc_aa




# # save images
# ggsave(paste(yaml$demographics$output_dir_figures, "/forest_continuous_ffr0.8.png"), forest_continuous_ffr0.8, width = 10, height = 5)

# ####################################################################################################################
# # binary explanatory variable
# categorical <- baseline %>%
#     mutate(
#         caa_slo = ifelse(caa_slo == "yes", 1, 0),
#         caa_elliptic = ifelse(caa_elliptic == "yes", 1, 0),
#         caa_pn = ifelse(caa_pn == "yes", 1, 0),
#         caa_im = ifelse(caa_im == "yes", 1, 0),
#         caa_angle = ifelse(caa_angle == "yes", 1, 0)
#     )

# # slit like ostium
# all <- logistic_regression_categorical(categorical, "ffr_0.8", "caa_slo")
# mdl <- all[[1]]
# odds_ratio <- all[[2]]
# lower_ci <- all[[3]]
# upper_ci <- all[[4]]
# summary(mdl) # significant

# # elliptic vessel shape
# all <- logistic_regression_categorical(categorical, "ffr_0.8", "caa_elliptic")
# mdl <- all[[1]]
# summary(mdl) # not significant

# # proximal narrowing
# all <- logistic_regression_categorical(categorical, "ffr_0.8", "caa_pn")
# mdl <- all[[1]]
# summary(mdl) # not significant

# # intramural course
# all <- logistic_regression_categorical(categorical, "ffr_0.8", "caa_im")
# mdl <- all[[1]]
# summary(mdl) # not significant

# # acute angle
# all <- logistic_regression_categorical(categorical, "ffr_0.8", "caa_angle")
# mdl <- all[[1]]
# summary(mdl) # not significant

# # forest plots
# vars = c("caa_slo", "caa_pn")
# forest_df <- data.frame(var = character(), odds_ratio = numeric(), lower_ci = numeric(), upper_ci = numeric(), stringsAsFactors = FALSE)
# for (var in vars){
#     all <- logistic_regression_categorical(categorical, "ffr_0.8", var)
#     mdl <- all[[1]]
#     odds_ratio <- all[[2]]
#     lower_ci <- all[[3]]
#     upper_ci <- all[[4]]
#     forest_df <- rbind(forest_df, data.frame(var = var, odds_ratio = odds_ratio, lower_ci = lower_ci, upper_ci = upper_ci))
# }


# combined_mdl <- glm(ffr_0.8 ~ ccta_ostial_a_bsa + ccta_ostial_elliptic + ccta_ostial_w_bsa + ccta_mla_a_bsa + ccta_mla_elliptic + ccta_mla_w_bsa, data = baseline, family = "binomial")
# vif(combined_mdl)
# #create vector of VIF values
# vif_values <- vif(combined_mdl_mdl)
# #create horizontal bar chart to display each VIF value
# barplot(vif_values, main = "VIF Values", horiz = TRUE, col = "steelblue")
# #add vertical line at 5
# abline(v = 5, lwd = 3, lty = 2)

# # plot correlation matrix
# library(corrplot)
# correlation_matrix <- cor(baseline %>% select(ccta_ostial_a_bsa, ccta_ostial_elliptic, ccta_ostial_w_bsa, ccta_mla_a_bsa, ccta_mla_elliptic, ccta_mla_w_bsa) %>% drop_na())
# corrplot(correlation_matrix, method = "color", type = "upper", order = "hclust", tl.col = "black", tl.srt = 45)

# # added variables to decrease the vif
# data <- baseline %>% select(ffr_0.8, ccta_ostial_a_bsa, ccta_ostial_elliptic, ccta_ostial_w_bsa, ccta_mla_a_bsa, ccta_mla_elliptic, ccta_mla_w_bsa)
# data <- data %>% drop_na()
# data <- data %>% mutate(
#     ostium_ellip = ccta_ostial_elliptic / ccta_ostial_w_bsa,
#     mla_ellip = ccta_mla_elliptic / ccta_mla_w_bsa
# )

# decreased_mdl <- glm(ffr_0.8 ~ ccta_ostial_a_bsa + ostium_ellip + ccta_mla_a_bsa + mla_ellip, data = data, family = "binomial")
# vif(decreased_mdl)