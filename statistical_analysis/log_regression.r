library(tidyverse)
library(broom)
library(ggplot2)
library(yardstick)

baseline <- readRDS("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/baseline.rds")

# quick data prepping
baseline <- baseline %>% 
    mutate(
        ccta_ostial_pn = ccta_ostial_a / ccta_dist_a,
        ccta_slo_num = ccta_ostial_elliptic / ccta_ostial_pn,
        ccta_quali = ifelse(ccta_photon == "yes", 3, ifelse(is.na(ccta_quali) & ccta_photon == "no", 2, ccta_quali)),
        ffr_0.8 = ifelse(ffr_0.8 == "yes", 1, 0)
    )

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
        !!sym(explanatory) := seq(min - buffer, max + buffer, length.out = 50)
    )

    prediction_data <- explanatory_data %>% 
        mutate(
            response = predict(mdl, explanatory_data, type = "response"),
            most_likely_outcome = round(response),
            odds_ratio = response / (1 - response),
            log_odds_ratio = log(odds_ratio),
            log_odds_ratio2 = predict(mdl, explanatory_data)
        )

    plot_log <- ggplot(data, aes(x = !!sym(explanatory), y = !!sym(response))) +
        geom_point() + 
        geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) +
        scale_x_continuous(breaks = seq(min - buffer, max + buffer, length.out = 50))

    plot_odds <- ggplot(prediction_data, aes(x = !!sym(explanatory), y = odds_ratio)) + 
        geom_point() + 
        geom_line() +
        geom_hline(yintercept = 1, linetype = "dashed") +
        scale_x_continuous(breaks = seq(30, 90, 10)) +
        scale_y_continuous(breaks = seq(0, 8, 1))

    actual_response <- data[[response]]
    predicted_response <- round(fitted(mdl))
    outcome <- table(predicted_response, actual_response)

    confusion <- conf_mat(outcome)
    plot_acc <- autoplot(confusion)

    return(list(mdl, prediction_data, plot_log, plot_odds, confusion, plot_acc))
}

# exploration of ccta_slo_num, ccta_mla_elliptic, ccta_imc_length, ccta_pn_dist, ccta_aa_degree
all <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_slo_num")
mdl <- all[[1]]$coefficients
prediction_data <- all[[2]]
plot_log <- all[[3]]
plot_odds <- all[[4]]
confusion <- all[[5]]
plot_acc <- all[[6]]



# exploration of ccta_mla_elliptic

# exploration of ccta_imc_length

# exploration of ccta_pn_dist

# exploration of ccta_aa_degree







# # logistic regression ccta_mla_a
# data <- baseline %>% 
#     select(ccta_pn_dist, ffr_0.8) %>%
#     drop_na() %>%
#     mutate(ffr_0.8 = ifelse(ffr_0.8 == "yes", 1, 0))

# mdl_mla_ffr <- glm(ffr_0.8 ~ ccta_pn_dist, data = data, family = binomial)

# explanatory_data <- tibble(
#     ccta_pn_dist = seq(0, 100, 2)
# )

# prediction_data <- explanatory_data %>% 
#     mutate(
#         ffr_0.8 = predict(mdl_mla_ffr, explanatory_data, type = "response"),
#         most_likely_outcome = round(ffr_0.8),
#         odds_ratio = ffr_0.8 / (1 - ffr_0.8),
#         log_odds_ratio = log(odds_ratio),
#         log_odds_ratio2 = predict(mdl_mla_ffr, explanatory_data)
#     )

# ggplot(data, aes(x = ccta_pn_dist, y = ffr_0.8)) + 
#     geom_point() + 
#     geom_smooth(method = "glm", method.args = list(family = "binomial"), se = FALSE) +
#     scale_x_continuous(breaks = seq(0, 100, 2))

# # plot odds ratio
# ggplot(prediction_data, aes(x = ccta_pn_dist, y = odds_ratio)) + 
#     geom_point() + 
#     geom_line() +
#     geom_hline(yintercept = 1, linetype = "dashed") +
#     scale_x_continuous(breaks = seq(30, 90, 10)) +
#     scale_y_continuous(breaks = seq(0, 8, 1))

# # plot log odds ratio
# ggplot(prediction_data, aes(x = ccta_pn_dist, y = log_odds_ratio)) + 
#     geom_point() + 
#     geom_line() +
#     geom_hline(yintercept = 0, linetype = "dashed") +
#     scale_x_continuous(breaks = seq(30, 90, 10)) +
#     scale_y_continuous(breaks = seq(-4, 4, 1))

# actual_response <- data$ffr_0.8
# predicted_response <- round(fitted(mdl_mla_ffr))
# outcome <- table(predicted_response, actual_response)

# confusion <- conf_mat(outcome)
# autoplot(confusion)

# summary(confusion, event_level = "second")
# # pull results
# # accuracy is TP + TN / (TP + TN + FP + FN)
# summary(confusion) %>% slice(1)
# # sensitivity is TP / (TP + FN)
# summary(confusion) %>% slice(3)
# # specificity is TN / (TN + FP)
# summary(confusion) %>% slice(4)

# # multiple logistic regression
# # new model
# data <- baseline %>% 
#     select(caa_slo, caa_elliptic, caa_im, caa_pn, caa_angle, ffr_0.8) %>%
#     drop_na() %>%
#     mutate(
#         ffr_0.8 = ifelse(ffr_0.8 == "yes", 1, 0)
#     )
# mdl_hrf_ffr <- glm(ffr_0.8 ~ caa_slo + caa_elliptic + caa_im + caa_pn + caa_angle, data = data, family = binomial)
# explanatory_data <- tibble(
#     caa_slo = levels(data$caa_slo),
#     caa_elliptic = levels(data$caa_elliptic),
#     caa_im = levels(data$caa_im),
#     caa_pn = levels(data$caa_pn),
#     caa_angle = levels(data$caa_angle)
# )

# prediction_data <- explanatory_data %>% 
#     mutate(
#         ffr_0.8 = predict(mdl_hrf_ffr, explanatory_data, type = "response"),
#     )

# actual_response <- data$ffr_0.8
# predicted_response <- round(fitted(mdl_hrf_ffr))
# outcome <- table(predicted_response, actual_response)
# confusion <- conf_mat(outcome)
# autoplot(confusion)
# # accuracy
# summary(confusion, event_level = "second") %>% slice(1)
# # sensitivity
# summary(confusion) %>% slice(3)
# # specificity
# summary(confusion) %>% slice(4)

# # multiple logistic regression
# # new model
# data <- baseline %>% 
#     select(ccta_quali, ccta_photon, ccta_slo_num, ccta_mla_elliptic, ccta_imc_length, ccta_pn_dist, ccta_aa_degree, ffr_0.81) %>%
#     drop_na() %>%
#     mutate(
#         ffr_0.81 = ifelse(ffr_0.81 == "yes", 1, 0)
#     )
# mdl_hrf_ffr <- glm(ffr_0.81 ~ ccta_slo_num + ccta_mla_elliptic + ccta_imc_length + ccta_pn_dist + ccta_aa_degree, data = data, family = binomial)
# explanatory_data <- tibble(
#     ccta_slo_num = seq(0.5, 25, 0.5),
#     ccta_mla_elliptic = seq(0, 4.9, 0.1),
#     ccta_imc_length = seq(0, 24.5, 0.5),
#     ccta_pn_dist = seq(5, 127.5, 2.5),
#     ccta_aa_degree = seq(0, 49, 1)
# )

# prediction_data <- explanatory_data %>% 
#     mutate(
#         ffr_0.81 = predict(mdl_hrf_ffr, explanatory_data, type = "response"),
#     )

# actual_response <- data$ffr_0.81
# predicted_response <- round(fitted(mdl_hrf_ffr))
# outcome <- table(predicted_response, actual_response)
# confusion <- conf_mat(outcome)
# autoplot(confusion)
# # accuracy
# summary(confusion, event_level = "second") %>% slice(1)
# # sensitivity
# summary(confusion) %>% slice(3)
# # specificity
# summary(confusion) %>% slice(4)

# ggplot(data, aes(x = ccta_pn_dist)) +
#     geom_histogram() +
#     facet_grid(rows = vars(ffr_0.81))

# str(mdl_hrf_ffr)