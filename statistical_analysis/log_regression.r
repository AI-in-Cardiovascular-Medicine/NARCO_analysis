library(tidyverse)
library(broom)
library(ggplot2)
library(yardstick)
library(pROC)
library(randomForest)
library(car)
library(readxl)

baseline <- readRDS("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/baseline.rds")
additional_ccta <- read_excel("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/additional_ccta_data.xlsx")

# keep only columns record_id, length_sinus, bogen_rca, bogen_lca, scalor_width, scalor_height
additional_ccta <- additional_ccta %>% 
    select(record_id, length_sinus, bogen_rca, bogen_lca, scalor_width, scalor_height)

# merge dataframes by record_id
baseline <- baseline %>% 
    left_join(additional_ccta, by = "record_id")

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
        ccta_stj_rca_bsa = ccta_stj_rca / bsa,
        bogen_rca_bsa = bogen_rca / bsa,
        ccta_stj_rca_scaled = ccta_stj_rca * scalor_height,
        bogen_rca_bsa_scaled = bogen_rca * scalor_height,
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

# rest of the variables
all_imc_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_imc_length_bsa")
mdl_imc_bsa_ffr <- all_imc_bsa[[1]] # not significant!!!
all_aa_bsa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_aa_degree_bsa")
mdl_aa_bsa_ffr <- all_aa_bsa[[1]] # not significant!!!

all_ost_y_pos <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_stj_rca_bsa")
mdl_ost_y_pos_ffr <- all_ost_y_pos[[1]] # not significant!!!
all_ost_y_pos_scaled <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_stj_rca_scaled")
mdl_ost_y_pos_scaled_ffr <- all_ost_y_pos_scaled[[1]] # not significant!!!

all_ost_x_pos <- simple_logistic_regression(baseline, "ffr_0.8", "bogen_rca_bsa")
mdl_ost_x_pos_ffr <- all_ost_x_pos[[1]] # not significant!!!
all_ost_x_pos_scaled <- simple_logistic_regression(baseline, "ffr_0.8", "bogen_rca_bsa_scaled")
mdl_ost_x_pos_scaled_ffr <- all_ost_x_pos_scaled[[1]] # not significant!!!

# forest plot
# bsa adjusted
forest_df_bsa <- data.frame(
    var = c("Ostial area (BSA) [mm2]", "Ostial elliptic ratio", "Ostial major axis (BSA) [mm]", "Ostial minor axis (BSA) [mm]", 
            "Ostial proximal narrowing [%]", "MLA (BSA) [mm2]", "ML elliptic ratio", "ML major axis (BSA) [mm]", "ML minor axis (BSA) [mm]", 
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
    geom_segment(aes(x = 6, xend = 7, y = 2, yend = 2), linetype = "dashed") +
    geom_segment(aes(x = 6, xend = 7, y = 8, yend = 8), linetype = "dashed") +
    theme_classic()

ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/forest_ccta_ffr.png", forest_cut_ffr_bsa, width = 10, height = 6)

# regression models #################################################################################################
# exploration of continuous variables that were significant
## ostium
### area
summary(mdl_ost_a_bsa_ffr) # significant
prediction_data_ost_a <- all_ost_a_bsa[[2]]
print(prediction_data_ost_a, n = 50)
plot_log_ost_a <- all_ost_a_bsa[[3]]
plot_odds_ost_a <- all_ost_a_bsa[[4]]
confusion_ost_a <- all_ost_a_bsa[[5]]
summary(confusion_ost_a, event_level = "second")
plot_acc_ost_a <- all_ost_a_bsa[[6]]

### elliptic
summary(mdl_ost_ellip_ffr) # significant
prediction_data_ost_ellip <- all_ost_ellip[[2]]
print(prediction_data_ost_ellip, n = 50)
plot_log_ost_ellip <- all_ost_ellip[[3]]
plot_odds_ost_ellip <- all_ost_ellip[[4]]

### width
summary(mdl_ost_w_bsa_ffr) # significant
prediction_data_ost_w <- all_ost_w_bsa[[2]]
print(prediction_data_ost_w, n = 50)
plot_log_ost_w <- all_ost_w_bsa[[3]]
plot_odds_ost_w <- all_ost_w_bsa[[4]]

### proximal narrowing
summary(mdl_ost_pn_ffr) # not significant
prediction_data_ost_pn <- all_ost_pn[[2]]
print(prediction_data_ost_pn, n = 50)
plot_log_ost_pn <- all_ost_pn[[3]]
plot_odds_ost_pn <- all_ost_pn[[4]]

## mla
### area
summary(mdl_mla_a_bsa_ffr) # significant
prediction_data_mla <- all_mla_a_bsa[[2]]
print(prediction_data_mla, n = 50)
plot_log_mla <- all_mla_a_bsa[[3]]
plot_odds_mla <- all_mla_a_bsa[[4]]
confusion_mla <- all_mla_a_bsa[[5]]
summary(confusion_mla, event_level = "second")
plot_acc_mla <- all_mla_a_bsa[[6]]

### elliptic
summary(mdl_mla_ellip_ffr) # significant
prediction_data_mla_ellip <- all_mla_ellip[[2]]
print(prediction_data_mla_ellip, n = 50)
plot_log_mla_ellip <- all_mla_ellip[[3]]
plot_odds_mla_ellip <- all_mla_ellip[[4]]

### width
summary(mdl_mla_w_bsa_ffr) # significant
prediction_data_mla_w <- all_mla_w_bsa[[2]]
print(prediction_data_mla, n = 50)
plot_log_mla_w <- all_mla_w_bsa[[3]]
plot_odds_mla_w <- all_mla_w_bsa[[4]]

### proximal narrowing
summary(mdl_mla_pn_ffr) # not significant
prediction_data_mla_pn <- all_mla_pn[[2]]
print(prediction_data_mla, n = 50)
plot_log_mla_pn <- all_mla_pn[[3]]
plot_odds_mla_pn <- all_mla_pn[[4]]

# roc analysis ######################################################################################################
## ostium
### area
data_ostium <- baseline %>% select(ffr_0.8, ccta_ostial_a_bsa) %>% drop_na()
roc_ost_a <- roc(data_ostium$ffr_0.8, mdl_ost_a_bsa_ffr$fitted.values)
cutoff_best <- coords(roc_ost_a, x = "best", best.method = "youden", ret = "all")
cutoff_trueneg <- coords(roc_ost_a, x = "best", best.method = "youden", best.weights = c(1, 12/50), ret = "all") # with weights for high true negative
cutoff_truepos <- coords(roc_ost_a, x = "best", best.method = "youden", best.weights = c(1, 38/50), ret = "all") # with weights for high true positive
plot_log_ost_a +
    geom_hline(yintercept = cutoff_trueneg$threshold, linetype = "dashed") + 
    geom_hline(yintercept = cutoff_best$threshold, linetype = "solid") + 
    geom_hline(yintercept = cutoff_truepos$threshold, linetype = "dashed")

sensitivity <- roc_ost_a$sensitivities
specificity <- 1 - roc_ost_a$specificities
distances <- sqrt((1 - specificity)^2 + sensitivity^2)
optimal_index <- which.min(distances)
optimal_cutoff <- roc_ost_a$thresholds[optimal_index]

x_youden <- 1 - cutoff_best$specificity
y_youden <- cutoff_best$sensitivity
x_weighted <- 1 - cutoff_truepos$specificity
y_weighted <- cutoff_truepos$sensitivity
x_min_dist <- 1 - roc_ost_a$specificities[optimal_index]
y_min_dist <- roc_ost_a$sensitivities[optimal_index]
fig_roc_ost_a <- ggroc(roc_ost_a, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Ostial area") +
    geom_point(aes(x = x_youden, y = y_youden, color = "youden"), size = 2) +
    geom_point(aes(x = x_weighted, y = y_weighted, color = "weighted"), size = 2) +
    scale_color_manual(values = c("youden" = "darkred", "weighted" = "darkblue")) +
    theme_classic()
ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/roc_ostial_a.png", fig_roc_ost_a, width = 6, height = 4)


### elliptic
data_ostium <- baseline %>% select(ffr_0.8, ccta_ostial_elliptic) %>% drop_na()
roc_ost_ellip <- roc(data_ostium$ffr_0.8, mdl_ost_ellip_ffr$fitted.values)
cutoff_best <- coords(roc_ost_ellip, x = "best", best.method = "youden", ret = "all")
cutoff_truepos <- coords(roc_ost_ellip, x = "best", best.method = "youden", best.weights = c(1, 38/50), ret = "all") # with weights for high true positive

x_youden <- 1 - cutoff_best$specificity
y_youden <- cutoff_best$sensitivity
x_weighted <- 1 - cutoff_truepos$specificity
y_weighted <- cutoff_truepos$sensitivity

fig_roc_ost_ellip <- ggroc(roc_ost_ellip, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Ostial elliptic") +
    geom_point(aes(x = x_youden, y = y_youden, color = "youden"), size = 2) +
    geom_point(aes(x = x_weighted, y = y_weighted, color = "weighted"), size = 2) +
    scale_color_manual(values = c("youden" = "darkred", "weighted" = "darkblue")) +
    theme_classic()
ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/roc_ost_ellip.png", fig_roc_ost_ellip, width = 6, height = 4)

### width
data_ostium <- baseline %>% select(ffr_0.8, ccta_ostial_w_bsa) %>% drop_na()
roc_ost_w <- roc(data_ostium$ffr_0.8, mdl_ost_w_bsa_ffr$fitted.values)
cutoff_best <- coords(roc_ost_w, x = "best", best.method = "youden", ret = "all")
cutoff_truepos <- coords(roc_ost_w, x = "best", best.method = "youden", best.weights = c(1, 38/50), ret = "all") # with weights for high true positive

x_youden <- 1 - cutoff_best$specificity
y_youden <- cutoff_best$sensitivity
x_weighted <- 1 - cutoff_truepos$specificity
y_weighted <- cutoff_truepos$sensitivity

fig_roc_ost_w <- ggroc(roc_ost_w, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Ostial minor axis") +
    geom_point(aes(x = x_youden, y = y_youden, color = "youden"), size = 2) +
    geom_point(aes(x = x_weighted, y = y_weighted, color = "weighted"), size = 2) +
    scale_color_manual(values = c("youden" = "darkred", "weighted" = "darkblue")) +
    theme_classic()
ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/roc_ost_w.png", fig_roc_ost_w, width = 6, height = 4)

### proximal narrowing
data_ostium <- baseline %>% select(ffr_0.8, ccta_ostial_pn) %>% drop_na()
roc_ost_pn <- roc(data_ostium$ffr_0.8, mdl_ost_pn_ffr$fitted.values)
cutoff_best <- coords(roc_ost_pn, x = "best", best.method = "youden", ret = "all")
cutoff_truepos <- coords(roc_ost_pn, x = "best", best.method = "youden", best.weights = c(1, 38/50), ret = "all") # with weights for high true positive

x_youden <- 1 - cutoff_best$specificity
y_youden <- cutoff_best$sensitivity
x_weighted <- 1 - cutoff_truepos$specificity
y_weighted <- cutoff_truepos$sensitivity

fig_roc_ost_pn <- ggroc(roc_ost_pn, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Ostial proximal narrowing") +
    geom_point(aes(x = x_youden, y = y_youden, color = "youden"), size = 2) +
    geom_point(aes(x = x_weighted, y = y_weighted, color = "weighted"), size = 2) +
    scale_color_manual(values = c("youden" = "darkred", "weighted" = "darkblue")) +
    theme_classic()
ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/roc_ost_pn.png", fig_roc_ost_pn, width = 6, height = 4)

## mla
data_mla <- baseline %>% select(ffr_0.8, ccta_mla_a_bsa) %>% drop_na()
### area
roc_mla <- roc(data_mla$ffr_0.8, mdl_mla_a_bsa_ffr$fitted.values)
cutoff_best <- coords(roc_mla, x = "best", best.method = "youden", ret = "all")
cutoff_truepos <- coords(roc_mla, x = "best", best.method = "youden", best.weights = c(1, 38/50), ret = "all") # with weights for high true positive

x_youden <- 1 - cutoff_best$specificity
y_youden <- cutoff_best$sensitivity
x_weighted <- 1 - cutoff_truepos$specificity
y_weighted <- cutoff_truepos$sensitivity

fig_roc_mla <- ggroc(roc_mla, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Minimal lumen area") +
    geom_point(aes(x = x_youden, y = y_youden, color = "youden"), size = 2) +
    geom_point(aes(x = x_weighted, y = y_weighted, color = "weighted"), size = 2) +
    scale_color_manual(values = c("youden" = "darkred", "weighted" = "darkblue")) +
    theme_classic()
ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/roc_mla.png", fig_roc_mla, width = 6, height = 4)

### width
data_mla <- baseline %>% select(ffr_0.8, ccta_mla_w_bsa) %>% drop_na()
roc_mla_w <- roc(data_mla$ffr_0.8, mdl_mla_w_bsa_ffr$fitted.values)
cutoff_best <- coords(roc_mla_w, x = "best", best.method = "youden", ret = "all")
cutoff_truepos <- coords(roc_mla_w, x = "best", best.method = "youden", best.weights = c(1, 38/50), ret = "all") # with weights for high true positive

x_youden <- 1 - cutoff_best$specificity
y_youden <- cutoff_best$sensitivity
x_weighted <- 1 - cutoff_truepos$specificity
y_weighted <- cutoff_truepos$sensitivity

fig_roc_mla_w <- ggroc(roc_mla_w, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("MLA minor axis") +
    geom_point(aes(x = x_youden, y = y_youden, color = "youden"), size = 2) +
    geom_point(aes(x = x_weighted, y = y_weighted, color = "weighted"), size = 2) +
    scale_color_manual(values = c("youden" = "darkred", "weighted" = "darkblue")) +
    theme_classic()
ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/roc_mla_w.png", fig_roc_mla_w, width = 6, height = 4)

### proximal narrowing
data_mla <- baseline %>% select(ffr_0.8, ccta_pn_dist) %>% drop_na()
roc_mla_pn <- roc(data_mla$ffr_0.8, mdl_mla_pn_ffr$fitted.values)
cutoff_best <- coords(roc_mla_pn, x = "best", best.method = "youden", ret = "all")
cutoff_truepos <- coords(roc_mla_pn, x = "best", best.method = "youden", best.weights = c(1, 38/50), ret = "all") # with weights for high true positive

x_youden <- 1 - cutoff_best$specificity
y_youden <- cutoff_best$sensitivity
x_weighted <- 1 - cutoff_truepos$specificity
y_weighted <- cutoff_truepos$sensitivity

fig_roc_mla_pn <- ggroc(roc_mla_pn, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("MLA proximal narrowing") +
    geom_point(aes(x = x_youden, y = y_youden, color = "youden"), size = 2) +
    geom_point(aes(x = x_weighted, y = y_weighted, color = "weighted"), size = 2) +
    scale_color_manual(values = c("youden" = "darkred", "weighted" = "darkblue")) +
    theme_classic()
ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/roc_mla_pn.png", fig_roc_mla_pn, width = 6, height = 4)

# test different thresholds
list_prediction_data <- list("CCTA ostial area (BSA)" = all_ost_a_bsa[[2]], 
                            "CCTA ostial elliptic" = all_ost_ellip[[2]],
                            "CCTA ostial minor axis (BSA)" = all_ost_w_bsa[[2]],
                            "CCTA ostial proximal narrowing" = all_ost_pn[[2]],
                            "CCTA MLA (BSA)" = all_mla_a_bsa[[2]],
                            "CCTA MLA width (BSA)" = all_mla_w_bsa[[2]],
                            "CCTA MLA proximal narrowing" = all_mla_pn[[2]])
list_roc <- list("CCTA ostial area (BSA)" = roc_ost_a, 
                "CCTA ostial elliptic" = roc_ost_ellip,
                "CCTA ostial width (BSA)" = roc_ost_w,
                "CCTA ostial proximal narrowing" = roc_ost_pn,
                "CCTA MLA (BSA)" = roc_mla,
                "CCTA MLA width (BSA)" = roc_mla_w,
                "CCTA MLA proximal narrowing" = roc_mla_pn)

roc_table_creater <- function(list_roc, list_prediction_data, wgt) {
    roc_table <- data.frame()
    for (i in 1:length(list_roc)) {
        roc <- list_roc[[i]]
        prediction_data <- list_prediction_data[[i]]
        best_coords <- coords(roc, x = "best", best.method = "youden", ret = "all")
        wgt_coords <- coords(roc, x = "best", best.method = "youden", best.weights = wgt, ret = "all")
        dist_coords <- coords(roc, x = "best", best.method = "closest.topleft", ret = "all")

        value_cutoff <- prediction_data %>% filter(prediction_data[,2] > best_coords$threshold)
        value_cutoff_wgt <- prediction_data %>% filter(prediction_data[,2] > wgt_coords$threshold)
        value_cutoff_opt <- prediction_data %>% filter(prediction_data[,2] > dist_coords$threshold)
        
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
            "Weighted cutoff value" = round(max(value_cutoff_wgt[,1]), 2),
            "Distance threshold" = round(dist_coords$threshold, 3),
            "Distance accuracy" = round(dist_coords$accuracy*100, 0),
            "Distance sensitivity" = round(dist_coords$sensitivity*100, 0),
            "Distance specificity" = round(dist_coords$specificity*100, 0),
            "Distance PPV" = round(dist_coords$ppv*100, 0),
            "Distance NPV" = round(dist_coords$npv*100, 0),
            "Distance cutoff value" = round(max(value_cutoff_opt[,1]), 2)
        ))
    }
    return(roc_table)
}

roc_table <- roc_table_creater(list_roc, list_prediction_data, c(1, 38/50))
write_csv(roc_table, "/statistical_analysis/data/roc_table.csv")