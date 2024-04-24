library(tidyverse)
library(broom)
library(ggplot2)
library(yardstick)
library(pROC)
library(randomForest)
library(car)
library(readxl)
library(lmtest)
library(corrplot)

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
        # ccta_ostial_pn = (1 - (ccta_ostial_a / ccta_dist_a)) * 100,
        # ccta_pn_dist = (1 - (ccta_mla_a / ccta_dist_a)) * 100,
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
        ccta_pn_dist = ifelse(ccta_pn_dist < 0, 0, ccta_pn_dist),
        ccta_ostial_pn = ifelse(ccta_ostial_pn < 0, 0, ccta_ostial_pn)
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

# all continuous models
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

# all categorical models
# ostium
mdl_slo_ffr <- glm(ffr_0.8 ~ caa_slo, data = baseline, family = "binomial")
summary(mdl_slo_ffr) # not significant!!! AIC 49.59
exp(coef(mdl_slo_ffr))[2]
exp(confint(mdl_slo_ffr))[2]
exp(confint(mdl_slo_ffr))[4]
lrtest(mdl_slo_ffr, mdl_ost_a_bsa_ffr) # mdl_ost_a_bsa_ffr is better
lrtest(mdl_slo_ffr, mdl_ost_ellip_ffr) # mdl_ost_ellip_ffr is better
mdl_ellip_ffr <- glm(ffr_0.8 ~ caa_elliptic, data = baseline, family = "binomial")
summary(mdl_ellip_ffr) # not significant!!!
mdl_pn_ffr <- glm(ffr_0.8 ~ caa_pn, data = baseline, family = "binomial")
summary(mdl_pn_ffr) # significant!!! AIC 53.184
summary(mdl_ost_pn_ffr) # significant!!! AIC 47.518
lrtest(mdl_pn_ffr, mdl_ost_pn_ffr) # mdl_ost_pn_ffr is better
# mdl_aa_ffr <- glm(ffr_0.8 ~ caa_angle, data = baseline, family = "binomial")
# summary(mdl_aa_ffr) # not significant!!! all acute take-off
# mdl_im_ffr <- glm(ffr_0.8 ~ caa_im, data = baseline, family = "binomial")
# summary(mdl_im_ffr) # not significant!!! all intramural


# forest plot
# bsa adjusted
forest_df_bsa <- data.frame(
    var = c("Ostial area (BSA) [mm2]", "Ostial elliptic ratio", "Ostial major axis (BSA) [mm]", "Ostial minor axis (BSA) [mm]", 
            "Ostial proximal narrowing [%]", "IM area (BSA) [mm2]", "IM elliptic ratio", "IM major axis (BSA) [mm]", "IM minor axis (BSA) [mm]", 
            "IM proximal narrowing [%]", "Intramural length (BSA) [mm]", "Take-off angle (BSA) [deg]"),
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

forest_df_bsa$var <- factor(forest_df_bsa$var, levels = rev(forest_df_bsa$var))

forest_cut <- forest_df_bsa %>%
    mutate(
        upper_ci = ifelse(upper_ci > 6, 6, upper_ci)
    )

forest_cut_ffr_bsa <- ggplot(data = forest_cut, aes(y = var, x = odds_ratio, xmin = lower_ci, xmax = upper_ci)) +
    geom_point() +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, 7, 1), labels = parse(text = seq(0, 7, 1))) +
    geom_vline(xintercept = 6, linetype = "solid", color = "white", lwd = 2) +
    geom_segment(aes(x = 6, xend = 7, y = 6, yend = 6), linetype = "dashed") +
    geom_segment(aes(x = 6, xend = 7, y = 11, yend = 11), linetype = "dashed") +
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
## ostium ###########################################################################################################
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
x_weighted_ost_a <- 1 - cutoff_truepos$specificity
y_weighted_ost_a <- cutoff_truepos$sensitivity
x_min_dist <- 1 - roc_ost_a$specificities[optimal_index]
y_min_dist <- roc_ost_a$sensitivities[optimal_index]
fig_roc_ost_a <- ggroc(roc_ost_a, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Ostial area") +
    geom_point(aes(x = x_youden, y = y_youden, color = "youden"), size = 2) +
    geom_point(aes(x = x_weighted_ost_a, y = y_weighted_ost_a, color = "weighted"), size = 2) +
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
x_weighted_ost_ellip <- 1 - cutoff_truepos$specificity
y_weighted_ost_ellip <- cutoff_truepos$sensitivity

fig_roc_ost_ellip <- ggroc(roc_ost_ellip, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Ostial elliptic") +
    geom_point(aes(x = x_youden, y = y_youden, color = "youden"), size = 2) +
    geom_point(aes(x = x_weighted_ost_ellip, y = y_weighted_ost_ellip, color = "weighted"), size = 2) +
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
x_weighted_ost_w <- 1 - cutoff_truepos$specificity
y_weighted_ost_w <- cutoff_truepos$sensitivity

fig_roc_ost_w <- ggroc(roc_ost_w, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Ostial minor axis") +
    geom_point(aes(x = x_youden, y = y_youden, color = "youden"), size = 2) +
    geom_point(aes(x = x_weighted_ost_w, y = y_weighted_ost_w, color = "weighted"), size = 2) +
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
x_weighted_ost_pn <- 1 - cutoff_truepos$specificity
y_weighted_ost_pn <- cutoff_truepos$sensitivity

fig_roc_ost_pn <- ggroc(roc_ost_pn, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Ostial proximal narrowing") +
    geom_point(aes(x = x_youden, y = y_youden, color = "youden"), size = 2) +
    geom_point(aes(x = x_weighted_ost_pn, y = y_weighted_ost_pn, color = "weighted"), size = 2) +
    scale_color_manual(values = c("youden" = "darkred", "weighted" = "darkblue")) +
    theme_classic()
ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/roc_ost_pn.png", fig_roc_ost_pn, width = 6, height = 4)

## mla #############################################################################################################
data_mla <- baseline %>% select(ffr_0.8, ccta_mla_a_bsa) %>% drop_na()
### area
roc_mla <- roc(data_mla$ffr_0.8, mdl_mla_a_bsa_ffr$fitted.values)
cutoff_best <- coords(roc_mla, x = "best", best.method = "youden", ret = "all")
cutoff_truepos <- coords(roc_mla, x = "best", best.method = "youden", best.weights = c(1, 38/50), ret = "all") # with weights for high true positive

x_youden <- 1 - cutoff_best$specificity
y_youden <- cutoff_best$sensitivity
x_weighted_mla <- 1 - 0.36842105
y_weighted_mla <- 1
threshold_mla <- 0.088438793

fig_roc_mla <- ggroc(roc_mla, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("Minimal lumen area") +
    geom_point(aes(x = x_youden, y = y_youden, color = "youden"), size = 2) +
    geom_point(aes(x = x_weighted_mla, y = y_weighted_mla, color = "weighted"), size = 2) +
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
x_weighted_mla_w <- 1 - 0.13157895
y_weighted_mla_w <- 1
threshold_mla_w <- 0.039378259

fig_roc_mla_w <- ggroc(roc_mla_w, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("MLA minor axis") +
    geom_point(aes(x = x_youden, y = y_youden, color = "youden"), size = 2) +
    geom_point(aes(x = x_weighted_mla_w, y = y_weighted_mla_w, color = "weighted"), size = 2) +
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
x_weighted_mla_pn <- 1 - 0.13157895
y_weighted_mla_pn <- 1
threshold_mla_pn <- 0.04657369

fig_roc_mla_pn <- ggroc(roc_mla_pn, legacy.axes = TRUE) +
    scale_color_manual(values = c("black", "blue")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    ggtitle("MLA proximal narrowing") +
    geom_point(aes(x = x_youden, y = y_youden, color = "youden"), size = 2) +
    geom_point(aes(x = x_weighted_mla_pn, y = y_weighted_mla_pn, color = "weighted"), size = 2) +
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


roc_table <- roc_table_creater(list_roc, list_prediction_data, c(1, 38/50))
write_csv(roc_table, "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/roc_table.csv")

# combined ROC plots
list_roc_ost <- list("CCTA ostial elliptic" = roc_ost_ellip,
                    "CCTA ostial proximal narrowing" = roc_ost_pn,
                    "CCTA ostial area (BSA)" = roc_ost_a, 
                    "CCTA ostial minor axis (BSA)" = roc_ost_w)

roc_ost_combined <- ggroc(list_roc_ost, legacy.axes = TRUE) +
    scale_color_manual(values = c("#1f1fb8", "#1c771c", "#fa9600", "#a31800", "black")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_point(aes(x = x_weighted_ost_a, y = y_weighted_ost_a, color = "cut-off"), size = 2) +
    geom_point(aes(x = x_weighted_ost_ellip, y = y_weighted_ost_ellip, color = "cut-off"), size = 2) +
    geom_point(aes(x = x_weighted_ost_w, y = y_weighted_ost_w, color = "cut-off"), size = 2) +
    geom_point(aes(x = x_weighted_ost_pn, y = y_weighted_ost_pn, color = "cut-off"), size = 2) +
    ggtitle("Ostial measurements") +
    theme_classic() +
    # add all AUCs with color
    annotate("text", x = 0.5, y = 0.1, label = paste0("AUC elliptic ratio: ", round(auc(roc_ost_ellip), 3)), size = 3, hjust = 0) +
    annotate("text", x = 0.5, y = 0.2, label = paste0("AUC prox. narrowing: ", round(auc(roc_ost_pn), 3)), size = 3, hjust = 0) +
    annotate("text", x = 0.5, y = 0.3, label = paste0("AUC area: ", round(auc(roc_ost_a), 3)), size = 3, hjust = 0) +
    annotate("text", x = 0.5, y = 0.4, label = paste0("AUC minor axis: ", round(auc(roc_ost_w), 3)), size = 3, hjust = 0)

ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/roc_ost_combined.png", roc_ost_combined, width = 8, height = 6)

# combined ROC plots
list_roc_ost <- list(
                    "CCTA intramural area (BSA)" = roc_mla, 
                    "CCTA intramural proximal narrowing" = roc_mla_pn,
                    "CCTA intramural minor axis (BSA)" = roc_mla_w)

roc_mla_combined <- ggroc(list_roc_ost, legacy.axes = TRUE) +
    scale_color_manual(values = c("#1f1fb8", "#fa9600", "#a31800", "black")) +
    theme(legend.position = "none") +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_point(aes(x = x_weighted_mla, y = y_weighted_mla, color = "cut-off"), size = 2) +
    geom_point(aes(x = x_weighted_mla_w, y = y_weighted_mla_w, color = "cut-off"), size = 2) +
    geom_point(aes(x = x_weighted_mla_pn, y = y_weighted_mla_pn, color = "cut-off"), size = 2) +
    ggtitle("Intramural measurements") +
    theme_classic() +
    # add all AUCs with color
    annotate("text", x = 0.5, y = 0.2, label = paste0("AUC prox. narrowing: ", round(auc(roc_mla_pn), 3)), size = 3, hjust = 0) +
    annotate("text", x = 0.5, y = 0.3, label = paste0("AUC area: ", round(auc(roc_mla), 3)), size = 3, hjust = 0) +
    annotate("text", x = 0.5, y = 0.4, label = paste0("AUC minor axis: ", round(auc(roc_mla_w), 3)), size = 3, hjust = 0)

ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/roc_mla_combined.png", roc_mla_combined, width = 8, height = 6)

# correlation matrix of all high risk features
high_risk_features <- baseline %>% select(ccta_ostial_a_bsa, ccta_mla_a_bsa, ccta_ostial_elliptic, ccta_mla_elliptic, ccta_ostial_h_bsa, ccta_mla_h_bsa, 
                                        ccta_ostial_w_bsa, ccta_mla_w_bsa, ccta_ostial_pn, ccta_pn_dist, 
                                        ccta_imc_length_bsa, ccta_aa_degree_bsa)

colnames(high_risk_features) <- c("OA (BSA)", "IM-LA (BSA)", "OA Elliptic", "IM-LA Elliptic", "OA major (BSA)", "IM-LA major (BSA)", 
                                  "OA minor (BSA)", "IM-LA minor (BSA)", "Ostial LN", "IM LN", "IM length (BSA)", "Take-off angle (BSA)")

correlation_matrix <- cor(high_risk_features, use = "pairwise.complete.obs")
correlation_matrix <- round(correlation_matrix, 2)
# save as excel
write.csv(correlation_matrix, "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/correlation_matrix.csv")

# plot correlation matrix
png(height = 800, width = 800, file = "C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/correlation_matrix.png")
fig_correlation_matrix <- 
    corrplot(correlation_matrix, 
            method = "color", 
            type = "upper", 
            tl.col = "black", 
            tl.srt = 45, 
            tl.cex = 1.1,
            addCoef.col = "black", 
            number.digits = 2)

dev.off()


# find exact value for cut-off
explanatory_data <- tibble(
    ccta_pn_dist = seq(23, 26, 1)
)

prediction_data <- explanatory_data %>% 
    mutate(
        response = predict(mdl_mla_pn_ffr, explanatory_data, type = "response"),
        most_likely_outcome = round(response),
        odds_ratio = response / (1 - response),
        log_odds_ratio = log(odds_ratio),
        log_odds_ratio2 = predict(mdl_mla_pn_ffr, explanatory_data)
    )

library(yardstick)

table <- baseline %>% mutate(ccta_ostial_a_bsa = ifelse(ccta_ostial_a_bsa >4.13, 0, 1)) %>% select(ccta_ostial_a_bsa, ffr_0.8) %>% drop_na()
table <- baseline %>% mutate(ccta_ostial_elliptic = ifelse(ccta_ostial_elliptic <1.86, 0, 1)) %>% select(ccta_ostial_elliptic, ffr_0.8) %>% drop_na()
table <- baseline %>% mutate(ccta_ostial_w_bsa = ifelse(ccta_ostial_w_bsa >0.94, 0, 1)) %>% select(ccta_ostial_w_bsa, ffr_0.8) %>% drop_na()
table <- baseline %>% mutate(ccta_ostial_pn = ifelse(ccta_ostial_pn <26, 0, 1)) %>% select(ccta_ostial_pn, ffr_0.8) %>% drop_na()
table <- baseline %>% mutate(ccta_mla_a_bsa = ifelse(ccta_mla_a_bsa >3.3, 0, 1)) %>% select(ccta_mla_a_bsa, ffr_0.8) %>% drop_na()
table <- baseline %>% mutate(ccta_mla_w_bsa = ifelse(ccta_mla_w_bsa >1.1, 0, 1)) %>% select(ccta_mla_w_bsa, ffr_0.8) %>% drop_na()
table <- baseline %>% mutate(ccta_pn_dist = ifelse(ccta_pn_dist <24, 0, 1)) %>% select(ccta_pn_dist, ffr_0.8) %>% drop_na()

outcome <- table(table)
confusion <- conf_mat(outcome)
summary(confusion, event_level = "second")

# repeat whole forest plot as normalized
# normalize variables in baseline
normalized <- baseline %>% mutate(
    ccta_ostial_a_bsa = -1 *(ccta_ostial_a_bsa - mean(ccta_ostial_a_bsa, na.rm = TRUE)) / sd(ccta_ostial_a_bsa, na.rm = TRUE),
    ccta_ostial_elliptic = (ccta_ostial_elliptic - mean(ccta_ostial_elliptic, na.rm = TRUE)) / sd(ccta_ostial_elliptic, na.rm = TRUE),
    ccta_ostial_w_bsa = -1 * (ccta_ostial_w_bsa - mean(ccta_ostial_w_bsa, na.rm = TRUE)) / sd(ccta_ostial_w_bsa, na.rm = TRUE),
    ccta_ostial_h_bsa = -1 * (ccta_ostial_h_bsa - mean(ccta_ostial_h_bsa, na.rm = TRUE)) / sd(ccta_ostial_h_bsa, na.rm = TRUE),
    ccta_ostial_pn = (ccta_ostial_pn - mean(ccta_ostial_pn, na.rm = TRUE)) / sd(ccta_ostial_pn, na.rm = TRUE),
    ccta_mla_a_bsa = -1 * (ccta_mla_a_bsa - mean(ccta_mla_a_bsa, na.rm = TRUE)) / sd(ccta_mla_a_bsa, na.rm = TRUE),
    ccta_mla_elliptic = (ccta_mla_elliptic - mean(ccta_mla_elliptic, na.rm = TRUE)) / sd(ccta_mla_elliptic, na.rm = TRUE),
    ccta_mla_w_bsa = -1 * (ccta_mla_w_bsa - mean(ccta_mla_w_bsa, na.rm = TRUE)) / sd(ccta_mla_w_bsa, na.rm = TRUE),
    ccta_mla_h_bsa = -1 * (ccta_mla_h_bsa - mean(ccta_mla_h_bsa, na.rm = TRUE)) / sd(ccta_mla_h_bsa, na.rm = TRUE),
    ccta_pn_dist = (ccta_pn_dist - mean(ccta_pn_dist, na.rm = TRUE)) / sd(ccta_pn_dist, na.rm = TRUE),
    ccta_imc_length_bsa = (ccta_imc_length_bsa - mean(ccta_imc_length_bsa, na.rm = TRUE)) / sd(ccta_imc_length_bsa, na.rm = TRUE),
    ccta_aa_degree_bsa = (ccta_aa_degree_bsa - mean(ccta_aa_degree_bsa, na.rm = TRUE)) / sd(ccta_aa_degree_bsa, na.rm = TRUE)
)


# ostium
all_ost_a_bsa <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_ostial_a_bsa")
mdl_ost_a_bsa_ffr <- all_ost_a_bsa[[1]] # significant
all_ost_ellip <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_ostial_elliptic")
mdl_ost_ellip_ffr <- all_ost_ellip[[1]] # significant
all_ost_h_bsa <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_ostial_h_bsa")
mdl_ost_h_bsa_ffr <- all_ost_h_bsa[[1]] # not significant!!!!
all_ost_w_bsa <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_ostial_w_bsa")
mdl_ost_w_bsa_ffr <- all_ost_w_bsa[[1]] # significant
all_ost_pn <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_ostial_pn")
mdl_ost_pn_ffr <- all_ost_pn[[1]] # significant

# mla
all_mla_a_bsa <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_mla_a_bsa")
mdl_mla_a_bsa_ffr <- all_mla_a_bsa[[1]] # significant
all_mla_ellip <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_mla_elliptic")
mdl_mla_ellip_ffr <- all_mla_ellip[[1]] # not significant!!!
all_mla_h_bsa <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_mla_h_bsa")
mdl_mla_h_bsa_ffr <- all_mla_h_bsa[[1]] # not significant!!!
all_mla_w_bsa <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_mla_w_bsa")
mdl_mla_w_bsa_ffr <- all_mla_w_bsa[[1]] # significant
all_mla_pn <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_pn_dist")
mdl_mla_pn_ffr <- all_mla_pn[[1]] # significant

# rest of the variables
all_imc_bsa <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_imc_length_bsa")
mdl_imc_bsa_ffr <- all_imc_bsa[[1]] # not significant!!!
all_aa_bsa <- simple_logistic_regression(normalized, "ffr_0.8", "ccta_aa_degree_bsa")
mdl_aa_bsa_ffr <- all_aa_bsa[[1]] # not significant!!!


forest_df_bsa <- data.frame(
    var = c("Ostial area (BSA) [mm2]", "Ostial elliptic ratio", "Ostial major axis (BSA) [mm]", "Ostial minor axis (BSA) [mm]", 
            "Ostial proximal narrowing [%]", "IM area (BSA) [mm2]", "IM elliptic ratio", "IM major axis (BSA) [mm]", "IM minor axis (BSA) [mm]", 
            "IM proximal narrowing [%]", "Intramural length (BSA) [mm]", "Take-off angle (BSA) [deg]"),
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

forest_df_bsa$var <- factor(forest_df_bsa$var, levels = rev(forest_df_bsa$var))

forest_cut <- forest_df_bsa %>%
    mutate(
        upper_ci = ifelse(upper_ci > 15, 15, upper_ci)
    )

forest_cut_ffr_bsa <- ggplot(data = forest_cut, aes(y = var, x = odds_ratio, xmin = lower_ci, xmax = upper_ci)) +
    geom_point() +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    scale_x_continuous(breaks = seq(0, 17, 1), labels = parse(text = seq(0, 17, 1))) +
    geom_vline(xintercept = 15, linetype = "solid", color = "white", lwd = 2) +
    theme_classic()

#completely new plot
# for ccta_ostial_a_bsa get min and max values, then create a seismic scale that has it's white bar at 4.13
# Get min and max values
min_val <- min(baseline$ccta_ostial_w_bsa, na.rm = TRUE)
max_val <- max(baseline$ccta_ostial_w_bsa, na.rm = TRUE)

# Create a data frame for the rectangle
rect_df <- data.frame(xmin = min_val, xmax = max_val, ymin = 0, ymax = 1)

# Create the plot
ggplot() +
    geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), 
                        fill = "white", color = "black") +
    geom_vline(xintercept = 0.94, linetype = "solid", color = "black") +
    coord_fixed(ratio = 1) +
    scale_x_continuous(breaks = seq(0, 2, 0.2), labels = parse(text = seq(0, 2, 0.2))) +
    theme_minimal()

####################################################################################################################
# distribution with different cut-offs
test_set <- list(coords(roc_ost_a), coords(roc_ost_ellip), 
            coords(roc_ost_w), coords(roc_ost_pn), 
            coords(roc_mla), coords(roc_mla_w), coords(roc_mla_pn))

prediction_data <- list(prediction_data_ost_a, prediction_data_ost_ellip, 
                        prediction_data_ost_w, prediction_data_ost_pn, 
                        prediction_data_mla, prediction_data_mla_w, prediction_data_mla_pn)

thresholds <- list()

for (i in 1:length(test_set)){
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

test <- thresholds[[1]][1]
test_prediction <- prediction_data[[1]]

test_prediction$ccta_ostial_a_bsa[which.min(abs(test_prediction$ffr_0.8 - test))]

mdl <- glm(ffr_0.8 ~ ccta_ostial_a_bsa, data = baseline, family = binomial)

explanatory_data <- tibble(
    ccta_ostial_a_bsa = seq(min(baseline$ccta_ostial_a_bsa, na.rm = TRUE), max(baseline$ccta_ostial_a_bsa, na.rm = TRUE), length.out = 1000)
)

prediction_data <- explanatory_data %>% 
    mutate(
        response = predict(mdl, explanatory_data, type = "response"),
        most_likely_outcome = round(response),
        odds_ratio = response / (1 - response),
        log_odds_ratio = log(odds_ratio),
        log_odds_ratio2 = predict(mdl, explanatory_data)
    )

value <- prediction_data$ccta_ostial_a_bsa[which.min(abs(prediction_data$response - test))]

table <- baseline %>% mutate(ccta_ostial_a_bsa = ifelse(ccta_ostial_a_bsa > value, 0, 1)) %>% select(ccta_ostial_a_bsa, ffr_0.8) %>% drop_na() %>% table()
conf_mat(table)
summary(confusion, event_level = "second")

table <- baseline %>% mutate(ccta_ostial_elliptic = ifelse(ccta_ostial_elliptic <1.86, 0, 1)) %>% select(ccta_ostial_elliptic, ffr_0.8) %>% drop_na()
table <- baseline %>% mutate(ccta_ostial_w_bsa = ifelse(ccta_ostial_w_bsa >0.94, 0, 1)) %>% select(ccta_ostial_w_bsa, ffr_0.8) %>% drop_na()
table <- baseline %>% mutate(ccta_ostial_pn = ifelse(ccta_ostial_pn <26, 0, 1)) %>% select(ccta_ostial_pn, ffr_0.8) %>% drop_na()
table <- baseline %>% mutate(ccta_mla_a_bsa = ifelse(ccta_mla_a_bsa >3.3, 0, 1)) %>% select(ccta_mla_a_bsa, ffr_0.8) %>% drop_na()
table <- baseline %>% mutate(ccta_mla_w_bsa = ifelse(ccta_mla_w_bsa >1.1, 0, 1)) %>% select(ccta_mla_w_bsa, ffr_0.8) %>% drop_na()
table <- baseline %>% mutate(ccta_pn_dist = ifelse(ccta_pn_dist <24, 0, 1)) %>% select(ccta_pn_dist, ffr_0.8) %>% drop_na()