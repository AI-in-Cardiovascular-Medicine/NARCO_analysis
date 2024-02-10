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
        ffr_0.8 = ifelse(ffr_0.8 == "yes", 1, 0),
        ffr_0.81 = ifelse(ffr_0.81 == "yes", 1, 0)
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

    tryCatch(
        {
            outcome <- table(explanatory, response)
            confusion <- conf_mat(outcome)
            plot_acc <- autoplot(confusion)
        },
        error = function(e) {
            print("Unfactor the explanatory variable first for confusion matrix and accuracy plot!")
        }
    )
    return(list(mdl, odds_ratio, lower_ci, upper_ci, confusion, plot_acc))
}

# regression models #################################################################################################
# exploration of continuous variables: ccta_slo_num, ccta_mla_elliptic, ccta_imc_length, ccta_pn_dist, ccta_aa_degree
# slit-like ostium combination of area and elliptic ratio with formula {slo_a / slo_elliptic}
all_slo <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_slo_num")
mdl_slo_ffr <- all_slo[[1]] # significant
summary(mdl_slo_ffr)
exp(coef(mdl_slo_ffr))
prediction_data_slo <- all_slo[[2]]
print(prediction_data_slo, n = 50)
plot_log_slo <- all_slo[[3]]
plot_log_slo 
plot_odds_slo <- all_slo[[4]]
plot_odds_slo
confusion_slo <- all_slo[[5]]
confusion_slo
summary(confusion_slo, event_level = "second")
plot_acc_slo <- all_slo[[6]]
plot_acc_slo

# test to split the variable again into ostial area and ostial ellipticity
data <- baseline %>% select(ccta_ostial_a, ccta_ostial_elliptic, ffr_0.8) %>% drop_na()
mdl_slo_multi <- glm(ffr_0.8 ~ ccta_ostial_a + ccta_ostial_elliptic, data = data, family = "binomial") # only ostial area is important
mdl_slo_multi_interaction <- glm(ffr_0.8 ~ ccta_ostial_a * ccta_ostial_elliptic, data = data, family = "binomial") # shows that interaction doesn't matter

all_ost_a <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_a")
mdl_ost_a_ffr <- all_ost_a[[1]]
summary(mdl_ost_a_ffr) # significant
prediction_data <- all[[2]]
print(prediction_data, n = 50)
plot_log_ost_a <- all_ost_a[[3]]
plot_log_ost_a
plot_odds_ost_a <- all_ost_a[[4]]
plot_odds_ost_a
confusion_ost_a <- all_ost_a[[5]]
confusion_ost_a
summary(confusion_ost_a, event_level = "second")
plot_acc_ost_a <- all_ost_a[[6]]
plot_acc_ost_a


# double check that elliptic is not significant
all_ost_ellip <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_elliptic")
mdl_ost_ellip_ffr <- all_ost_ellip[[1]]
summary(mdl_ost_ellip_ffr)
prediction_data_ost_ellip <- all_ost_ellip[[2]]
print(prediction_data_ost_ellip, n = 50)


# ostium proximal narrowing
all_ost_pn <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_pn")
mdl_ost_pn_ffr <- all_ost_pn[[1]]
summary(mdl_ost_pn_ffr) # not significant
prediction_data_ost_pn <- all_ost_pn[[2]]
print(prediction_data_ost_pn, n = 50)


# mla
all_mla <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_a")
mdl_mla_ffr <- all_mla[[1]]
summary(mdl_mla_ffr) # significant
prediction_data_mla <- all_mla[[2]]
print(prediction_data_mla, n = 50)
plot_log_mla <- all_mla[[3]]
plot_log_mla
plot_odds_mla <- all_mla[[4]]
plot_odds_mla
confusion_mla <- all_mla[[5]]
confusion_mla
summary(confusion_mla, event_level = "second")
plot_acc_mla <- all_mla[[6]]
plot_acc_mla


# exploration of ccta_mla_elliptic
all_ellip <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_elliptic")
mdl_ellip <- all_ellip[[1]] 
summary(mdl_ellip) # not significant
prediction_data_ellip <- all_ellip[[2]]
print(prediction_data_ellip, n = 50)
plot_log_ellip <- all_ellip[[3]]
plot_odds_ellip <- all_ellip[[4]]
plot_odds_ellip
confusion_ellip <- all_ellip[[5]]
summary(confusion, event_level = "second")
plot_acc <- all_ellip[[6]]


# exploration of ccta_pn_dist
all_pn <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_pn_dist")
mdl_pn_ffr <- all_pn[[1]]
summary(mdl_pn_ffr) # not significant
prediction_data_pn <- all_pn[[2]]
print(prediction_data_pn, n = 50)
plot_log_pn <- all_pn[[3]]
plot_log_pn
plot_odds_pn <- all_pn[[4]]
plot_odds_pn
confusion_pn <- all_pn[[5]]
summary(confusion_pn, event_level = "second")
plot_acc_pn <- all_pn[[6]]
plot_acc_pn


# exploration of ccta_imc_length
all_imc <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_imc_length")
mdl_imc_ffr <- all_imc[[1]]
summary(mdl_imc_ffr) # not significant
prediction_data_imc <- all_imc[[2]]
print(prediction_data_imc, n = 50)
plot_log_imc <- all_imc[[3]]
plot_log_imc
plot_odds_imc <- all_imc[[4]]
plot_odds_imc
confusion_imc <- all_imc[[5]]
confusion_imc
summary(confusion_imc, event_level = "second")
plot_acc_imc <- all_imc[[6]]
plot_acc_imc


# exploration of ccta_aa_degree
all_aa <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_aa_degree")
mdl_aa_ffr <- all_aa[[1]]
summary(mdl_aa_ffr) # not significant
prediction_data_aa <- all_aa[[2]]
print(prediction_data_aa, n = 50)
plot_log_aa <- all_aa[[3]]
plot_log_aa
plot_odds_aa <- all_aa[[4]]
plot_odds_aa
confusion_aa <- all_aa[[5]]
summary(confusion_aa, event_level = "second")
plot_acc_aa <- all_aa[[6]]
plot_acc_aa

# forest plot
# looping through models is pain!
forest_df <- data.frame(
    var = c("ostial_area", "ostial_elliptic", "ostial_pn", "mla", "mla_elliptic", "pn_dist", "imc_length", "aa_degree"),
    odds_ratio = c(exp(mdl_ost_a_ffr$coefficients[2]), exp(mdl_ost_ellip_ffr$coefficients[2]), exp(mdl_ost_pn_ffr$coefficients[2]), exp(mdl_mla_ffr$coefficients[2]), exp(mdl_ellip$coefficients[2]), exp(mdl_pn_ffr$coefficients[2]), exp(mdl_imc_ffr$coefficients[2]), exp(mdl_aa_ffr$coefficients[2])),
    lower_ci = c(exp(confint(mdl_ost_a_ffr))[2], exp(confint(mdl_ost_ellip_ffr))[2], exp(confint(mdl_ost_pn_ffr))[2], exp(confint(mdl_mla_ffr))[2], exp(confint(mdl_ellip))[2], exp(confint(mdl_pn_ffr))[2], exp(confint(mdl_imc_ffr))[2], exp(confint(mdl_aa_ffr))[2]),
    upper_ci = c(exp(confint(mdl_ost_a_ffr))[4], exp(confint(mdl_ost_ellip_ffr))[4], exp(confint(mdl_ost_pn_ffr))[4], exp(confint(mdl_mla_ffr))[4], exp(confint(mdl_ellip))[4], exp(confint(mdl_pn_ffr))[4], exp(confint(mdl_imc_ffr))[4], exp(confint(mdl_aa_ffr))[4])
)
ggplot(data = forest_df, aes(y = var, x = odds_ratio, xmin = lower_ci, xmax = upper_ci)) +
    geom_point() +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    theme_minimal()

# quick check how much model would change with cut-off 0.81
mdl_ost_a_ffr_081 <- simple_logistic_regression(baseline, "ffr_0.81", "ccta_ostial_a", family = "binomial")
mdl_ost_ellip_ffr_081 <- simple_logistic_regression(baseline, "ffr_0.81", "ccta_ostial_elliptic", family = "binomial")
mdl_ost_pn_ffr_081 <- simple_logistic_regression(baseline, "ffr_0.81", "ccta_ostial_pn", family = "binomial")
mdl_mla_ffr_081 <- simple_logistic_regression(baseline, "ffr_0.81", "ccta_mla_a", family = "binomial")
mdl_ellip_081 <- simple_logistic_regression(baseline, "ffr_0.81", "ccta_mla_elliptic", family = "binomial")
mdl_pn_ffr_081 <- simple_logistic_regression(baseline, "ffr_0.81", "ccta_pn_dist", family = "binomial")
mdl_imc_ffr_081 <- simple_logistic_regression(baseline, "ffr_0.81", "ccta_imc_length", family = "binomial")
mdl_aa_ffr_081 <- simple_logistic_regression(baseline, "ffr_0.81", "ccta_aa_degree", family = "binomial")
forest_df0.81 <- data.frame(
    var = c("ostial_area", "ostial_elliptic", "ostial_pn", "mla", "mla_elliptic", "pn_dist", "imc_length", "aa_degree"),
    odds_ratio = c(exp(mdl_ost_a_ffr_081[[1]]$coefficients[2]), exp(mdl_ost_ellip_ffr_081[[1]]$coefficients[2]), exp(mdl_ost_pn_ffr_081[[1]]$coefficients[2]), exp(mdl_mla_ffr_081[[1]]$coefficients[2]), exp(mdl_ellip_081[[1]]$coefficients[2]), exp(mdl_pn_ffr_081[[1]]$coefficients[2]), exp(mdl_imc_ffr_081[[1]]$coefficients[2]), exp(mdl_aa_ffr_081[[1]]$coefficients[2])),
    lower_ci = c(exp(confint(mdl_ost_a_ffr_081[[1]]))[2], exp(confint(mdl_ost_ellip_ffr_081[[1]]))[2], exp(confint(mdl_ost_pn_ffr_081[[1]]))[2], exp(confint(mdl_mla_ffr_081[[1]]))[2], exp(confint(mdl_ellip_081[[1]]))[2], exp(confint(mdl_pn_ffr_081[[1]]))[2], exp(confint(mdl_imc_ffr_081[[1]]))[2], exp(confint(mdl_aa_ffr_081[[1]]))[2]),
    upper_ci = c(exp(confint(mdl_ost_a_ffr_081[[1]]))[4], exp(confint(mdl_ost_ellip_ffr_081[[1]]))[4], exp(confint(mdl_ost_pn_ffr_081[[1]]))[4], exp(confint(mdl_mla_ffr_081[[1]]))[4], exp(confint(mdl_ellip_081[[1]]))[4], exp(confint(mdl_pn_ffr_081[[1]]))[4], exp(confint(mdl_imc_ffr_081[[1]]))[4], exp(confint(mdl_aa_ffr_081[[1]]))[4])
)
ggplot(data = forest_df0.81, aes(y = var, x = odds_ratio, xmin = lower_ci, xmax = upper_ci)) +
    geom_point() +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    theme_minimal()


####################################################################################################################
# binary explanatory variable
categorical <- baseline %>%
    mutate(
        caa_slo = ifelse(caa_slo == "yes", 1, 0),
        caa_elliptic = ifelse(caa_elliptic == "yes", 1, 0),
        caa_pn = ifelse(caa_pn == "yes", 1, 0),
        caa_im = ifelse(caa_im == "yes", 1, 0),
        caa_angle = ifelse(caa_angle == "yes", 1, 0)
    )

# slit like ostium
all <- logistic_regression_categorical(categorical, "ffr_0.8", "caa_slo")
mdl <- all[[1]]
odds_ratio <- all[[2]]
lower_ci <- all[[3]]
upper_ci <- all[[4]]
summary(mdl) # significant

# elliptic vessel shape
all <- logistic_regression_categorical(categorical, "ffr_0.8", "caa_elliptic")
mdl <- all[[1]]
summary(mdl) # not significant

# proximal narrowing
all <- logistic_regression_categorical(categorical, "ffr_0.8", "caa_pn")
mdl <- all[[1]]
summary(mdl) # not significant

# intramural course
all <- logistic_regression_categorical(categorical, "ffr_0.8", "caa_im")
mdl <- all[[1]]
summary(mdl) # not significant

# acute angle
all <- logistic_regression_categorical(categorical, "ffr_0.8", "caa_angle")
mdl <- all[[1]]
summary(mdl) # not significant

# forest plots
vars = c("caa_slo", "caa_pn")
forest_df <- data.frame(var = character(), odds_ratio = numeric(), lower_ci = numeric(), upper_ci = numeric(), stringsAsFactors = FALSE)
for (var in vars){
    all <- logistic_regression_categorical(categorical, "ffr_0.8", var)
    mdl <- all[[1]]
    odds_ratio <- all[[2]]
    lower_ci <- all[[3]]
    upper_ci <- all[[4]]
    forest_df <- rbind(forest_df, data.frame(var = var, odds_ratio = odds_ratio, lower_ci = lower_ci, upper_ci = upper_ci))
}

ggplot(data = forest_df, aes(y = var, x = odds_ratio, xmin = lower_ci, xmax = upper_ci)) +
    geom_point() +
    geom_errorbarh(height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    theme_minimal()