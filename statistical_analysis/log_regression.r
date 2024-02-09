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
        scale_y_continuous(breaks = seq(0, 8, 1))

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
            confusion <- "No confusion matrix"
            plot_acc <- "No accuracy plot"
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

# exploration of ccta_slo_num, ccta_mla_elliptic, ccta_imc_length, ccta_pn_dist, ccta_aa_degree
all <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_slo_num")
mdl <- all[[1]] # significant
summary(mdl)
prediction_data <- all[[2]]
print(prediction_data, n = 50)
plot_log <- all[[3]]
plot_log
plot_odds <- all[[4]]
plot_odds
confusion <- all[[5]]
summary(confusion, event_level = "second")
plot_acc <- all[[6]]
plot_acc

# test to split the variable again into ostial area and ostial ellipticity
data <- baseline %>% select(ccta_ostial_a, ccta_ostial_elliptic, ffr_0.8) %>% drop_na()
mdl <- glm(ffr_0.8 ~ ccta_ostial_a + ccta_ostial_elliptic, data = data, family = "binomial") # only ostial area is important
mdl <- glm(ffr_0.8 ~ ccta_ostial_a * ccta_ostial_elliptic, data = data, family = "binomial") # shows that interaction doesn't matter
all <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_ostial_a")
mdl <- all[[1]]
summary(mdl) # significant
prediction_data <- all[[2]]
print(prediction_data, n = 50)
plot_log <- all[[3]]
plot_log
plot_odds <- all[[4]]
plot_odds
confusion <- all[[5]]
confusion
summary(confusion, event_level = "second")
plot_acc <- all[[6]]
plot_acc



# exploration of ccta_mla_elliptic
all <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_mla_elliptic")
mdl <- all[[1]] 
summary(mdl) # not significant
prediction_data <- all[[2]]
plot_log <- all[[3]]
plot_odds <- all[[4]]
confusion <- all[[5]]
summary(confusion, event_level = "second")
plot_acc <- all[[6]]


# exploration of ccta_imc_length
all <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_imc_length")
mdl <- all[[1]]
summary(mdl) # not significant
prediction_data <- all[[2]]
print(prediction_data, n = 50)
plot_log <- all[[3]]
plot_log
plot_odds <- all[[4]]
plot_odds
confusion <- all[[5]]
confusion
summary(confusion, event_level = "second")
plot_acc <- all[[6]]
plot_acc

# exploration of ccta_pn_dist
all <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_pn_dist")
mdl <- all[[1]]
summary(mdl) # not significant
prediction_data <- all[[2]]
print(prediction_data, n = 50)
plot_log <- all[[3]]
plot_log
plot_odds <- all[[4]]
plot_odds
confusion <- all[[5]]
summary(confusion, event_level = "second")
plot_acc <- all[[6]]
plot_acc

# exploration of ccta_aa_degree
all <- simple_logistic_regression(baseline, "ffr_0.8", "ccta_aa_degree")
mdl <- all[[1]]
summary(mdl) # not significant
prediction_data <- all[[2]]
print(prediction_data, n = 50)
plot_log <- all[[3]]
plot_log
plot_odds <- all[[4]]
plot_odds
confusion <- all[[5]]
summary(confusion, event_level = "second")
plot_acc <- all[[6]]
plot_acc

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

# forest plot
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