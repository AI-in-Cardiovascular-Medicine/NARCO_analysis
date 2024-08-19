library(survival)
library(survminer)
library(ggplot2)
library(tidyverse)
library(readxl)

followup <- read_excel("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/dataframes/follow_up.xlsx")
followup$ae_mace_type <- factor(followup$ae_mace_type, levels = c(0, 1, 2, 3, 4), labels = c("Death", "VT", "MI", "Emergency revascularization", "ICD"))
baseline <- readRDS("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/data/baseline.rds")
index <- baseline %>% filter(!is.na(inv_ivusrest_mla)) %>% select(record_id, patient_id, ffr_0.8)
patho <- index %>% filter(ffr_0.8 == "yes")
normal <- index %>% filter(ffr_0.8 == "no")
followup <- followup %>% filter(record_id %in% index$record_id)
followup$ffr_0.8 <- index$ffr_0.8
followup$timeto_censor_y
followup$mace

median <- median(followup$timeto_censor_y, na.rm = TRUE)
quantile <- quantile(followup$timeto_censor_y, probs = c(0.25, 0.75), na.rm = TRUE)

surv_fit <- survfit(Surv(timeto_censor_y, mace) ~ ffr_0.8, data = followup)
gg <- ggsurvplot(surv_fit, pval = TRUE, conf.int = FALSE, risk.table = TRUE, 
                   legend.title = "FFR < 0.8", legend.labs = c("No", "Yes"), 
                   palette = "jco", xlab = "Time (years)", ylab = "Survival probability", 
                   legend = "right", xlim = c(0,5))

# Save plot
# ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/figures/survival_curve.jpg", plot = print(gg), width = 3.3, height = 2.2, dpi = 1000 )
ggsave("C:/WorkingData/Documents/2_Coding/Python/NARCO_analysis/statistical_analysis/figures/survival_curve.jpg", 
       plot = gg$plot, width = 6, height = 4, dpi = 1000)

print(paste0("Median survival time: ", median, " (", quantile[1], "-", quantile[2], ") years"))

followup %>%
    select(record_id, timeto_censor_y, mace, ffr_0.8, ae_mace_type, ae_mace_type_2, ae_mace_type_3, ae_mace_type_4) %>% print(n = 50)

# test of proportion for 1/38 and 2/12
prop.test(c(1, 2), c(38, 12), correct = FALSE)

followup <- followup %>%
  mutate(
    fu_symptoms = if_any(
      # Select all columns that start with "sym_fu___" and exclude "sym_fu___0" variants
      select(., starts_with("sym_fu___")) %>% select(-sym_fu___0, -sym_fu___0_y3, -sym_fu___0_y5) %>% colnames(),
      ~ . == 1
    ) * 1 # Convert logical TRUE/FALSE to 1/0
  )

followup <- followup %>%
  mutate(
    fu_nosymptoms = if_any(
      starts_with("sym_fu___0"),
      ~ . == 1
    ) * 1 # Convert logical TRUE/FALSE to 1/0
  )

followup <- followup %>%
    mutate(
        fu_any_symptoms = ifelse(fu_symptoms == 1, 1, ifelse(fu_nosymptoms == 1, 0, NA))
    )

baseline <- baseline %>% filter(record_id %in% index$record_id)

# calculate date diff between baseline$inv_date and followup$pf_date_fu
baseline$inv_date <- as.Date(baseline$inv_date, format = "%Y-%m-%d")
followup$pf_date_fu <- as.Date(followup$pf_date_fu, format = "%Y-%m-%d")
baseline$followup_date_diff <- as.numeric(difftime(followup$pf_date_fu, baseline$inv_date, units = "days"))
baseline$followup_date_diff <- round(baseline$followup_date_diff / 365.25, 2)
median(baseline$followup_date_diff, na.rm = TRUE)
quantile(baseline$followup_date_diff, probs = c(0.25, 0.75), na.rm = TRUE)
max(baseline$followup_date_diff, na.rm = TRUE)
min(baseline$followup_date_diff, na.rm = TRUE)