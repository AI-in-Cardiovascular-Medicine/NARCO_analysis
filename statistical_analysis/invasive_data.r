library(tidyverse)
library(ggplot2)
library(ggpubr)
library(yaml)

yaml <- yaml.load_file(file.path(getwd(), "config.yaml"))
output_dir <- yaml$demographics$output_dir_figures
# Load the data
baseline <- readRDS(paste0(yaml$demographics$output_dir_data,"/baseline.rds"))
baseline <- baseline %>% mutate(
    inv_aomean_change = inv_dobu_aomean - inv_rest_aomean,
    inv_aosys_change = inv_dobu_aosys - inv_rest_aosys,
    inv_aodia_change = inv_dobu_aodia - inv_rest_aodia,
    inv_hr_change = inv_dobu_hr - inv_rest_hr,
    inv_ffr_change = inv_ffrado - inv_ffrdobu,
    inv_ivus_change = inv_ivusrest_mla - inv_ivusdobu_mla
)

# prep data for plotting ffr
n <- nrow(baseline)
inv_ffrado <- pull(baseline, inv_ffrado)
inv_ffrdobu <- pull(baseline, inv_ffrdobu)
method_1 <- rep("FFR adenosine", n)
method_2 <- rep("FFR dobutamine", n)
ffr_df <- data.frame(ffr_method = c(method_1, method_2), 
                    ffr_value = c(inv_ffrado, inv_ffrdobu))
colnames(ffr_df) <- c("ffr_method", "ffr_value")

plot_ffr_change <- ggpaired(ffr_df, x = "ffr_method", y = "ffr_value", color = "ffr_method", 
              line.color = "#bebebec7", line.size = 0.4, palette = c("#193bac", "#ebc90b")) + 
  stat_compare_means(paired = TRUE) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  scale_y_continuous(breaks = seq(0.4, 1, 0.1))

scatter_ffr <- ggplot(ffr_df, aes(x = ffr_method, y = ffr_value)) +
  geom_hline(yintercept = 0.8, linetype = "dashed", color = "red") +
  geom_boxplot(color = c("#193bac", "#ebc90b")) +
  geom_jitter(width = 0.05, color = "#bebebec7") +
  theme_classic() +
  scale_y_continuous(breaks = seq(0.4, 1, 0.05)) + 
  ylab("FFR value")

############################################################################################################
# check that area change was correlated with ffr change
# FFR
plot_ffr_meanao <- ggplot(baseline, aes(x = inv_aomean_change, y = inv_ffr_change)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Mean pressure change", y = "FFR change") +
  theme_classic()

plot_ffr_sysao <- ggplot(baseline, aes(x = inv_aosys_change, y = inv_ffr_change)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "Systolic pressure change", y = "FFR change") +
    theme_classic()

plot_ffr_diaao <- ggplot(baseline, aes(x = inv_aodia_change, y = inv_ffr_change)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "Diastolic pressure change", y = "FFR change") +
    theme_classic()

plot_ffr_hr <- ggplot(baseline, aes(x = inv_hr_change, y = inv_ffr_change)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "Heart rate change", y = "FFR change") +
    theme_classic()

# save plots
ggsave(paste0(output_dir,"/plot_ffr_change.png"), plot_ffr_change, width = 6, height = 4)
ggsave(paste0(output_dir,"/scatter_ffr.png"), scatter_ffr, width = 4, height = 4)
ggsave(paste0(output_dir,"/plot_ffr_meanao.png"), plot_ffr_meanao, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ffr_sysao.png"), plot_ffr_sysao, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ffr_diaao.png"), plot_ffr_diaao, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ffr_hr.png"), plot_ffr_hr, width = 6, height = 4)

baseline_new <- baseline %>% filter(!is.na(inv_ivusrest_mla))
shapiro.test(baseline_new$inv_ffrado)
med_ffrado <- median(baseline_new$inv_ffrado, na.rm = TRUE)
quant_ffrado <- quantile(baseline_new$inv_ffrado, probs = c(0.25, 0.75), na.rm = TRUE)
shapiro.test(baseline_new$inv_ffrdobu)
med_ffrdobu <- median(baseline_new$inv_ffrdobu, na.rm = TRUE)
quant_ffrdobu <- quantile(baseline_new$inv_ffrdobu, probs = c(0.25, 0.75), na.rm = TRUE)
print(paste0("Median FFR adenosine: ", med_ffrado, " (", quant_ffrado[1], "-", quant_ffrado[2], ")"))
print(paste0("Median FFR dobutamine: ", med_ffrdobu, " (", quant_ffrdobu[1], "-", quant_ffrdobu[2], ")"))
print(wilcox.test(baseline_new$inv_ffrado, baseline_new$inv_ffrdobu, paired = TRUE))