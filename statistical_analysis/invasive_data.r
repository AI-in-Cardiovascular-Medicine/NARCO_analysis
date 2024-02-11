library(tidyverse)
library(ggplot2)
library(ggpubr)
library(yaml)

yaml <- yaml.load_file(file.path(getwd(), "config.yaml"))
output_dir <- yaml$demographics$output_dir_figures
# Load the data
baseline <- readRDS(paste0(yaml$demographics$output_dir_data,"/baseline.rds"))

shapiro.test(baseline$inv_ffrado)
shapiro.test(baseline$inv_ffrdobu)
wilcox.test(baseline$inv_ffrado, baseline$inv_ffrdobu, paired = TRUE)

shapiro.test(baseline$inv_ivusrest_mla)
shapiro.test(baseline$inv_ivusdobu_mla)
t.test(baseline$inv_ivusrest_mla, baseline$inv_ivusdobu_mla, paired = TRUE)

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

# prep data for plotting ivus
n <- nrow(baseline)
inv_ivusrest <- pull(baseline, inv_ivusrest_mla)
inv_ivusdobu <- pull(baseline, inv_ivusdobu_mla)
method_1 <- rep("IVUS rest", n)
method_2 <- rep("IVUS dobutamine", n)
ivus_df <- data.frame(ivus_method = c(method_1, method_2), 
                    ivus_mla = c(inv_ivusrest, inv_ivusdobu))
colnames(ivus_df) <- c("ivus_method", "ivus_mla")

plot_ivus_change <- ggpaired(ivus_df, x = "ivus_method", y = "ivus_mla", color = "ivus_method", 
              line.color = "#bebebec7", line.size = 0.4, palette = c("#193bac", "#ebc90b")) + 
  stat_compare_means(paired = TRUE) +
  scale_y_continuous(breaks = seq(1, 20, 1))

# check that area change was correlated with ffr change
baseline <- baseline %>% mutate(
    inv_aomean_change = inv_dobu_aomean - inv_rest_aomean,
    inv_aosys_change = inv_dobu_aosys - inv_rest_aosys,
    inv_aodia_change = inv_dobu_aodia - inv_rest_aodia,
    inv_hr_change = inv_dobu_hr - inv_rest_hr,
    inv_ffr_change = inv_ffrado - inv_ffrdobu,
    inv_ivus_change = inv_ivusrest_mla - inv_ivusdobu_mla
)

# IVUS MLA
plot_ivus_meanao <- ggplot(baseline, aes(x = inv_aomean_change, y = inv_ivus_change)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  labs(x = "Mean pressure change", y = "IVUS MLA change") +
  theme_classic()

plot_ivus_sysao <- ggplot(baseline, aes(x = inv_aosys_change, y = inv_ivus_change)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "Systolic pressure change", y = "IVUS change") +
    theme_classic()

plot_ivus_diaao <- ggplot(baseline, aes(x = inv_aodia_change, y = inv_ivus_change)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "Diastolic pressure change", y = "IVUS change") +
    theme_classic()

plot_ivus_hr <- ggplot(baseline, aes(x = inv_hr_change, y = inv_ivus_change)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    labs(x = "Heart rate change", y = "IVUS change") +
    theme_classic()

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
ggsave(paste0(output_dir, "/plot_ffr_change.png"), plot_ffr_change, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ivus_change.png"), plot_ivus_change, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ivus_meanao.png"), plot_ivus_meanao, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ivus_sysao.png"), plot_ivus_sysao, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ivus_diaao.png"), plot_ivus_diaao, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ivus_hr.png"), plot_ivus_hr, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ffr_meanao.png"), plot_ffr_meanao, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ffr_sysao.png"), plot_ffr_sysao, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ffr_diaao.png"), plot_ffr_diaao, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ffr_hr.png"), plot_ffr_hr, width = 6, height = 4)