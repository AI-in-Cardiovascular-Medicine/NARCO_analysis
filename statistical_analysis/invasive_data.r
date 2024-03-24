library(tidyverse)
library(ggplot2)
library(ggpubr)
library(yaml)

yaml <- yaml.load_file(file.path(getwd(), "config.yaml"))
output_dir <- yaml$demographics$output_dir_figures
# Load the data
baseline <- readRDS(paste0(yaml$demographics$output_dir_data,"/baseline.rds"))


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

########################################################################################
# prep data for plotting ivus mla
n <- nrow(baseline)
inv_ivusrest <- pull(baseline, inv_ivusrest_mla)
inv_ffrado <- pull(baseline, inv_ffrado)
inv_ivusdobu <- pull(baseline, inv_ivusdobu_mla)
inv_ffrdobu <- pull(baseline, inv_ffrdobu)
method_1 <- rep("IVUS rest", n)
method_2 <- rep("IVUS stress", n)
ivus_df <- data.frame(ivus_method = c(method_1, method_2), 
                    ivus_mla = c(inv_ivusrest, inv_ivusdobu),
                    ffr_value = c(inv_ffrado, inv_ffrdobu))
ivus_df <- ivus_df %>% mutate(ffr_value = ifelse(ffr_value <= 0.8, 1, 0))
ivus_df$ffr_value <- factor(ivus_df$ffr_value, levels = c(0, 1), labels = c("FFR > 0.8", "FFR <= 0.8"))
colnames(ivus_df) <- c("ivus_method", "ivus_mla", "ffr_value")

plot_ivus_mla <- ggpaired(ivus_df, ivus_df, x = "ivus_method", y = "ivus_mla", color = "ivus_method", 
                        line.color = "#bebebec7", line.size = 0.4, palette = c("red", "#bebebec7","#193bac", "#ebc90b")) +
  geom_point(aes(color = ffr_value)) +
  stat_compare_means(paired = TRUE) +
  scale_y_continuous(breaks = seq(1, 20, 1)) +
  ylab(expression(paste("IVUS MLA [", mm^2, "]")))

ivus_df$ivus_method <- factor(ivus_df$ivus_method, levels = c("IVUS rest", "IVUS stress"))
scatter_mla <- ggplot(ivus_df, aes(x = ivus_method, y = ivus_mla)) +
  geom_boxplot(color = c("#193bac", "#ebc90b")) +
  geom_jitter(width = 0.05, color = "#bebebec7") +
  theme_classic() +
  scale_y_continuous(breaks = seq(1, 20, 1)) + 
  ylab(expression(paste("IVUS MLA [", mm^2, "]")))

# ggplot(ivus_df, aes(x = ivus_method, y = ivus_mla)) +
#   geom_boxplot(color = c("#193bac", "#ebc90b"), width = 0.5) +
#   geom_violin(fill = "grey", alpha = 0.3) +
#   geom_jitter(width = 0.05, color = "#bebebec7") +
#   theme_classic() +
#   scale_y_continuous(breaks = seq(1, 20, 1)) + 
#   ylab(expression(paste("IVUS MLA [", mm^2, "]")))

# ivus ostial a
n <- nrow(baseline)
inv_ivusrest <- pull(baseline, inv_ivusrest_ostial_a)
inv_ffrado <- pull(baseline, inv_ffrado)
inv_ivusdobu <- pull(baseline, inv_ivusdobu_ostial_a)
inv_ffrdobu <- pull(baseline, inv_ffrdobu)
method_1 <- rep("IVUS rest", n)
method_2 <- rep("IVUS stress", n)
ivus_df <- data.frame(ivus_method = c(method_1, method_2), 
                    ivus_mla = c(inv_ivusrest, inv_ivusdobu),
                    ffr_value = c(inv_ffrado, inv_ffrdobu))
ivus_df <- ivus_df %>% mutate(ffr_value = ifelse(ffr_value <= 0.8, 1, 0))
ivus_df$ffr_value <- factor(ivus_df$ffr_value, levels = c(0, 1), labels = c("FFR > 0.8", "FFR <= 0.8"))
colnames(ivus_df) <- c("ivus_method", "ivus_ostial_a", "ffr_value")

plot_ivus_ostial_a <- ggpaired(ivus_df, ivus_df, x = "ivus_method", y = "ivus_ostial_a", color = "ivus_method", 
                        line.color = "#bebebec7", line.size = 0.4, palette = c("red", "#bebebec7","#193bac", "#ebc90b")) +
  geom_point(aes(color = ffr_value)) +
  stat_compare_means(paired = TRUE) +
  scale_y_continuous(breaks = seq(1, 20, 1)) +
  ylab(expression(paste("IVUS ostial area [", mm^2, "]")))

ivus_df$ivus_method <- factor(ivus_df$ivus_method, levels = c("IVUS rest", "IVUS stress"))
scatter_ostial_a <- ggplot(ivus_df, aes(x = ivus_method, y = ivus_ostial_a)) +
  geom_boxplot(color = c("#193bac", "#ebc90b")) +
  geom_jitter(width = 0.05, color = "#bebebec7") +
  theme_classic() +
  scale_y_continuous(breaks = seq(1, 20, 1)) + 
  ylab(expression(paste("IVUS ostial area [", mm^2, "]")))

########################################################################################
# data for mla ellipticity
n <- nrow(baseline)
inv_ivusrest_ellip <- pull(baseline, inv_ivusrest_mla_ellip)
inv_ffrado <- pull(baseline, inv_ffrado)
inv_ivusdobu_ellip <- pull(baseline, inv_ivusdobu_mla_ellip)
inv_ffrdobu <- pull(baseline, inv_ffrdobu)
method_1 <- rep("IVUS rest", n)
method_2 <- rep("IVUS dobutamine", n)
ivus_df <- data.frame(ivus_method = c(method_1, method_2), 
                    ivus_mla_ellip = c(inv_ivusrest_ellip, inv_ivusdobu_ellip),
                    ffr_value = c(inv_ffrado, inv_ffrdobu))
ivus_df <- ivus_df %>% mutate(ffr_value = ifelse(ffr_value <= 0.8, 1, 0))
ivus_df$ffr_value <- factor(ivus_df$ffr_value, levels = c(0, 1), labels = c("FFR > 0.8", "FFR <= 0.8"))
colnames(ivus_df) <- c("ivus_method", "ivus_mla_ellipticity", "ffr_value")

plot_ivus_mla_ellip <- ggpaired(ivus_df, ivus_df, x = "ivus_method", y = "ivus_mla_ellipticity", color = "ivus_method", line.color = "#bebebec7", line.size = 0.4, palette = c("red", "#bebebec7","#193bac", "#ebc90b")) +
  geom_point(aes(color = ffr_value)) +
  stat_compare_means(paired = TRUE) +
  scale_y_continuous(breaks = seq(1, 20, 1)) +
  ylab("IVUS MLA elliptic ratio")

ivus_df$ivus_method <- factor(ivus_df$ivus_method, levels = c("IVUS rest", "IVUS dobutamine"))
scatter_mla_ellip <- ggplot(ivus_df, aes(x = ivus_method, y = ivus_mla_ellipticity)) +
  geom_boxplot(color = c("#193bac", "#ebc90b")) +
  geom_jitter(width = 0.05, color = "#bebebec7") +
  theme_classic() +
  scale_y_continuous(breaks = seq(1, 20, 1)) + 
  ylab("IVUS MLA elliptic ratio")

# data for ostial ellipticity
n <- nrow(baseline)
inv_ivusrest_ellip <- pull(baseline, inv_ivusrest_ostial_ellip)
inv_ffrado <- pull(baseline, inv_ffrado)
inv_ivusdobu_ellip <- pull(baseline, inv_ivusdobu_ostial_ellip)
inv_ffrdobu <- pull(baseline, inv_ffrdobu)
method_1 <- rep("IVUS rest", n)
method_2 <- rep("IVUS dobutamine", n)
ivus_df <- data.frame(ivus_method = c(method_1, method_2), 
                    ivus_mla = c(inv_ivusrest_ellip, inv_ivusdobu_ellip),
                    ffr_value = c(inv_ffrado, inv_ffrdobu))
ivus_df <- ivus_df %>% mutate(ffr_value = ifelse(ffr_value <= 0.8, 1, 0))
ivus_df$ffr_value <- factor(ivus_df$ffr_value, levels = c(0, 1), labels = c("FFR > 0.8", "FFR <= 0.8"))
colnames(ivus_df) <- c("ivus_method", "ivus_ostial_ellipticity", "ffr_value")

plot_ivus_ostial_ellip <- ggpaired(ivus_df, ivus_df, x = "ivus_method", y = "ivus_ostial_ellipticity", color = "ivus_method", line.color = "#bebebec7", line.size = 0.4, palette = c("red", "#bebebec7","#193bac", "#ebc90b")) +
  geom_point(aes(color = ffr_value)) +
  stat_compare_means(paired = TRUE) +
  scale_y_continuous(breaks = seq(1, 20, 1)) +
  ylab("IVUS ostial elliptic ratio")

ivus_df$ivus_method <- factor(ivus_df$ivus_method, levels = c("IVUS rest", "IVUS dobutamine"))
scatter_ostial_ellip <- ggplot(ivus_df, aes(x = ivus_method, y = ivus_ostial_ellipticity)) +
  geom_boxplot(color = c("#193bac", "#ebc90b")) +
  geom_jitter(width = 0.05, color = "#bebebec7") +
  theme_classic() +
  scale_y_continuous(breaks = seq(1, 20, 1)) + 
  ylab("IVUS ostial elliptic ratio")

############################################################################################################
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
ggsave(paste0(output_dir,"/scatter_ffr.png"), scatter_ffr, width = 4, height = 4)
ggsave(paste0(output_dir,"/plot_ivus_mla.png"), plot_ivus_mla, width = 6, height = 4)
ggsave(paste0(output_dir,"/scatter_mla.png"), scatter_mla, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ivus_ostial_a.png"), plot_ivus_ostial_a, width = 6, height = 4)
ggsave(paste0(output_dir,"/scatter_ostial_a.png"), scatter_ostial_a, width = 4, height = 4)
ggsave(paste0(output_dir,"/plot_ivus_mla_ellip.png"), plot_ivus_mla_ellip, width = 6, height = 4)
ggsave(paste0(output_dir,"/scatter_mla_ellip.png"), scatter_mla_ellip, width = 5, height = 5)
ggsave(paste0(output_dir,"/plot_ivus_ostial_ellip.png"), plot_ivus_ostial_ellip, width = 6, height = 4)
ggsave(paste0(output_dir,"/scatter_ostial_ellip.png"), scatter_ostial_ellip, width = 4, height = 4)
ggsave(paste0(output_dir,"/plot_ivus_meanao.png"), plot_ivus_meanao, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ivus_sysao.png"), plot_ivus_sysao, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ivus_diaao.png"), plot_ivus_diaao, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ivus_hr.png"), plot_ivus_hr, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ffr_meanao.png"), plot_ffr_meanao, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ffr_sysao.png"), plot_ffr_sysao, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ffr_diaao.png"), plot_ffr_diaao, width = 6, height = 4)
ggsave(paste0(output_dir,"/plot_ffr_hr.png"), plot_ffr_hr, width = 6, height = 4)

baseline_new <- baseline %>% filter(!is.na(inv_ivusrest_mla))
shapiro.test(baseline_new$inv_ffrado)
median(baseline_new$inv_ffrado, na.rm = TRUE)
quantile(baseline_new$inv_ffrado, probs = c(0.25, 0.75), na.rm = TRUE)
shapiro.test(baseline_new$inv_ffrdobu)
median(baseline_new$inv_ffrdobu, na.rm = TRUE)
quantile(baseline_new$inv_ffrdobu, probs = c(0.25, 0.75), na.rm = TRUE)
wilcox.test(baseline_new$inv_ffrado, baseline_new$inv_ffrdobu, paired = TRUE)

shapiro.test(baseline_new$inv_rest_hr)
shapiro.test(baseline_new$inv_dobu_hr)
median(baseline_new$inv_rest_hr, na.rm = TRUE)
quantile(baseline_new$inv_rest_hr, probs = c(0.25, 0.75), na.rm = TRUE)
median(baseline_new$inv_dobu_hr, na.rm = TRUE)
quantile(baseline_new$inv_dobu_hr, probs = c(0.25, 0.75), na.rm = TRUE)
wilcox.test(baseline_new$inv_rest_hr, baseline_new$inv_dobu_hr, paired = TRUE)

shapiro.test(baseline_new$inv_rest_aosys)
shapiro.test(baseline_new$inv_dobu_aosys)
median(baseline_new$inv_rest_aosys, na.rm = TRUE)
quantile(baseline_new$inv_rest_aosys, probs = c(0.25, 0.75), na.rm = TRUE)
median(baseline_new$inv_dobu_aosys, na.rm = TRUE)
quantile(baseline_new$inv_dobu_aosys, probs = c(0.25, 0.75), na.rm = TRUE)
wilcox.test(baseline_new$inv_rest_aosys, baseline_new$inv_dobu_aosys, paired = TRUE)

shapiro.test(baseline_new$inv_rest_aodia)
shapiro.test(baseline_new$inv_dobu_aodia)
mean(baseline_new$inv_rest_aodia, na.rm = TRUE)
sd(baseline_new$inv_rest_aodia, na.rm = TRUE)
mean(baseline_new$inv_dobu_aodia, na.rm = TRUE)
sd(baseline_new$inv_dobu_aodia, na.rm = TRUE)
t.test(baseline_new$inv_rest_aodia, baseline_new$inv_dobu_aodia, paired = TRUE)

shapiro.test(baseline_new$inv_rest_aomean)
shapiro.test(baseline_new$inv_dobu_aomean)
mean(baseline_new$inv_rest_aomean, na.rm = TRUE)
sd(baseline_new$inv_rest_aomean, na.rm = TRUE)
mean(baseline_new$inv_dobu_aomean, na.rm = TRUE)
sd(baseline_new$inv_dobu_aomean, na.rm = TRUE)
t.test(baseline_new$inv_rest_aomean, baseline_new$inv_dobu_aomean, paired = TRUE)

shapiro.test(baseline_new$inv_ivusrest_mla)
median(baseline_new$inv_ivusrest_mla, na.rm = TRUE)
quantile(baseline_new$inv_ivusrest_mla, probs = c(0.25, 0.75), na.rm = TRUE)
shapiro.test(baseline_new$inv_ivusdobu_mla)
median(baseline_new$inv_ivusdobu_mla, na.rm = TRUE)
quantile(baseline_new$inv_ivusdobu_mla, probs = c(0.25, 0.75), na.rm = TRUE)
wilcox.test(baseline_new$inv_ivusrest_mla, baseline_new$inv_ivusdobu_mla, paired = TRUE)

shapiro.test(baseline_new$inv_ivusrest_mla_ellip)
median(baseline_new$inv_ivusrest_mla_ellip, na.rm = TRUE)
quantile(baseline_new$inv_ivusrest_mla_ellip, probs = c(0.25, 0.75), na.rm = TRUE)
shapiro.test(baseline_new$inv_ivusdobu_mla_ellip)
median(baseline_new$inv_ivusdobu_mla_ellip, na.rm = TRUE)
quantile(baseline_new$inv_ivusdobu_mla_ellip, probs = c(0.25, 0.75), na.rm = TRUE)
wilcox.test(baseline_new$inv_ivusrest_mla_ellip, baseline_new$inv_ivusdobu_mla_ellip, paired = TRUE)

shapiro.test(baseline_new$inv_ivusrest_mla_w)
median(baseline_new$inv_ivusrest_mla_w, na.rm = TRUE)
quantile(baseline_new$inv_ivusrest_mla_w, probs = c(0.25, 0.75), na.rm = TRUE)
shapiro.test(baseline_new$inv_ivusdobu_mla_w)
median(baseline_new$inv_ivusdobu_mla_w, na.rm = TRUE)
quantile(baseline_new$inv_ivusdobu_mla_w, probs = c(0.25, 0.75), na.rm = TRUE)
wilcox.test(baseline_new$inv_ivusrest_mla_w, baseline_new$inv_ivusdobu_mla_w, paired = TRUE)

shapiro.test(baseline_new$inv_ivusrest_mla_h)
shapiro.test(baseline_new$inv_ivusdobu_mla_h)
mean(baseline_new$inv_ivusrest_mla_h, na.rm = TRUE)
mean(baseline_new$inv_ivusdobu_mla_h, na.rm = TRUE)
sd(baseline_new$inv_ivusrest_mla_h, na.rm = TRUE)
sd(baseline_new$inv_ivusdobu_mla_h, na.rm = TRUE)
t.test(baseline_new$inv_ivusrest_mla_h, baseline_new$inv_ivusdobu_mla_h, paired = TRUE)

shapiro.test(baseline_new$inv_ivusrest_ostial_a)
shapiro.test(baseline_new$inv_ivusdobu_ostial_a)
median(baseline_new$inv_ivusrest_ostial_a, na.rm = TRUE)
quantile(baseline_new$inv_ivusrest_ostial_a, probs = c(0.25, 0.75), na.rm = TRUE)
median(baseline_new$inv_ivusdobu_ostial_a, na.rm = TRUE)
quantile(baseline_new$inv_ivusdobu_ostial_a, probs = c(0.25, 0.75), na.rm = TRUE)
wilcox.test(baseline_new$inv_ivusrest_ostial_a, baseline_new$inv_ivusdobu_ostial_a, paired = TRUE)

shapiro.test(baseline_new$inv_ivusrest_ostial_w)
shapiro.test(baseline_new$inv_ivusdobu_ostial_w)
mean(baseline_new$inv_ivusrest_ostial_w, na.rm = TRUE)
sd(baseline_new$inv_ivusrest_ostial_w, na.rm = TRUE)
mean(baseline_new$inv_ivusdobu_ostial_w, na.rm = TRUE)
sd(baseline_new$inv_ivusdobu_ostial_w, na.rm = TRUE)
t.test(baseline_new$inv_ivusrest_ostial_w, baseline_new$inv_ivusdobu_ostial_w, paired = TRUE)

shapiro.test(baseline_new$inv_ivusrest_ostial_h)
shapiro.test(baseline_new$inv_ivusdobu_ostial_h)
mean(baseline_new$inv_ivusrest_ostial_h, na.rm = TRUE)
sd(baseline_new$inv_ivusrest_ostial_h, na.rm = TRUE)
mean(baseline_new$inv_ivusdobu_ostial_h, na.rm = TRUE)
sd(baseline_new$inv_ivusdobu_ostial_h, na.rm = TRUE)
t.test(baseline_new$inv_ivusrest_ostial_h, baseline_new$inv_ivusdobu_ostial_h, paired = TRUE)

shapiro.test(baseline_new$inv_ivusrest_ostial_ellip)
shapiro.test(baseline_new$inv_ivusdobu_ostial_ellip)
median(baseline_new$inv_ivusrest_ostial_ellip, na.rm = TRUE)
quantile(baseline_new$inv_ivusrest_ostial_ellip, probs = c(0.25, 0.75), na.rm = TRUE)
median(baseline_new$inv_ivusdobu_ostial_ellip, na.rm = TRUE)
quantile(baseline_new$inv_ivusdobu_ostial_ellip, probs = c(0.25, 0.75), na.rm = TRUE)
wilcox.test(baseline_new$inv_ivusrest_ostial_ellip, baseline_new$inv_ivusdobu_ostial_ellip, paired = TRUE)

