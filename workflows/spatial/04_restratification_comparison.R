# Script: 04_restratification_comparison.R
# Purpose: Restratification comparison plots for spatially defined groups
# Workflow: spatial
# Required inputs: high_low_mid_tex_mreg.csv
# Main outputs: Proportion_restratification.pdf
# Prerequisites: Requires a precomputed slide-level summary table.

library(ggpubr)
library(ggplot2)
library(dplyr)
rm(list = ls())
color_fov <- c("#F484AE", "#6DCCFD")
names(color_fov) <- c("High","Low")

slide_proportion <- read.csv("high_low_mid_tex_mreg.csv", row.names = 1)

new_df <- slide_proportion %>%
  group_by(Sample) %>%
  mutate(
    q25 = quantile(bacterial_por, 0.25, na.rm = TRUE),
    q75 = quantile(bacterial_por, 0.75, na.rm = TRUE),
    Group = case_when(
      bacterial_por < q25 ~ "Low",
      bacterial_por > q75 ~ "High",
      TRUE ~ "Mid"
    )
  ) %>%
  ungroup() %>%
  select(-q25, -q75)  

slide_proportion <- subset(new_df, subset = Group %in% c("High","Low"))
slide_proportion$Group <- factor(slide_proportion$Group, levels = c("High","Low"))
result <- slide_proportion %>%
  group_by(Group, Sample) %>%
  summarise(
    mean_CD8_Tex_por = mean(CD8_Tex_por, na.rm = TRUE),
    mean_mregDC_por = mean(mregDC_por, na.rm = TRUE),
    .groups = 'drop' 
  )


valid_samples <- result %>%
  group_by(Sample) %>%
  filter(n_distinct(Group) == 2) %>%
  pull(Sample) %>%
  unique()

paired_data <- result %>%
  filter(Sample %in% valid_samples)
paired_data$Group <- factor(paired_data$Group, levels = c("High","Low"))

p1 <- ggplot(paired_data, aes(x = Group, y = mean_CD8_Tex_por)) +
  geom_line(aes(group = Sample), color = "gray60", alpha = 0.6, size = 0.5) +
  geom_point(aes(color = Group), size = 3, alpha = 1) +
  stat_compare_means(paired = TRUE, label = "p.format", 
                     method = "wilcox.test") +
  labs(x = "Group", 
       y = "mean_CD8_Tex_por") +
  scale_color_manual(values = color_fov) +
  theme_bw()+
  theme(aspect.ratio = 2)

p2 <- ggplot(paired_data, aes(x = Group, y = mean_mregDC_por)) +
  geom_line(aes(group = Sample), color = "gray60", alpha = 0.6, size = 0.5) +
  geom_point(aes(color = Group), size = 3, alpha = 1) +
  stat_compare_means(paired = TRUE, label = "p.format", 
                     method = "wilcox.test") +
  labs(x = "Group", 
       y = "mean_mregDC_por") +
  scale_color_manual(values = color_fov) +
  theme_bw()+
  theme(aspect.ratio = 2)

library(patchwork)
p1+p2
ggsave("Proportion_restratification.pdf", width = 6, height = 4)