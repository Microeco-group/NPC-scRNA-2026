# Script: 02_beta_diversity.R
# Purpose: Beta-diversity and ordination analysis
# Workflow: microbiome
# Required inputs: physeq.end.rds; metadata.txt
# Main outputs: Beta_PCoA_uunifrac.pdf
# Prerequisites: Requires a prepared phyloseq object and sample metadata.

rm(list=ls())#clear Global Environment

library(vegan)
library(ggplot2)
library(ggrepel)
library(patchwork)
library(dplyr)
physeq.end <- readRDS("physeq.end.rds")
jaccard <- distance(physeq.end,method = "jaccard", binary = TRUE)

pcoa <- cmdscale (jaccard,eig=TRUE)  
pc12 <- pcoa$points[,1:2]
pc <- round(pcoa$eig/sum(pcoa$eig)*100,digits=2)

pc12 <- as.data.frame(pc12)
pc12$samples <- row.names(pc12)
head(pc12)

group <- read.delim('metadata.txt',header = T, stringsAsFactors = F)
colnames(group) <- c("samples","group")
df3 <- merge(pc12,group,by="samples")
head(df3)

shapes <- c(16, 17, 18, 22)
color=c("#e780a6","#71c3ee")

df3_means <- df3 %>%
  group_by(group) %>%
  summarise(mean_V1 = mean(V1), mean_V2 = mean(V2), .groups = 'drop')
df3 <- df3 %>%
  left_join(df3_means, by = "group")


p6 <- ggplot(data = df3, aes(x = V1, y = V2, color = group, shape = group)) +
  stat_ellipse(data = df3, geom = "polygon", level = 0.95, linetype = 2, size = 0, type = "t",
               aes(fill = group), alpha = 0.2) +  
  scale_fill_manual(values = c("#e780a6", "#71c3ee")) +  
  scale_color_manual(values = c("#e780a6", "#71c3ee")) +  
  theme(axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12, angle = 90),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 10),
        panel.grid = element_blank()) + 
  theme_bw() +  
  geom_point(size = 2.8) +  
  theme(panel.grid = element_blank(),
        aspect.ratio = 1) +
  geom_vline(xintercept = 0, lty = "dashed") + 
  geom_hline(yintercept = 0, lty = "dashed") + 
  labs(x = paste0("PC1 ", pc[1], "%"),
       y = paste0("PC2 ", pc[2], "%")) +  
  geom_segment(aes(x = mean_V1, y = mean_V2, xend = V1, yend = V2, color = group), 
               size = 0.5, alpha = 0.6, linetype = "solid")+  scale_color_manual(values = color) +
  scale_fill_manual(values = color)+
  theme(axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12,angle=90),
        axis.text.y=element_text(size=10),
        axis.text.x=element_text(size=10),
        panel.grid=element_blank(),
        aspect.ratio = 1)
p6

plot_build <- ggplot_build(p6)
x_range <- plot_build$layout$panel_scales_x[[1]]$range$range
x_range_min <- x_range[1]
x_range_max <- x_range[2]


plot_build <- ggplot_build(p6)
y_range <- plot_build$layout$panel_scales_y[[1]]$range$range
y_range_min <- y_range[1]
y_range_max <- y_range[2]


p7 <- ggplot(data = df3) +
  geom_density(aes(x = V1, fill=group), 
               color = 'black', alpha = 0.66, position = 'identity') +
  scale_linetype_manual(values = 'solid')+
  scale_x_continuous(limits = c(x_range_min, x_range_max)) 
p7

p8 <- p7 +
  scale_fill_manual(values=c("#e780a6","#71c3ee")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'right',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())
p8

p9 <- ggplot(data = df3) +
  geom_density(aes(x = V2, fill=group),
               color = 'black', alpha = 0.66, position = 'identity') +
  scale_linetype_manual(values = 'solid') +
  scale_fill_manual(values=c("#e780a6","#71c3ee")) +
  theme_classic() +
  theme(legend.title = element_blank(),
        legend.position = 'right',
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  scale_x_continuous(limits = c(y_range_min, y_range_max))+ 
  coord_flip() 
p9


pp <- p8 + plot_spacer() + p6 + p9 +
  plot_layout(guides = 'collect') &
  theme(legend.position='right')
pp
pp + plot_layout(widths = c(5, 1), heights = c(1,5))
pp_end <- pp + plot_layout(widths = c(5, 1), heights = c(1,5))

Adonis <- adonis2(jaccard~group, data=df3,permutations = 999)
Adonis

ggsave("Beta_PCoA_uunifrac.pdf", pp_end , width = 6, height = 6, units = "in")
