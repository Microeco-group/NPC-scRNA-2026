# Script: 01_alpha_diversity.R
# Purpose: Alpha-diversity estimation and plotting
# Workflow: microbiome
# Required inputs: otu_table_end.csv
# Main outputs: Plot_alpha.pdf
# Prerequisites: Run from this workflow directory so relative paths resolve correctly.

library(otuSummary)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ggsignif)
library(gghalves)

otu <- read.csv("otu_table_end.csv", header = T, row.names = 1, check.names = F )

alpha <- alphaDiversity(
  otutab = otu,         
  siteInCol = TRUE,          
  threshold = 1,                
  percent = FALSE,       
  write = F          
)

head(alpha)

data <- alpha$allBio
data$group <- c(rep("High",7),rep("Low",7))
head(data)
#write.csv(data,"16Sdata_alpha.csv")

if(T){
  mytheme <- theme(
    axis.title = element_text(size = 12,color ="black"), 
    axis.text = element_text(size= 12,color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1 ),
    legend.text = element_text(size= 12),
    legend.title= element_text(size= 12)
  ) }

ggplot(data, aes(x = group, y = observed)) + 
  geom_boxplot(aes(fill = group, color = group),  
               position = position_dodge(0.5), 
               width = 0.5, 
               outlier.alpha = 0) + 
  scale_fill_manual(values = scales::alpha(c("#e780a6", "#71c3ee"), alpha = 0.6)) + 
  scale_color_manual(values = c("#e780a6", "#71c3ee")) + 
  theme_bw() + 
  mytheme + 
  geom_signif(comparisons = list(c("High","Low")), 
              map_signif_level = F, 
              test = wilcox.test, 
              size=0.8,color="black")+
  theme(plot.margin = unit(c(0.5, 0.5, 0.5, 1), "cm"),) + 
  geom_jitter(position = position_jitterdodge(
    jitter.width = 0.4, 
    dodge.width = 0.85  
  ), 
  aes(fill = group), 
  color = "white",   
  size = 2,        
  shape = 21,       
  stroke = 0)+
  scale_y_log10()+ 
  labs(x="")+
  theme(legend.position = 'none')

ggsave("Plot_alpha.pdf",width = 6, height = 4)
