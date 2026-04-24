# Script: 03_mcp_counter_abundance.R
# Purpose: MCPcounter immune abundance estimation
# Workflow: bulk-rna
# Required inputs: Bulk_batch_count.csv
# Main outputs: MCPcounter_Results.csv; box_MCPcounter.pdf
# Prerequisites: Typically run after batch correction.

rm(list = ls())
library(MCPcounter)

expr_matrix <- read.csv("Bulk_batch_count.csv", row.names = 1)
exp <- data.frame(apply(expr_matrix, 2, function(x){round(log2(x+1),2)})) #统一标准化
exp <- as.matrix(exp)
res<- MCPcounter.estimate( 
  exp, 
  featuresType='HUGO_symbols',  
  probesets=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/probesets.txt"),sep='\t',stringsAsFactors=FALSE,colClasses='character'),
  genes=read.table(curl("http://raw.githubusercontent.com/ebecht/MCPcounter/master/Signatures/genes.txt"),sep='\t',stringsAsFactors=FALSE,header=TRUE,colClasses='character',check.names=FALSE))
head(res[,1:4],n=3)
write.csv(res, "MCPcounter_Results.csv")

#res <- read.csv("MCPcounter_Results.csv", row.names = 1)
results <- data.frame(t(res))
results$sample <- rownames(results)
results <- data.frame(results,row.names = NULL)
results$group <- c(rep('High',40),rep('Low',62))
library(tidyr)
results <- pivot_longer(data = results,cols = colnames(results)[1:9],names_to = "cell.type",values_to = 'value')
results <- as.data.frame(results)
head(results,n=3)

library(ggplot2)
library(ggpubr)
library(ggsci)
library(lemon)

if(T){
  mytheme <- theme(
    axis.title = element_text(size = 12,color ="black"), 
    axis.text = element_text(size= 12,color = "black"),
    axis.text.x = element_text(angle = 45, hjust = 1 ),
    legend.text = element_text(size= 12),
    legend.title= element_text(size= 12)
  ) }

box_MCPcounter <- ggplot(results, aes(x = cell.type, y = value)) + 
  labs(y = "MCP Counter Score", x = NULL) +  
  geom_boxplot(aes(fill = group, color = group), 
               position = position_dodge(0.5), 
               width = 0.5, 
               outlier.alpha = 0) + 
  scale_fill_manual(values = scales::alpha(c("#F484AE", "#6DCCFD"), alpha = 0.6)) + 
  scale_color_manual(values = c("#F484AE", "#6DCCFD")) + 
  theme_bw() + 
  mytheme + 
  stat_compare_means(aes(group = group),
                     label = "p.format",
                     method = "wilcox.test",
                     hide.ns = FALSE) + 
  geom_jitter(position = position_jitterdodge(
    jitter.width = 0.25, 
    dodge.width = 0.5  
  ), 
  aes(fill = group), 
  color = "white",   
  size = 2,        
  shape = 21,       
  stroke = 0)+
  scale_y_continuous(expand=c(0.05,0.05))
box_MCPcounter
ggsave("box_MCPcounter.pdf",box_MCPcounter,height=4,width=6.5)
