# Script: 03_lefse_analysis.R
# Purpose: LEfSe differential taxa analysis
# Workflow: microbiome
# Required inputs: otu_table_end.csv; tax_table_end.csv; metadata.txt
# Main outputs: lefse_diff.csv; lefse.pdf
# Prerequisites: Requires OTU, taxonomy, and sample metadata tables.

rm(list=ls())

library(ggplot2)
library(microeco)
library(tidyverse)
library(magrittr)
library(ggtree)

feature_table <- read.csv("otu_table_end.csv",row.names = 1, check.names = F)
sample_table <- read.delim('metadata.txt',row.names = 1)
tax_table <- read.csv("tax_table_end.csv",row.names = 1)
tax_table <- tax_table[,1:7]
colnames(tax_table) <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
tax_table <- tidy_taxonomy(tax_table)

tax_table <- tidy_taxonomy(tax_table)
head(feature_table)[1:6,1:6]; head(sample_table)[1:6, ]; head(tax_table)[,1:2]

dataset <- microtable$new(sample_table = sample_table,
                          otu_table = feature_table, 
                          tax_table = tax_table)
dataset

lefse <- trans_diff$new(dataset = dataset, 
                        method = "lefse", 
                        group = "group", 
                        alpha = 0.05, 
                        p_adjust_method = "none",
                        taxa_level = "Genus",
                        lefse_subgroup = NULL)
df <- lefse$res_diff
df <- as.data.frame(lefse$res_diff)
write.csv(df,"lefse_diff.csv")

head(lefse$res_diff)
lefse$plot_diff_bar(use_number = 1:10, 
                    width = 0.8, 
                    group_order = c("High", "Low"),
                    keep_prefix = F) +
  scale_fill_manual(values = c("#e780a6","#71c3ee"))+
  scale_colour_manual(values = c("#e780a6","#71c3ee"))

g1 <- lefse$plot_diff_bar(use_number = 1:10, 
                          width = 0.8, 
                          group_order = c("High", "Low"),
                          keep_prefix = F) +
  scale_fill_manual(values = c("#e780a6","#71c3ee"))+
  scale_colour_manual(values = c("#e780a6","#71c3ee"))
g1 <- g1 + theme(legend.position = "none")
g1

g2 <- lefse$plot_diff_abund(group_order = c("High", "Low"), select_taxa = lefse$plot_diff_bar_taxa, add_sig = T, add_sig_label = "Significance",keep_prefix = F)+
  scale_fill_manual(values = c("#71c3ee","#e780a6"))+
  scale_colour_manual(values = c("#71c3ee","#e780a6"))
g2 <- g2 + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
g2
g <- g1+g2
g

ggsave("lefse.pdf",g,width = 8,height = 4,units = "in")