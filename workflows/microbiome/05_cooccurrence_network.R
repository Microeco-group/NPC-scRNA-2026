# Script: 05_cooccurrence_network.R
# Purpose: Co-occurrence network construction by group
# Workflow: microbiome
# Required inputs: otu_table_end.csv; tax_table_end.csv
# Main outputs: Co_Network_High&Low.pdf
# Prerequisites: Requires OTU and taxonomy tables.

rm(list = ls())
library(tidyverse)
library(ggNetView)
library(patchwork)
library(ggplot2)
library(microeco)
library(tidyverse)

otu <- read.csv("otu_table_end.csv",row.names = 1, check.names = F)

otu_All <- otu[rowSums(otu) != 0, ]
rel_ab <- prop.table(as.matrix(otu_All), 2)
filtered_otu_table <- rel_ab[rowMeans(rel_ab) > 0.0001, ]

otu_High <- filtered_otu_table[,1:7]
otu_High <- otu_High[rowSums(otu_High) != 0, ]
otu_Low <- filtered_otu_table[,8:14]
otu_Low <- otu_Low[rowSums(otu_Low) != 0, ]

tax <- subset(tax, subset = OTUID %in% rownames(filtered_otu_table))

tax_table <- read.csv("tax_table_end.csv",row.names = 1)
tax_table <- tax_table[,1:7]
colnames(tax_table) <- c("Domain","Phylum","Class","Order","Family","Genus","Species")
tax_table <- tidy_taxonomy(tax_table)
tax <- data.frame("OTUID" = rownames(tax_table))
tax[2:8] <- tax_table[,1:7]
tax <- subset(tax, subset = OTUID %in% rownames(filtered_otu_table))

r_threshold <- 0.7

mycolors <- c(
  '#E0D4CA', '#CCC9E6', '#D6E7A3', '#F3B1A0', '#53A85F', '#F1BB72', '#E95C59',
  '#57C3F3', '#E5D2DD', '#AB3282', '#91D0BE', '#625D9E', '#BD956A', '#C5DEBA',
  '#CCE0F5', '#F7F398', '#C1E6F3')
names(mycolors) <- sort(unique(tax$Phylum))

obj_High <- build_graph_from_mat(
  mat = otu_High,
  transfrom.method = "none",
  method = "WGCNA",
  cor.method = "spearman",
  r.threshold = r_threshold,
  p.threshold = 0.01,
  node_annotation = tax,
  proc = "BH"
)
obj_Low <- build_graph_from_mat(
  mat = otu_Low,
  transfrom.method = "none",
  method = "WGCNA",
  cor.method = "spearman",
  r.threshold = r_threshold,
  p.threshold = 0.01,
  node_annotation = tax,
  proc = "BH"
)
p_High <- ggNetView(
  graph_obj = obj_High,
  layout = "gephi",
  layout.module = "adjacent",
  group.by = "Phylum",
  center = FALSE,
  shrink = 0.9,
  linealpha = 0.35,
  linecolor = "#d9d9d9") + 
  scale_fill_manual(
    values = mycolors) + 
  scale_size(range = c(2, 10),limits = c(0,70))
p_Low <- ggNetView(
  graph_obj = obj_Low,
  layout = "gephi",
  layout.module = "adjacent",
  group.by = "Phylum",
  center = FALSE,
  shrink = 0.9,
  linealpha = 0.35,
  linecolor = "#d9d9d9") + 
  scale_fill_manual(
    values = mycolors) + 
  scale_size(range = c(2, 10),limits = c(0,70))
p_High+p_Low+plot_layout(guides = "collect", widths = c(1,0.5))
ggsave("Co_Network_High&Low.pdf", width = 14, height = 14)
