# Script: 06_correlation_analysis.R
# Purpose: Correlation analysis between bacterial features and biomarkers
# Workflow: microbiome
# Required inputs: Cor_IBL&ThreeGenus_merge.csv
# Main outputs: Cor_IBL&ThreeGenus_rel.pdf
# Prerequisites: Requires a prepared merged analysis table.

rm(list = ls())

library(ggplot2)

data_merge <- read.csv("Cor_IBL&ThreeGenus_merge.csv", row.names = 1)
colnames(data_merge)
data_select <- subset(data_merge, type == "abs")
data_select <- subset(data_merge, type == "rel")

data_select$project <- factor(data_select$project, levels = c("NPCcohort1",
                                                              "NPCcohort3_qPCR",
                                                              "NPCcohort3_16S",
                                                              "NPCcohort4_qPCR",
                                                              "NPCcohort4_16S",
                                                              "COAD",
                                                              "ESCA",
                                                              "HNSC",
                                                              "READ",
                                                              "STAD"))

data_select$Variable <- factor(data_select$Variable, levels = c("Fusobacterium_abs",
                                                                "Campylobacter_abs",
                                                                "Parvimonas_abs"))
data_select$Variable <- factor(data_select$Variable, levels = c("Fusobacterium_rel",
                                                                "Campylobacter_rel",
                                                                "Parvimonas_rel"))


ggplot(data_select, aes(x = project, y = cor, fill = IBL_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("#f484ae", "#815aa8","#2ba8a5")) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylim(0,1)+
  facet_wrap(~ Variable)+
  geom_text(aes(label = paste0("P=", format(pvalue, scientific = TRUE, digits = 2)),
                y = 0), 
            position = position_dodge(width = 0.6),
            vjust = 0.5, 
            color = "black",
            angle = 90, 
            hjust = 0)

ggplot(data_select, aes(x = Variable, y = cor, fill = Variable)) +
  geom_bar(stat = "identity", width = 0.6) +
  scale_fill_manual(values = c("#f484ae", "#815aa8","#2ba8a5")) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  ylim(0,0.85)+
  facet_wrap(~ project, nrow = 1)+
  geom_text(aes(label = paste0("P=", format(pvalue, scientific = TRUE, digits = 2)),
                y = 0),  
            position = position_dodge(width = 0.6),
            vjust = 0.5, 
            color = "black",
            angle = 90,
            hjust = 0)
ggsave("Cor_IBL&ThreeGenus_rel.pdf", width = 10, height = 3)
