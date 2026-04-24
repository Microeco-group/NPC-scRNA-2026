# Script: 04_gsea_pathway_analysis.R
# Purpose: Pathway enrichment analysis from differential expression results
# Workflow: bulk-rna
# Required inputs: differ.csv or the exported DEG result table.
# Main outputs: GSEA_result.pdf and pathway result tables.
# Prerequisites: Requires a ranked differential expression result table.

library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(ggplot2)
library(GseaVis)

resdata <- read.csv("differ.csv", row.names = 1)

symbol  <- resdata$Gene
entrez  <- bitr(symbol, 
                fromType  = "SYMBOL",
                toType  = "ENTREZID",
                OrgDb  = "org.Hs.eg.db")
head(entrez)

genelist  <- resdata$log2FoldChange
names(genelist) <- resdata$Gene

genelist  <- genelist[names(genelist) %in% entrez[,1]]
names(genelist) <- entrez[match(names(genelist),entrez[,1]),2]
length(genelist)
head(genelist)
genelist <- sort(genelist, decreasing = T)
head(genelist)
KEGG_ges  <- gseKEGG(
  geneList  = genelist,
  organism  = "hsa",
  minGSSize  = 10,
  maxGSSize  = 500,
  pvalueCutoff  = 0.05,
  pAdjustMethod  = "BH",
  verbose  = FALSE,
  eps  = 0
)

KEGG_ges  <- setReadable(KEGG_ges,
                         OrgDb  = org.Hs.eg.db,
                         keyType  = "ENTREZID")
KEGG_ges_result  <- KEGG_ges@result

library(GseaVis)
terms <- c("hsa04064")
gseaNb(object = KEGG_ges,
       geneSetID = terms,
       subPlot = 3,
       htCol = c("#6DCCFD","#F484AE"),
       curveCol = c("#F484AE"),
       addPval = T,
       pDigit = 100)
ggsave("GSEA_result.pdf", width = 4, height = 4)


