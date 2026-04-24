# Script: 07_ssgsea_scoring.R
# Purpose: ssGSEA scoring for custom signatures
# Workflow: bulk-rna
# Required inputs: Primary_exp.csv; signature.csv
# Main outputs: ssGSEA_result.csv
# Prerequisites: Independent from the DEG branch but requires expression data and signature definitions.

library(GSVA)
RNA_matrix <- read.csv("Primary_exp.csv")
mymatrix<-as.matrix(RNA_matrix)

Gene_name1="signature"  
geneset <- read.csv("signature.csv")
Markers <- geneset$gene
Markers <- subset(Markers, subset = Markers != "")
mysymbol <- data.frame(Gene_set= Gene_name1 ,Gene_symbol= Markers)

colnames(mysymbol)<-c("Gene_set","Gene_symbol")
head(mysymbol)
table(mysymbol$Gene_set)

type <- unique(mysymbol$Gene_set)
type
gs <- list()
for (i in type){
  tmp <- mysymbol$Gene_symbol[which(mysymbol$Gene_set == i)]
  tmp <- list(tmp)
  gs <- c(gs,tmp)
}
names(gs) <- type
gs

es.dif <- gsva(mymatrix, gs, method = "ssgsea", kcdf = "Gaussian", 
               ssgsea.norm = T, mx.diff=TRUE, verbose=FALSE, parallel.sz=20)

select_table<-data.frame(sample=colnames(RNA_matrix),
                         gene1=es.dif[1,])

write.csv(select_table,"ssGSEA_result.csv")