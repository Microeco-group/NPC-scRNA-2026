# Script: 15_cellphonedb_visualization.R
# Purpose: CellPhoneDB result export and visualization
# Workflow: scRNA
# Required inputs: row.txt; colum.txt; CellPhoneDB result tables.
# Main outputs: npc_project_merge.txt; npc_project_filt.txt; CellphoneDB_bubble.pdf
# Prerequisites: Requires CellPhoneDB outputs prepared outside this repository.

rm(list = ls())

library(reshape2)
library(plyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

pvalues <-read.table(paste0("./pvalues.txt"),sep = "\t",row.names = 1, header = T, stringsAsFactors = F,check.names = F)
signification_means <- read.table(paste0("./means.txt"),sep = "\t",row.names = 1, header = T, stringsAsFactors = F,check.names = F)

pvalues_npc <- pvalues[,c(1,11:length(pvalues))]
pvalues_npc_1 <- melt(pvalues_npc,id.var=colnames(pvalues_npc)[1])
colnames(pvalues_npc_1) <- c("interacting_pair","cell_type","pvalues")

signification_means_npc <- signification_means[,c(1,11:length(signification_means))]
signification_means_npc_1 <- melt(signification_means_npc,id.var=colnames(signification_means_npc)[1])
colnames(signification_means_npc_1) <- c("interacting_pair","cell_type","signification_means")

cp <- (pvalues_npc_1)[,1]
cell_type <- pvalues_npc_1[,2]
pvalues <- pvalues_npc_1[,3]
signification_means <- signification_means_npc_1[,3]
cp_DB <- data.frame(cp,cell_type,pvalues,signification_means)
write.table(cp_DB,file =paste0("npc_project_merge.txt") , sep = "\t", row.names = F, quote = F)

cp_DB <- na.omit(cp_DB)
cp_DB <- cp_DB[cp_DB$pvalues < 0.05,]
write.table(cp_DB,file =paste0("npc_project_filt.txt") , sep = "\t", row.names = F, quote = F)

npc_project_merge <- read.table(file = "npc_project_merge.txt",header = T,sep = "\t",stringsAsFactors = F)

dcat_dt <- dcast(data = npc_project_merge,cp ~ cell_type, value.var = "signification_means")
melt_dt <- melt(dcat_dt,value.name = "cp")
colnames(melt_dt)[3] <- "signification_means"
melt_dt$variable <- as.character(melt_dt$variable)
melt_dt$cp_cell_type <- paste0(melt_dt$cp,melt_dt$variable)

dcat_dt1 <- dcast(data = npc_project_merge,cp ~ cell_type, value.var = "pvalues")
melt_dt1 <- melt(dcat_dt1,value.name = "cp")
colnames(melt_dt1)[3] <- "pvalues"
melt_dt1$variable <- as.character(melt_dt1$variable)
melt_dt1$cp_cell_type <- paste0(melt_dt1$cp,melt_dt1$variable)

melt_dt$p_value <- 1
melt_dt$p_value <- mapvalues(melt_dt$cp_cell_type , from = melt_dt1$cp_cell_type, to = melt_dt1$pvalues)
melt_dt$p_value[melt_dt$p_value %in% melt_dt$cp_cell_type] <- 1
melt_dt$p_value <- as.numeric(melt_dt$p_value)
melt_dt$p_value[melt_dt$p_value == 0] <- 0.01
melt_dt$signification_means[melt_dt$signification_means > 2.5] <- 2.5

cp <- read.delim("./row.txt", check.names = FALSE, header = T)
cp <- unique(cp)
cp$interaction <- factor(cp$interaction, levels = cp$interaction)
cell_type <- read.delim("./colum.txt", check.names = FALSE, header = F)
cell_type <- unique(cell_type)

Plot_data <- melt_dt[melt_dt$cp %in% cp$interaction,]
Plot_data <- Plot_data[Plot_data$variable %in% cell_type$V1, ]
Plot_data$cp <- factor(Plot_data$cp, levels = cp$interaction)

cols <- colorRampPalette(brewer.pal(10, "RdBu"))(20) %>% rev()
ggplot(Plot_data, aes(x = variable, y = cp)) + geom_point(aes(size = -log10(p_value),color = signification_means)) +
  scale_colour_gradientn(colours = cols) +
  coord_flip()+
  theme_bw() +
  theme(axis.text.x=element_text(angle = 90,hjust = 1))
ggsave("CellphoneDB_bubble.pdf", height = 4, width = 12.5)
