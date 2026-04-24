# Script: 11_scmetabolism_analysis.R
# Purpose: Metabolic pathway scoring for cDC cells
# Workflow: scRNA
# Required inputs: cDC.rds; access to biomaRt/Ensembl for gene mapping.
# Main outputs: cDC_scMetabolism_type.rds and downstream visualizations.
# Prerequisites: Requires an annotated cDC subset and online biomaRt access unless cached locally.

rm(list = ls())
library(Seurat)
library(dplyr)
library(monocle)
options(stringsAsFactors=FALSE)
library(reticulate)
library(SingleR)
library(scMetabolism)
library(ggplot2)
library(rsvd)
library(gridExtra)
library(httr)
library(biomaRt)
library(stringr)
library(conflicted)
library(dplyr)
library(tidyr)
library(tibble)
library(circlize)
library(RColorBrewer)

cDC <- readRDS("cDC.rds")
Idents(cDC) <- cDC$celltype3
sce <- cDC

ensembl <- with_config(config(httr::timeout(60)), {
  useMart("ensembl", dataset = "hsapiens_gene_ensembl")
})
conversion <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),
                    filters = "ensembl_gene_id", 
                    values = rownames(sce),
                    mart = ensembl)
sce_sub <- sce[rownames(sce) %in% conversion$ensembl_gene_id[!duplicated(conversion$external_gene_name)],]
conflicts_prefer(biomaRt::select)
RenameGenesSeurat <- function(obj, newnames) { 
  print("Run this before integration. It only changes obj@assays$RNA@counts, @data, and @scale.data.") 
  RNA <- obj@assays$RNA 
  if (nrow(RNA) == length(newnames)) { 
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]] <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]] <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]] <- newnames 
  } else { 
    print("Unequal gene sets: nrow(RNA) != nrow(newnames)") } 
  obj@assays$RNA <- RNA 
  return(obj)} 
sce_sub <- RenameGenesSeurat(obj = sce_sub, newnames = conversion$external_gene_name[!duplicated(conversion$external_gene_name)])

results <- list()
for(var in unique(Idents(sce_sub))) {
  var_name <- paste0('sce_', str_replace_all(var, "\\s+", "_")) 
  sce_subset <- subset(sce_sub, idents = var) %>% 
    sc.metabolism.Seurat(method = "VISION", imputation = FALSE, metabolism.type = "KEGG") 
  results[[var_name]] <- sce_subset
}
saveRDS(results, file = 'cDC_scMetabolism_type.rds')

conflicts_prefer(dplyr::filter(), dplyr::select())
sce_sub <- sce_sub %>% sc.metabolism.Seurat(method = "VISION", imputation = FALSE, metabolism.type = "KEGG", ncores = 25)
metabolism.matrix <- sce_sub@assays$METABOLISM$score
metabolism.matrix <- t(metabolism.matrix)
metabolism.matrix <- as.data.frame(metabolism.matrix)
metadata <- sce_sub@meta.data
metabolism.matrix$celltype3 <- metadata$celltype3

library(dplyr)
aver_result <- metabolism.matrix %>%
  group_by(celltype3) %>%
  summarise(across(everything(), mean, na.rm = TRUE))
aver_result <- aver_result[,2:ncol(aver_result)]
aver_result <- t(aver_result)
aver_result_df <- as.data.frame(aver_result)

filtered_data_cDC1 <- subset(aver_result_df, V1 > V2 & V1 > V3)
filtered_data_cDC1 <- filtered_data_cDC1[order(filtered_data_cDC1$V1, decreasing = T),]
colnames(filtered_data_cDC1) <- c("cDC1","cDC2","mregDC")
filtered_data_cDC1 <- filtered_data_cDC1[1:6,]

filtered_data_cDC2 <- subset(aver_result_df, V2 > V3 & V2 > V1)
filtered_data_cDC2 <- filtered_data_cDC2[order(filtered_data_cDC2$V2, decreasing = T),]
colnames(filtered_data_cDC2) <- c("cDC1","cDC2","mregDC")
filtered_data_cDC2 <- filtered_data_cDC2[1:6,]

filtered_data_mregDC <- subset(aver_result_df, V3 > V2 & V3 > V1)
filtered_data_mregDC <- filtered_data_mregDC[order(filtered_data_mregDC$V3, decreasing = T),]
colnames(filtered_data_mregDC) <- c("cDC1","cDC2","mregDC")
filtered_data_mregDC <- filtered_data_mregDC[1:6,]

filtered_data_merge <- rbind(filtered_data_cDC1,filtered_data_cDC2,filtered_data_mregDC)

library(colorRamp2)
library(ComplexHeatmap)
htdf <- t(scale(t(filtered_data_merge),scale = T,center = T))
col_fun = colorRamp2::colorRamp2(c(-2, 0, 2), c("#6DCCFD", "white", "#F484AE"))
col_celltype3 <- c("#eed785","#d4b595","#efe0e7")
names(col_celltype3) <- colnames(filtered_data_merge)
column_ha = HeatmapAnnotation(cluster = colnames(htdf),
                              col = list(cluster = col_celltype3))
pdf("Heatmap_DC_scMetabolism.pdf", width = 6.5, height = 6)
Heatmap(htdf,
        name = "Metabolism Activity (Z-score)",
        cluster_columns = F,cluster_rows = F,
        row_names_side = 'left',
        border = T,
        rect_gp = gpar(col = "white", lwd = 1),
        column_names_side = 'top',
        column_names_rot = 45,
        top_annotation = column_ha,
        col = col_fun,
        width = 3, 
        height = 3 
)
dev.off()


filtered_data <- subset(aver_result_df, V3 > V2 & V3 > V1)
filtered_data <- filtered_data[order(filtered_data$V3, decreasing = T),]
colnames(filtered_data) <- c("cDC1","cDC2","mregDC")
filtered_data <- filtered_data[1:6,]



data_long <- filtered_data %>%
  rownames_to_column("Pathway") %>%
  pivot_longer(cols = -Pathway, 
               names_to = "CellType", 
               values_to = "Activity")
data_long$CellType <- factor(data_long$CellType,levels = c("cDC1","cDC2","mregDC"))

df <- data_long
colnames(df) <- c("Category","Group","Value")
df$Category <- factor(df$Category, levels = unique(df$Category))
group_name <- unique(df$Group)
group_colors <- c("#eed785","#d4b595","#efe0e7")
names(group_colors) <- group_name

cat_name <- unique(df$Category)
cat_num <-length(cat_name)
cat_colors <- brewer.pal(n = cat_num, name ="Set2")
names(cat_colors)<- cat_name

pdf("plot_circos_test.pdf",width=7,height=7)
circos.clear()
circos.par(
  start.degree =90,
  gap.degree =15,
  track.margin =c(0,0.12),
  cell.padding =c(0,0,0,0)
)
circos.initialize(
  factors = df$Category,
  xlim =c(0,1)
)

circos.track(
  factors = df$Category,
  ylim =c(0,1),
  track.height =0.1,
  bg.col = cat_colors,
  panel.fun =function(x,y){
    circos.text(
      CELL_META$xcenter,
      CELL_META$ylim[2]+0.8,
      CELL_META$sector.index,
      facing ="bending.inside",
      cex =1.1,
      adj =c(0.5,0)
    )
  }
)

circos.trackPlotRegion(
  factors = df$Category,
  ylim = c(0, max(df$Value) * 1.15), 
  track.height = 0.5,
  bg.col = adjustcolor(cat_colors, alpha.f = 0),
  panel.fun = function(x, y) {
    sector_data <- df[df$Category == CELL_META$sector.index, ]
    group_num <- length(unique(sector_data$Group))
    current_max <- max(sector_data$Value)
    max_value <- current_max * 1 
    at <- pretty(c(0, 0.75), n = 3)  
    at <- at[at > 0 & at < 0.75] 
    for(a in at) {
      circos.lines(
        CELL_META$xlim,
        c(a, a),  
        lty = 3,
        col = "grey50",
        lwd = 0.8
      )
      circos.text(
        CELL_META$cell.xlim[1] - mm_h(1.5),
        a,  
        labels = a,
        facing = "clockwise",
        adj = c(0.5, 0.5), 
        cex = 0.8,
        col = "grey50"
      )
    }
    n_bars <- nrow(sector_data)
    bar_width <- 0.75 / group_num  
    gap_width <- 0.25 / group_num  
    for(i in 1:n_bars) {
      x_left <- (i - 1) * (bar_width + gap_width) + gap_width / 2
      x_right <- x_left + bar_width
      bar_height <- sector_data$Value[i]  
      circos.rect(
        xleft = x_left,
        ybottom = 0,
        xright = x_right,
        ytop = bar_height,  
        col = group_colors[sector_data$Group[i]], 
        border = "white",
        lwd = 1.2
      )
    }
  }
)
dev.off()
