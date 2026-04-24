# Script: 10_cellchat_communication.R
# Purpose: CellChat communication comparison between groups
# Workflow: scRNA
# Required inputs: cellchat_high.rds; cellchat_low.rds
# Main outputs: cellchat.pdf_1; netVisual_diffInteraction.pdf; netVisual_heatmap.pdf; netVisual_bubble_possible_all.pdf
# Prerequisites: Requires group-specific CellChat objects prepared upstream.

rm(list = ls())
library(CellChat)
library(patchwork)
cellchat_high <- readRDS("cellchat_high.rds")
cellchat_low <- readRDS("cellchat_low.rds")

object.list <- list(Low = cellchat_low, High = cellchat_high)
cellchat <- mergeCellChat(object.list,  add.names = names(object.list))
cellchat

col_cellchat <- c("#eed785","#d4b595","#efe0e7","#eebed6")
names(col_cellchat) <- c("c39_cDC1_CLEC9A",
                         "c40_cDC2_CD1C",
                         "c41_mregDC_LAMP3",
                         "c42_pDC_LILRA4")
a=cellchat_low@net$weight
b=cellchat_high@net$weight
res= data.frame(edit_subtype=colnames(a),
                strength_difference=b['c07_CD8_Tex_PDCD1',]/a['c07_CD8_Tex_PDCD1',])
res <- subset(res, edit_subtype %in% c("c39_cDC1_CLEC9A",
                                       "c40_cDC2_CD1C",
                                       "c41_mregDC_LAMP3",
                                       "c42_pDC_LILRA4"))
data <- res[order(res$strength_difference, decreasing = TRUE), ]
data$edit_subtype <- factor(data$edit_subtype, levels = rev(unique(data$edit_subtype)))

g1 <- ggplot(data, aes(x = strength_difference, y = reorder(edit_subtype, strength_difference))) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_segment(
    aes(x = 1, xend = strength_difference, y = edit_subtype, yend = edit_subtype, color = edit_subtype),
    linewidth = 1.5) +
  geom_point(
    aes(size = 2, color = edit_subtype),
    size = 4  ) +
  scale_color_manual(values = col_cellchat) +
  theme_bw() + 
  labs(y = '', x = 'FoldChange of Strength Difference') +
  theme(
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(hjust = 0.5)) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1)))

a=cellchat_low@net$count
b=cellchat_high@net$count
res= data.frame(edit_subtype=colnames(a),
                strength_difference=b['c07_CD8_Tex_PDCD1',]/a['c07_CD8_Tex_PDCD1',])
res <- subset(res, edit_subtype %in% c("c39_cDC1_CLEC9A",
                                       "c40_cDC2_CD1C",
                                       "c41_mregDC_LAMP3",
                                       "c42_pDC_LILRA4"))
data <- res[order(res$strength_difference, decreasing = TRUE), ]
data$edit_subtype <- factor(data$edit_subtype, levels = rev(unique(data$edit_subtype)))

g2 <- ggplot(data, aes(x = strength_difference, y = reorder(edit_subtype, strength_difference))) + 
  geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
  geom_segment(
    aes(x = 1, xend = strength_difference, y = edit_subtype, yend = edit_subtype, color = edit_subtype),
    linewidth = 1.5) +
  geom_point(
    aes(size = 2, color = edit_subtype),
    size = 4  
  ) +
  scale_color_manual(values = col_cellchat) +
  theme_bw() + 
  labs(y = '', x = 'FoldChange of Number Difference') +
  theme(
    axis.text = element_text(size = 12),
    axis.title.x = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12),
    axis.text.y = element_text(hjust = 0.5)) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.1)))

g1+g2+plot_layout(guides = "collect", nrow = 1)
ggsave("cellchat.pdf_1", width = 9, height = 3)


pdf("netVisual_diffInteraction.pdf", width = 10, height = 8)
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_diffInteraction(cellchat, 
                          weight.scale = TRUE, 
                          comparison = c(1, 2),
                          color.use = col_cellchat,
                          color.edge = c("#fd9aa0", "#6dccfd"))
netVisual_diffInteraction(cellchat, 
                          weight.scale = TRUE, 
                          measure = "weight", 
                          comparison = c(1, 2),
                          color.use = col_cellchat,
                          color.edge = c("#fd9aa0", "#6dccfd"))
dev.off()

pdf("netVisual_heatmap.pdf", width = 9, height = 5)
gg1 <- netVisual_heatmap(cellchat, 
                         comparison = c(1, 2),
                         color.use = col_cellchat,
                         color.heatmap = c("#6dccfd", "#fd9aa0"))
gg2 <- netVisual_heatmap(cellchat, 
                         measure = "weight", 
                         comparison = c(1, 2),
                         color.use = col_cellchat,
                         color.heatmap = c("#6dccfd", "#fd9aa0"))
gg1 + gg2
dev.off()


num.link <- sapply(object.list, function(x) {
  rowSums(x@net$count) + colSums(x@net$count) - diag(x@net$count)
})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(object.list)) {
  a <- netAnalysis_computeCentrality(object.list[[i]]) 
  gg[[i]] <- netAnalysis_signalingRole_scatter(a, 
                                               title = names(object.list)[i], 
                                               color.use = col_cellchat,
                                               weight.MinMax = weight.MinMax)+
    theme(aspect.ratio = 1)+
    scale_x_continuous(limits = c(0, 0.5)) + 
    scale_y_continuous(limits = c(0, 0.6))
}
patchwork::wrap_plots(plots = gg)
ggsave("netAnalysis_signalingRole_scatter.pdf", width = 8, height = 4)


netVisual_bubble(cellchat, sources.use = 4, targets.use = c(5:11), comparison = c(1, 2), angle.x = 45)
netVisual_bubble(cellchat, sources.use = c(2:5), targets.use = 1, comparison = c(1, 2), angle.x = 45)
display <- subset(aa, subset = annotation == "Cell-Cell Contact")
gg1 <- netVisual_bubble(
  cellchat, 
  sources.use = c("c39_cDC1_CLEC9A"  ,   "c40_cDC2_CD1C"  ,  "c41_mregDC_LAMP3"  ,  "c42_pDC_LILRA4" ), 
  targets.use = c("c07_CD8_Tex_PDCD1"), 
  #signaling = unique(display$pathway_name),
  comparison = c(1, 2), 
  max.dataset = 2, 
  title.name = "Increased signaling in High", 
  angle.x = 45, 
  remove.isolate = TRUE,
  color.text = c("#6dccfd", "#fd9aa0")
)
ggsave("netVisual_bubble_possible_all.pdf",width = 6,height = 9)