# Script: 05_decoupleR_tf_activity.R
# Purpose: Transcription-factor activity inference with decoupleR
# Workflow: bulk-rna
# Required inputs: differ.csv or the exported DEG result table.
# Main outputs: TF_diff_decoupleR.pdf and decoupleR-derived activity tables.
# Prerequisites: Requires differential expression statistics.


library(decoupleR)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(pheatmap)
library(ggrepel)
resdata <- read.csv("differ.csv", row.names = 1)

counts <- resdata[,9:16]
rownames(counts) <- resdata$Gene

design <- tibble(
  sample = c("moDC_1", "moDC_2", "moDC_3","moDC_4",
             "mrDC_1", "mrDC_2", "mrDC_3","mrDC_4"),
  condition = c("moDC", "moDC", "moDC","moDC",
                "mrDC", "mrDC", "mrDC","mrDC")
)

deg <- resdata[,c("log2FoldChange","stat","pvalue")]
rownames(deg) <- resdata$Gene

contrast_acts <- run_ulm(mat=deg[, 'stat', drop=FALSE], net=net, .source='source', .target='target',
                         .mor='mor', minsize =   5)
contrast_acts

f_contrast_acts <- contrast_acts %>%
  mutate(rnk = NA)

# Rank the scores based on sign (positive or negative)
msk <- f_contrast_acts$score > 0
f_contrast_acts[msk, 'rnk'] <- rank(-f_contrast_acts[msk, 'score'])
f_contrast_acts[!msk, 'rnk'] <- rank(-abs(f_contrast_acts[!msk, 'score']))

colnames(f_contrast_acts)
TF_mregDC <- subset(f_contrast_acts, subset = score >5 & p_value < 0.05)

n_tfs <- 25
tfs <- f_contrast_acts %>%
  arrange(rnk) %>%
  head(n_tfs) %>%
  pull(source)

f_contrast_acts <- f_contrast_acts %>%
  filter(source %in% tfs)

ggplot(f_contrast_acts, aes(x = reorder(source, score), y = score)) +
  geom_bar(aes(fill = score), stat = "identity") +
  scale_fill_gradient2(
    low = "#71c3ee", high = "#e780a6", mid = "white", midpoint = 0
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(face = "bold", size = 12),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
    axis.text.y = element_text(size = 10, face = "bold")
  ) 

ggsave("TF_diff_decoupleR.pdf", width = 9, height = 3.5)