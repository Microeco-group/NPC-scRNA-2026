# Script: 06_survival_analysis.R
# Purpose: Survival analysis from precomputed risk metrics
# Workflow: bulk-rna
# Required inputs: survival_result.csv
# Main outputs: survival_plot.pdf
# Prerequisites: Requires a table with survival time, censoring status, and grouping variables.

rm(list = ls())
library(survival)
library(survminer)

exp <- read.csv("survival_result.csv")   
rownames(exp) <- exp[,1]

res.cut <- surv_cutpoint(exp, 
                                time = "OStime", 
                                event = "Death",
                                variables = "score"  
)
summary(res.cut)
plot(res.cut, "score", palette = "npg") 
res.cat <- surv_categorize(res.cut)
fit1 <- survfit(Surv(OStime, Death) ~ score, data = res.cat)
p1 <- ggsurvplot(fit1,
                 pval = TRUE,  
                 pval.method = T,
                 conf.int = F,  
                 risk.table = "absolute",  
                 risk.table.col="strata",
                 risk.table.y.text.col = T,
                 risk.table.y.text =F,
                 risk.table.pos="out",#in 
                 break.time.by = 6,
                 xlab = "Time in Month", 
                 ylab = "OS",
                 ncensor.plot = F, 
                 legend.title="",
                 legend.labs =c("high", "low"),   
                 palette = c("#e780a6","#71c3ee"), 
                 ggtheme = theme_bw(),
                 tables.theme=theme_cleantable()
)
p1
ggsave("survival_plot.pdf", width = 5, height = 5)