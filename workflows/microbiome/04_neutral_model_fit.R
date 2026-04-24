# Script: 04_neutral_model_fit.R
# Purpose: Neutral community model fitting
# Workflow: microbiome
# Required inputs: physeq.end.rds
# Main outputs: neutral_fit.pdf
# Prerequisites: Requires a phyloseq object with abundance and taxonomy data.

rm(list = ls())

library(ape)
library(phyloseq)
library(MicEco)
library(ggplot2)
library(dplyr)
library(ggimage)

physeq.great <- readRDS("physeq.end.rds")
physeq.great.genus <- tax_glom(physeq.great, "g")

otu <- as.data.frame(t(physeq.great.genus@otu_table@.Data))
res = neutral.fit(otu)
m = res[[1]][1]
r2 = res[[1]][3]
out = res[[2]]

p1 = ggplot()+
  geom_line(data = out,aes(x=log(p),y=freq.pred),size = 1.2,linetype = 1)+
  geom_line(data = out,aes(x=log(p),y=Lower),size = 1.2,linetype = 2)+
  geom_line(data = out,aes(x=log(p),y=Upper),size = 1.2,linetype = 2)+
  xlab("log10(mean relative abundance)")+ylab("Occurrence frequency")+
  coord_fixed(ratio = 10)

out2 = mutate(out, group=NA)
out2$group[out[,2]<out[,4]]="Below" 
out2$group[out[,2]>out[,5]]="Above" 
out2$group[(out[,2]>=out[,4])&(out[,2]<=out[,5])]="In"
mycols<-c("#e780a6","#71c3ee","grey")

p2 <- p1 + geom_point(data = out2,aes(x=log(p),y=freq,color = group, shape = differ),size = 5,alpha = 0.75)+
  scale_shape_manual(values = c(1,16))+
  scale_colour_manual(values = mycols)+
  annotate("text",x=-11.5,y=0.95,label=paste("m = ",round(m,3),sep=''),size=6)+
  annotate("text",x=-11.5,y=1,label=paste("R2 = ",round(r2,3),sep=''),size=6)

plot_theme = theme(panel.background=element_blank(),
                   panel.grid=element_blank(),
                   axis.line.x=element_line(size=.5, colour="black"),
                   axis.line.y=element_line(size=.5, colour="black"),
                   axis.ticks=element_line(color="black"),
                   axis.text=element_text(color="black", size=15),
                   legend.position="right",   ##right
                   legend.background=element_blank(),
                   legend.key=element_blank(),
                   legend.text= element_text(size=15),
                   text=element_text(family="sans", size=15),
                   aspect.ratio = 1)

p3 = p2 + plot_theme

Below = nrow(out2[out2[,6]== "Below",])
In = nrow(out2[out2[,6]== "In",])
Above = nrow(out2[out2[,6]== "Above",])
group <- c('Above','Below','Med')
nums <- c(Above,Below,In)
data_summary <- data.frame(group = group, nums = nums)
data_summary$percentage <- round(data_summary$nums / sum(data_summary$nums) * 100 , 1)
data_summary$label <- paste(data_summary$group,paste(data_summary$percentage, "%", sep = ''))

p4 <- ggplot(data_summary, aes(x = "", y = nums, fill = group))+
  geom_bar(stat = "identity", width = 1)+
  coord_polar(theta = "y")+
  scale_fill_manual(values = c("#e780a6","#71c3ee","grey"), labels = data_summary$label)+
  theme_void()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text = element_text(size=12))
p4
p_final <- p3 + geom_subview(subview = p4, x=-11, y=0.5,w=4.5,h=4.5)
p_final
ggsave("neutral_fit.pdf", height = 8, width = 8)
