# Script: 13_meta_program_analysis.R
# Purpose: Meta-program aggregation from NMF factors
# Workflow: scRNA
# Required inputs: Genes_nmf_w_basis.rds
# Main outputs: Robust_NMF_MP.pdf; MP_list.rds
# Prerequisites: Run after the NMF workflow or equivalent basis generation.

#Gavish, A. et al. Hallmarks of transcriptional intratumour heterogeneity across a thousand tumours. Nature 618, 598-606 (2023). 
#https://doi.org/10.1038/s41586-023-06130-4

rm(list = ls())
library(reshape2)
library(NMF)
library(ggplot2)
library(scales)
source("custom_magma.R")
source("robust_nmf_programs.R")   

Genes_nmf_w_basis <- readRDS("Genes_nmf_w_basis.rds")

intra_min_parameter <- 30
intra_max_parameter <- 10
inter_min_parameter <- 10

nmf_programs          <- lapply(Genes_nmf_w_basis, function(x) apply(x, 2, function(y) names(sort(y, decreasing = T))[1:50]))
nmf_filter_ccle       <- robust_nmf_programs(nmf_programs, intra_min = intra_min_parameter, intra_max = intra_max_parameter, inter_filter=T, inter_min = inter_min_parameter)  
nmf_programs          <- lapply(nmf_programs, function(x) x[, is.element(colnames(x), nmf_filter_ccle),drop=F])
nmf_programs          <- do.call(cbind, nmf_programs)

nmf_intersect        <- apply(nmf_programs , 2, function(x) apply(nmf_programs , 2, function(y) length(intersect(x,y)))) 
nmf_intersect_hc     <- hclust(as.dist(50-nmf_intersect), method="average") 
nmf_intersect_hc     <- reorder(as.dendrogram(nmf_intersect_hc), colMeans(nmf_intersect))
nmf_intersect        <- nmf_intersect[order.dendrogram(nmf_intersect_hc), order.dendrogram(nmf_intersect_hc)]

Min_intersect_initial <- 10    
Min_intersect_cluster <- 5    
Min_group_size        <- 2    

Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE)

Cluster_list              <- list()  
MP_list                   <- list()
k                         <- 1
Curr_cluster              <- c()

nmf_intersect_original    <- nmf_intersect

while (Sorted_intersection[1]>Min_group_size) {  
  Curr_cluster <- c(Curr_cluster , names(Sorted_intersection[1]))
  Genes_MP                    <- nmf_programs[,names(Sorted_intersection[1])] 
  nmf_programs                <- nmf_programs[,-match(names(Sorted_intersection[1]) , colnames(nmf_programs))] 
  Intersection_with_Genes_MP  <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE)
  NMF_history                 <- Genes_MP 
  
  while ( Intersection_with_Genes_MP[1] >= Min_intersect_cluster) {  
    Curr_cluster  <- c(Curr_cluster , names(Intersection_with_Genes_MP)[1])
    Genes_MP_temp   <- sort(table(c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]])), decreasing = TRUE)  
    Genes_at_border <- Genes_MP_temp[which(Genes_MP_temp == Genes_MP_temp[50])]  
    
    if (length(Genes_at_border)>1){
      Genes_curr_NMF_score <- c()
      for (i in Curr_cluster) {
        curr_study           <- paste( strsplit(i , "[.]")[[1]][1 : which(strsplit(i , "[.]")[[1]] == "RDS")]   , collapse = "."  )
        Q                    <- Genes_nmf_w_basis[[curr_study]][ match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]])))[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))]   ,i] 
        names(Q)             <- names(Genes_at_border[!is.na(match(names(Genes_at_border),toupper(rownames(Genes_nmf_w_basis[[curr_study]]))))])
        Genes_curr_NMF_score <- c(Genes_curr_NMF_score,  Q )
      }
      Genes_curr_NMF_score_sort <- sort(Genes_curr_NMF_score , decreasing = TRUE)
      Genes_curr_NMF_score_sort <- Genes_curr_NMF_score_sort[unique(names(Genes_curr_NMF_score_sort))]   
      Genes_MP_temp             <- c(names(Genes_MP_temp[which(Genes_MP_temp > Genes_MP_temp[50])]) , names(Genes_curr_NMF_score_sort))
    } else {
      Genes_MP_temp <- names(Genes_MP_temp)[1:50] 
    }
    
    NMF_history     <- c(NMF_history , nmf_programs[,names(Intersection_with_Genes_MP)[1]]) 
    Genes_MP        <- Genes_MP_temp[1:50]
    nmf_programs    <- nmf_programs[,-match(names(Intersection_with_Genes_MP)[1] , colnames(nmf_programs))]  
    Intersection_with_Genes_MP <- sort(apply(nmf_programs, 2, function(x) length(intersect(Genes_MP,x))) , decreasing = TRUE) 
  }
  
  Cluster_list[[paste0("Cluster_",k)]] <- Curr_cluster
  MP_list[[paste0("MP_",k)]]           <- Genes_MP
  k <- k+1
  nmf_intersect             <- nmf_intersect[-match(Curr_cluster,rownames(nmf_intersect) ) , -match(Curr_cluster,colnames(nmf_intersect) ) ] 
  Sorted_intersection       <-  sort(apply(nmf_intersect , 2, function(x) (length(which(x>=Min_intersect_initial))-1)  ) , decreasing = TRUE) 
  Curr_cluster <- c()
  print(dim(nmf_intersect)[2])
}

inds_sorted <- c()

for (j in 1:length(Cluster_list)){
  inds_sorted <- c(inds_sorted , match(Cluster_list[[j]] , colnames(nmf_intersect_original)))
}
inds_new <- c(inds_sorted   ,   which(is.na( match(1:dim(nmf_intersect_original)[2],inds_sorted)))) 
nmf_intersect_meltI_NEW <- reshape2::melt(nmf_intersect_original[inds_new,inds_new]) 

ggplot(data = nmf_intersect_meltI_NEW, aes(x=Var1, y=Var2, fill=100*value/(100-value), color=100*value/(100-value))) + 
  geom_tile() + 
  scale_color_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)") +                                
  scale_fill_gradient2(limits=c(2,25), low=custom_magma[1:111],  mid =custom_magma[112:222], high = custom_magma[223:333], midpoint = 13.5, oob=squish, name="Similarity\n(Jaccard index)")  +
  theme( axis.ticks = element_blank(), panel.border = element_rect(fill=F), panel.background = element_blank(),  axis.line = element_blank(), axis.text = element_text(size = 11), axis.title = element_text(size = 12), legend.title = element_text(size=11), legend.text = element_text(size = 10), legend.text.align = 0.5, legend.justification = "bottom") + 
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank()) + 
  theme(axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  guides(fill = guide_colourbar(barheight = 4, barwidth = 1))+
  coord_fixed(ratio = 1)
ggsave("Robust_NMF_MP.pdf", height = 6, width = 6)
MP_list <-  do.call(cbind, MP_list)
saveRDS(MP_list,"MP_list.rds")