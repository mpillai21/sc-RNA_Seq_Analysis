#!/usr/bin/ 

##This script allows us to perform QC on multiple samples at the same time without merging the samples

#Load the packages
library(tidyverse)
library(magrittr)
library(janitor)
library(Seurat)
library(Matrix)
library(DropletUtils)
library(stringdist)
library(igraph)
library(ggseqlogo)
library(DoubletFinder)
library(ggrepel)
library(celldex)
library(SeuratDisk)
library(readr)

#vignettes used:
#sc-RNA: https://satijalab.org/seurat/articles/pbmc3k_tutorial.html


#Set working directory
setwd("/opt/localdata/msomada/gdm")

#Load the samples
sample1.data <-  Read10X(data.dir = "//outs/filtered_feature_bc_matrix") #path to the cellranger output
OO1 <- CreateSeuratObject(counts = sample1.data, project = "OO1")
sample2.data <-  Read10X(data.dir = "//outs/filtered_feature_bc_matrix")
OO2 <- CreateSeuratObject(counts = sample2.data, project = "OO2")

#Make a list of the objects 
data = c(OO1, OO2)
dataset = c("OO1", "OO2")
data <- setNames(data, dataset)
saveRDS(data, "./Objects/data.rds")

#Quality Control
#Calculate the mitopercent and ribopercent 
data = lapply(data, FUN= function(x){
  # association between UMI and genes
  x@meta.data  %<>% mutate(log10genePerUMI = log10(nFeature_RNA)/log10(nCount_RNA))
  # mitopercent
  x[['mitopercent']] = PercentageFeatureSet(object = x, pattern = "^MT-")#change the pattern as per the gene name format for the organism being analyzed.
  # proportion of ribosome genes in each cell
  x[['ribopercent']] =PercentageFeatureSet(object = x, pattern = "^RP[SL][[:digit:]]")
  return(x)
  
})
# umi counts
umi_count = lapply(data , function(x) {
  RidgePlot(x,
            features = "nCount_RNA",
            group.by = "orig.ident") +
    guides(fill = "none") +
    scale_fill_viridis_d() + 
    xlim(c(0,15000))
}
)
p = cowplot::plot_grid(umi_count[[1]], umi_count[[2]], ncol = 2)
p

p = cowplot::plot_grid(umi_count[[1]],  
                       umi_count[[2]], ncol = 2)
p
ggsave("./umi_plot.pdf", p, units = "in", height = 8, width = 11, dpi = 72)



# percent mito
mito_sample = mapply(FUN =function(x,y){
  data = x@meta.data
  plot = data %>% ggplot(aes(x='', y=mitopercent)) + geom_violin() + 
    geom_jitter(width = 0.3, alpha = 0.3)+ylim(c(0,60)) +
    ggtitle(y)
  return(plot)
},
data, as.list(dataset),
SIMPLIFY=F
)
mito <- egg::ggarrange(plots = mito_sample, ncol = 2)
ggsave("./mito_plot.pdf", mito, units = "in", height = 5, width = 8, dpi = 72)

#n_gene plot
ngene_sample = mapply(FUN = function(x,y){
  data = x@meta.data
  plot = data %>% ggplot(aes(x= '', y= nFeature_RNA)) + geom_violin() +
    geom_jitter(width =0.2, alpha = 0.3) +
    xlab(y)
  return(plot)
},
data, as.list(dataset),
SIMPLIFY=F
)
ngene = egg::ggarrange(plots =ngene_sample , ncol =3)
ggsave("./ngene_plot.pdf", ngene, units = "in", height = 5, width = 8, dpi = 72)

#Count vs Feature Scatter Plot
pointLoggeneUMI_sample = mapply(FUN = function(x,y){
  data = x@meta.data
  plot = data %>% ggplot(aes(nCount_RNA, nFeature_RNA)) + geom_point() +
    ggtitle(y)
  return(plot)
},
data, as.list(dataset),
SIMPLIFY=F
)
scatter = egg::ggarrange(plots =pointLoggeneUMI_sample , ncol =3)
ggsave("./pointLoggeneUMI_sample_plot.pdf",scatter, units = "in", height = 5, width = 8, dpi = 72)


# add log10gene_per_umi
data = lapply(data, function(x){
  new = x@meta.data %>% mutate(log10gene_per_umi = log10(nFeature_RNA)/log10(nCount_RNA)) %>% 
    dplyr::select(log10gene_per_umi)
  rownames(new)  = rownames(x@meta.data)
  x = AddMetaData(x, new)
  return(x)
  
})

# filter out by mito, minimum gens, log linear 
data = lapply(data, function(x){
  x = subset(x, subset =( (mitopercent < 5)
      (nFeature_RNA >= 100) &
      (nCount_RNA >= 100) &
      (log10gene_per_umi >0.85)
  )
  ) 
}
)
saveRDS(data, "./data_after_qc.rds")
