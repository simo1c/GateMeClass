library(CATALYST)
library(cowplot)
library(flowCore)
library(diffcyt)
library(scater)
library(SingleCellExperiment)
library(ggpubr)
library(ggrepel)
library(scMerge)
library(mclust)
library(moments)

setwd("flow_cytometry_test")

r1 <- list.files(pattern = "*.fcs")
sce_list <- list()

for(i in 1:length(r1)){
  
  fcs_files <- r1[i]
  patient <- sapply(strsplit(fcs_files,"_"), `[`, 1)
  cell <- sapply(strsplit(fcs_files,"_"), `[`, 2)
  cell <- gsub(".fcs", "", cell)
  
  fs <- read.flowSet(fcs_files, transformation = F, truncate_max_range = F)
  col <- c("FS-A", "SS-A", "FL1-A", "FL10-A", "FL2-A", "FL3-A", "FL4-A", "FL5-A", "FL6-A", "FL7-A", "FL8-A", "FL9-A")
  fs <- fs[, col]
  fs
  
  m0 <- as.matrix(fs@frames[[fcs_files]]@exprs)
  x <- fs[[1]]
  comp_list <- spillover(x)
  comp <- comp_list$`$SPILLOVER`
  
  fs_comp <- compensate(x, comp)
  logic <- estimateLogicle(fs_comp, channels = colnames(fs), m = 4.8)
  fs_comp2 <- transform(fs_comp, logic)
  colnames(fs_comp2) <- colnames(fs)
  colnames(fs_comp2) <- c(colnames(fs)[1:2], "CD57", "CD19","CD45RA","CD8","CD56","CD4","CD27","CD45","CD3", "CD16")
  
  antigen <- colnames(fs_comp2)
  class <- rep("type", length(antigen))
  panel <- data.frame(fcs_colname = antigen, antigen = antigen, marker_class = class)
  md <- data.frame(file_name = fcs_files, sample_id = patient, condition  = patient, patient_id = cell)
  sce <- prepData(fs_comp2, panel = panel, md = md, transform = F, FACS = TRUE)
  rownames(sce) <- antigen
  sce
  
  sce@assays@data$exprs <- sce@assays@data$counts
  
  if(i == 1){
    sce_tot <- sce
  }else{
    sce_tot <- cbind(sce_tot, sce)
  }
}

set.seed(1)
sce_control <- runUMAP(sce_tot, exprs_values = "exprs")

plotDR(sce_control, dr = "UMAP", color_by = "patient_id") +
  theme_classic() +
  xlab("UMAP 1") +
  ylab("UMAP 2") +
  theme(legend.title = element_blank())










