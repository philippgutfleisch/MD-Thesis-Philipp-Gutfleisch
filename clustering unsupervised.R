# Clustering unsupervised


library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tidyverse)
library(readxl) 
library(corrplot)
library(cluster)


column_annotations <- TCGA_clinical_Data

matrix <- significant.gsva.scores.pathway.hallmark

row_annotations <- 

#Preparing heatmap matrix
heatmap <- matrix %>% # matrix is a tibble (because you imported it with the tool above) but i want it as a data frame (to set row names)
        as.data.frame(.)
#heatmap <- heatmap[match(GSVA_Pathway_Activity_results$pathways, rownames(heatmap)),]
heatmap <- as.matrix.data.frame(heatmap)
#heatmap <- t(heatmap)
#heatmap <- log(heatmap+1)
heatmap <- scale(heatmap)
#heatmap <- t(heatmap)

#Preparing Annotations
clinical <- column_annotations %>% # It's a tibble but you want a data frame
        as.data.frame(.)
row.names(clinical) <- clinical$Patient_ID
clinical$Patient_ID <- NULL

clinical_annotations = HeatmapAnnotation(df = clinical,
                                         which = "column",
                                         col = list(ICS_risk = c("high" = "red", "low" = "blue"),
                                                    HPV = c("Negative" = "red", "Positive" = "blue")),
                                         annotation_name_side = "right",
                                         na_col = "white")

genes_anno <- row_annotations %>%
        as.data.frame(.)
genes_anno <- genes_anno%>%
        column_to_rownames("Gene_Symbol")

gene_annotation <-  HeatmapAnnotation(df = genes_anno,
                                    which = "row",
                                    col = list("up/down" = c("up" = "red", "down" = "blue" )),
                                    annotation_name_side = "top")

#color range for heatmap values
col_rna = colorRamp2(c(-2, 0, 2), c("blue", "#EEEEEE", "red"), space = "LAB")

#heatmap
HM_unsupervised <- ComplexHeatmap::Heatmap(heatmap,
                                  col = col_rna,
                                  name = "GSVA score",
                                  row_names_gp = grid::gpar(fontsize=7),
                                  
                                  
                                  cluster_rows = T,
                                  cluster_row_slices = T,
                                  clustering_distance_rows = "euclidean",  # "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
                                  clustering_method_rows = "ward.D2", #"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid" 
                                  #row_split = 2,
                                  
                                  
                                  cluster_columns = T,
                                  cluster_column_slices = T,
                                  clustering_distance_columns = "euclidean", # "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski", "pearson", "spearman", "kendall"
                                  clustering_method_columns = "ward.D2",  #"ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"
                                  column_split = 2, 
                                  
                                  #column_title =  c("A", "B"),
                                  column_title_side = "bottom",
                                  show_column_names = F,
                                  show_row_names = F, 
                                  show_row_dend = T,
                                  show_column_dend = T,
                                  use_raster = TRUE,
                                  top_annotation = clinical_annotations)
                                  #left_annotation = gene_annotation)

tiff(filename = "GSVA unsupervised heatmap immunogenic .tiff",width = 32, height = 25,units = "cm",pointsize = 12,res=400)

draw(HM_unsupervised)
dev.off()
