# kmeans clustering

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)
library(tidyverse)
library(readxl) # for importing Excel files
library(corrplot)
library(cluster)
library(WriteXLS)

View(column_annotations)

View(matrix)

row_annotations <- ICS_risk_DEGs_Philipp_Counts

#Preparing heatmap matrix
heatmap <- matrix %>% # matrix is a tibble (because you imported it with the tool above) but i want it as a data frame (to set row names)
        as.data.frame(.)
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

#col_annotations
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
col_rna = colorRamp2(c(-2,0,2), c("blue", "#EEEEEE", "red"), space = "LAB", transparency = 0)

#heatmap
set.seed(123) #to get reproducible kmeans results
HM_kmeans <- Heatmap(heatmap,
                     name = "GSVA score",
                     column_km = 2,
                     column_km_repeats = 1000,
                     cluster_rows = T,
                     row_km = 2,
                     row_km_repeats = 1000,
                     col = col_rna,
                     row_names_gp = grid::gpar(fontsize=6),
                     column_names_gp = grid::gpar(fontsize=6),
                     column_title_side = "bottom",
                     show_column_names = F,
                     show_row_names = F, 
                     show_row_dend = T,
                     show_column_dend = T,
                     use_raster = TRUE,
                     top_annotation = clinical_annotations)
                     #left_annotation = gene_annotation)

tiff(filename = "GSVA kmeans heatmap immunogenic .tiff",width = 32,height = 25,units = "cm",pointsize = 12,res=400)

draw(HM_kmeans)
dev.off()



library(factoextra)
library(cluster)

df <- matrix
sum(is.na(df))
df <- t(df)
#df <- log(df+1)
df <- scale(df)
df <- t(df)



#elbow plot
fviz_nbclust(df, kmeans, method = "wss")

#calculate gap statistic based on number of clusters
gap_stat <- clusGap(df,
                    FUN = kmeans,
                    set.seed(123),
                    K.max = 10,
                    B = 500)

#plot number of clusters vs. gap statistic
fviz_gap_stat(gap_stat)

#make this example reproducible
set.seed(123)

#perform k-means clustering 
km <- kmeans(df, centers = 2)

#view results
km

#plot results of final k-means model
fviz_cluster(km, data = df)                     

