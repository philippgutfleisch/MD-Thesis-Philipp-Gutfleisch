# unsupervised Heatmap

BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)
library(circlize)
library(cluster)
library(RColorBrewer)
library(tidyverse)
library(readxl)
library(WriteXLS)



# improt data & files

column_annotations <- DEGs_clinical_ICS_risk_sorted_low_to_high #Sheet: ColAnnotations
View(column_annotations)

matrix <- DEGs_matrix_TPM_sorted_up_to_down_and_low_to_high_risk#Sheet: Heatmap
matrix <- matrix%>%
        column_to_rownames("Gene_Symbol")
View(matrix)  

row_annotations <- DEGs_rowannotations_up_to_down_regulated
View(row_annotations)


#Preparing heatmap matrix
heatmap <- matrix %>% # matrix is a tibble (because you improted it with the tool above) but i want it as a data frame (to set row names)
        as.data.frame(.)
heatmap <- as.matrix.data.frame(heatmap)
heatmap <- t(heatmap)
heatmap <- log(heatmap+1)
heatmap <- scale(heatmap)
heatmap <- t(heatmap)
#heatmap <- scale(heatmap)
#heatmap <- t(heatmap)

#Preparing Annotations
clinical <- column_annotations %>% # It's a tibble but you want a data frame
        as.data.frame(.)
clinical <- clinical%>%
        select(ICS_risk)

clinical_annotations = HeatmapAnnotation(df = clinical,
                                         which = "column",
                                         col = list(
                                                 ICS_risk = c("high"="red", "low" = "blue")
                                         ),
                                         annotation_name_side = "right",
                                         na_col = "white")


genes_anno <- row_annotations %>%
        as.data.frame(.)
genes_anno <- genes_anno%>%
        column_to_rownames("Gene_Symbol")

gene_annotation = HeatmapAnnotation(df = genes_anno,
                                    which = "row",
                                    col = list("up/down" = c("up" = "red", "down" = "blue" )),
                                    annotation_name_side = "top")

#color range for heatmap values
col_fun = colorRamp2(c(-2, 0, 2), c("blue", "white", "red"), space = "LAB", transparency = 0)

#heatmap
HM <- ComplexHeatmap::Heatmap(heatmap,
                              col = col_fun,
                              cluster_rows = F,
                              name = "Expression",
                              row_names_gp = grid::gpar(fontsize=7),
                              cluster_columns = F,
                              column_title_side = "top",
                              show_column_names = F,
                              show_row_names = T, 
                              show_row_dend = F,
                              use_raster = TRUE, 
                              top_annotation = clinical_annotations,
                              left_annotation = gene_annotation) 
tiff(filename = "DEG supervised heatmap.tiff",width = 32,height = 25,units = "cm",pointsize = 12,res=400)
draw(HM)
dev.off()

#Colum/Row Sample Selection
hm <- HM

hm <-  draw(hm)

gl <- column_order(hm)
gln <- list()
gldf <- data.frame()
for(i in 1:length(gl)){
        agg <- colnames(heatmap)[gl[[i]]]
        gln[[i]] <- agg
        agldf <- data.frame(Patient_ID=agg,Kmean=i)
        gldf <- rbind(gldf,agldf)
}
names(gln) <- names(gl)
gln

# export to excel
cluster1 <- data.frame(gln[[1]]) 
write.xlsx(cluster1,"cluster1.xlsx")

cluster2 <- data.frame(gln[[2]])
write.xlsx(cluster2,"cluster2.xlsx")

cluster3 <- data.frame(gln[[3]])
write.xlsx(cluster3,"cluster3.xlsx")

cluster4 <- data.frame(gln[[4]])
write.xlsx(cluster4,"cluster4.xlsx")

cluster5 <- data.frame(gln[[5]])
write.xlsx(cluster5,"cluster5.xlsx")









col = list(
        Gender = c("Male" = "blue", "Female" = "red"),
        Alcohol = c("No" = "blue", "Yes" = "red"),
        Smoking = c("No" = "blue", "Yes" = "red"),
        HPV = c("Positive" = "blue", "Negative" = "red"),
        Subsite_groups = c("OPSCC" = "blue", "non_OPSCC" = "red"),
        Grading_groups = c ("G1-2" = "blue", "G3" = "red"),
        ESR1_related_risk = c("low" = "blue", "high" = "red"),
        OS_5y_event = c("0" = "blue", "1" = "red")
        

