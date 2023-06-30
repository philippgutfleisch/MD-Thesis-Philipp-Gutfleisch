install.packages(c("tidyverse",
                 "WriteXLS",
                 "circlize",
                 "cluster",
                 "RColorBrewer",
                 "readxl",
                 "corrplot",
                 "org.Hs.eg.db",
                 "ConsensusClusterPlus",
                 "survminer",
                 "survival",
                 "DESeq2",
                 "edgeR",
                 "VennDiagram",
                 "forestplot",
                 "AnnotationDbi",
                 "clusterProfiler",
                 "GSVA",
                 "GSEABase",
                 "GSVAdata",
                 "ggpubr",
                 "forcats",
                 "psych",
                 "car",
                 "lsr",
                 "maxstat",
                 "corrplot",
                 "glmnet",
                 "magick",
                 "BiocManager",
                 "devtools"
                 ))

BiocManager::install("TCGAbiolinks", force = T)
BiocManager::install(version = "3.17")
BiocManager::install("maftools")
BiocManager::install("GSVA", force = T)
BiocManager::install("GSVAdata")
BiocManager::install("limma")
BiocManager::install("ComplexHeatmap")

