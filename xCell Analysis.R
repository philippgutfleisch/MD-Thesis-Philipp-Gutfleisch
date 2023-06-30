#xCell Analysis

# load xCell package
#devtools::install_github('dviraran/xCell')
library(xCell)
library(tidyverse)
library(xlsx)


df <- df_without_NA_duplicates

# Gene Symbols have to be rownames
df <- df%>%
        column_to_rownames("Gene_ID")

xCell_result <- xCellAnalysis(df,
                        rnaseq = F # RNAseq Data?
                        )

xCell_result <- as.data.frame(xCell_result)

write.xlsx(xCell_result, "GSE65858 xCell results.xlsx")
