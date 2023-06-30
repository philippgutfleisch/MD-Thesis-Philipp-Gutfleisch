library(tidyverse)
library(plyr)

#load datasets

df <- HS_CPTAC_HNSCC_RNAseq_RSEM_UQ_log2_Tumor.cct

#drop NAs in Gene Symbols

df <- df%>%
        drop_na(Gene_ID)
sum(is.na.data.frame(df$Gene_ID))

#search for duplicates
sum(duplicated(df$Gene_ID))

# filter out singles and duplicates
#filters singles
df_single_genes <- df%>%
        group_by(Gene_ID)%>%
        filter(n() == 1)%>%
        ungroup()
sum(duplicated(df_single_genes$Gene_ID))

#filters duplicates and computes mean of duplicates
df_duplicated_genes <- df%>%
        group_by(Gene_ID)%>%
        filter(n()>1)%>%
        summarise_all(mean)%>%
        ungroup()

sum(duplicated(df_duplicated_genes$Gene_ID))

# create new df with no NAs & no duplicates
df_without_NA_duplicates <- rbind(df_single_genes,df_duplicated_genes)

#recheck for NAs & duplicates
sum(duplicated(df_without_NA_duplicates$Gene_ID))
sum(is.na.data.frame(df))

#modify df how you want it
#df_without_NA_duplicates <- df_without_NA_duplicates%>%
        #column_to_rownames("Gene_ID")

#export as .txt for xCell/Cibersortx
write.table(ILM_to_ID, "E-MTAB-8588 ILM to Gene ID.txt" , quote = F, sep = "\t")    


