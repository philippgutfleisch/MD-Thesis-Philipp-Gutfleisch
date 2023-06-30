library(tidyverse)


#BiocManager::install("illuminaHumanv4.db")
library("illuminaHumanv4.db")


probeID <- GSE65858.RAW.DATA%>%
        dplyr::select(ID_REF) # function select is masked from illuminaHumanv4.db
probeID_vector <- probeID$ID_REF
  
Gene_ID_df <- data.frame(Gene=unlist(mget(x = probeID_vector,envir = illuminaHumanv4SYMBOL)))
Gene_ID_df <- Gene_ID_df%>%
        rownames_to_column("ID_REF")%>%
        rename("Gene_ID" = "Gene")

GSE <- merge(Gene_ID_df, GSE65858.RAW.DATA, all = T)





write.table(GSE,"GSE65858.txt", sep = "\t", quote = F, row.names = F)


# to use dplyr package command select, detach illumina db
detach("package:illuminaHumanv4.db", unload = TRUE)
detach("package:org.Hs.eg.db", unload = TRUE)
detach("package:AnnotationDbi", unload = TRUE) # kommt warning wird trotzdem detached

#prepare df for xCell analysis
GSE_xCell <- GSE%>%
        select(!ID_REF)
