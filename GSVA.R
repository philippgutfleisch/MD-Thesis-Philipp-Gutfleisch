# GSVA

BiocManager::install("GSVA", force = T)
BiocManager::install("GSVAdata", force = T)

library(tidyverse)
library(GSVA)
library(GSEABase)
library(GSVAdata)

#create list of gene sets of interest
#import ICS_Gene_Sets_up&down 
gene_sets <- GSVA_HNSC_TLS_HPVneg

#-----when doing repeated analysis start here-------------------
# expression matrix, genes in rows, samples in cols
matrix <- TCGA_HNSC_TPM_proteincoding_all_cases
clinical <- TCGA_clinical_Data

gsva_matrix <- as.matrix(matrix%>%
        column_to_rownames("Gene_ID"))
gsva_clinical <- as.data.frame(clinical%>%
        dplyr::select(Patient_ID, ICS_risk))


# check, if all Genes of Genelist are in Matrix
# interpret as follows: if nothing comes: all genes are in matrix
#                       if genes come: those genes are not in matrix
#                               -> check for synonyms!!!
in_up <- ICS_Gene_Sets_up$DEGs_up_regulated %in% rownames(gsva_matrix)
in_up <- cbind(ICS_Gene_Sets_up, in_up)
in_up%>%filter(in_up == F) 

in_down <- ICS_Gene_Sets_down$DEGs_down_regulated %in% rownames(gsva_matrix)
in_down <- cbind(ICS_Gene_Sets_down,in_down)
in_down%>%filter(in_down == F)


#check the synonymes:
sum(rownames(gsva_matrix) == "LOC284254")

#rename Genes, that are different between geva_matrix & up/down
row.names(gsva_matrix)[row.names(gsva_matrix) == "FAM46C"] <- "TENT5C"

#recheck for wrong genenames
in_up <- ICS_Gene_Sets_up$DEGs_up_regulated %in% rownames(gsva_matrix)
in_up <- cbind(ICS_Gene_Sets_up, in_up)
in_up%>%filter(in_up == F) 

in_down <- ICS_Gene_Sets_down$DEGs_down_regulated %in% rownames(gsva_matrix)
in_down <- cbind(ICS_Gene_Sets_down,in_down)
in_down%>%filter(in_down == F)

#safe new gsva_matrix
write.table(gsva_matrix, "new_gsva_matrix.txt", row.names = T, sep = "\t", quote = F)
#-------------------------------------------------------------------
# run GSVA
gsva_results <- gsva(gsva_matrix, gene_sets)
row.names(gsva_results)[row.names(gsva_results) == "DEGs_up_regulated"] <- "GSVA_score_up"
row.names(gsva_results)[row.names(gsva_results) == "DEGs_down_regulated"] <- "GSVA_score_down"
gsva_results <- as.data.frame(t(gsva_results))

# next step only if Patient_ID = rownames
gsva_results <- gsva_results%>%
        rownames_to_column("Patient_ID")

# merge GSVA and ICS risk for violinplots
gsva_results <- merge(gsva_results,gsva_clinical)
sum(is.na(gsva_results))

WriteXLS(gsva_results, "GSVA_HNSC_TLS.xlsx")


#--------------------------------------------------------------------
# VIOLINPLOTS & T-TEST

library(tidyverse)
library(ggpubr)
library(forcats)
library(psych)
library(WriteXLS)
library(car)
library(lsr)

#load data

df_violin <- gsva_results

sum(is.na(df_violin))

df_violin <- na.omit(df_violin)

#------ GSVA score UP regulated ------------------------------------
#violinplot
ggplot(data = df_violin, aes(x = ICS_risk ,fill = ICS_risk ,
                             y = Gene_ID))+
        geom_violin()+
        geom_boxplot(width=0.3)+ 
        theme_classic()+
        ggtitle("HIPO_HNSC_up")+
        #ylim(-2,2)+
        scale_fill_manual(values = c("red","blue", "yellow"))


#t.test
y <- df_violin$ICS_risk
x <- df_violin$Gene_ID

describeBy(x,y) # auf sd schauen und schauen ob in etwa gleich. 
# Levene Test berechnen

#Levene Test:
leveneTest(x,y) # wenn p<0.05 -> Varianz verschieden -> var.equal = F
                # wenn p>0.05 -> Varianz gleich -> var.equal = T

# t-test:
t.test(x ~ y, var.equal = T)

# Cohens D für Effektstärke:
cohensD(x~y)    # cohensD ab 0.2 -> kleiner Effekt
                # cohensD ab 0.5 -> mittlerer Effekt
                # cohensD ab 0.8 -> großer Effekt


#violinplot
ggplot(data = df_violin, aes(x = ICS_risk ,fill = ICS_risk ,
                             y = GSVA_score_down))+
        geom_violin()+
        geom_boxplot(width=0.3)+ 
        theme_classic()+
        ggtitle("HIPO_HNSC_down")+
        #ylim(0,2)+
        scale_fill_manual(values = c("red","blue", "yellow"))


#t.test
y <- df_violin$ICS_risk
x <- df_violin$GSVA_score_down

describeBy(x,y) # auf sd schauen und schauen ob in etwa gleich. 
# Levene Test berechnen
#Levene Test:
leveneTest(x,y) # wenn p<0.05 -> Varianz verschieden -> var.equal = F
                # wenn p>0.05 -> Varianz gleich -> var.equal = T

# t-test:
t.test(x ~ y, var.equal = F)

# Cohens D für Effektstärke:
cohensD(x~y)    # cohensD ab 0.2 -> kleiner Effekt
        # cohensD ab 0.5 -> mittlerer Effekt
        # cohensD ab 0.8 -> großer Effekt





