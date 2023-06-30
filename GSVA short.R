#GSVA short
library(tidyverse)
library(GSVA)
library(GSEABase)
library(GSVAdata)
library(WriteXLS)

#load dfs
gene_sets <- TLS_Gensignaturen%>%
        drop_na()
gene_sets <- as.list(gene_sets)
gsva_matrix <- GSE65858_expression_primary_cases%>%
        column_to_rownames("Gene_ID")
gsva_matrix <- as.matrix(gsva_matrix)
gsva_clinical <- GSE65858%>%
        dplyr::select(Patient_ID, ICS_risk)

#check if all genes of genesignature are in matrix
gene_sets$Gene_ID %in% rownames(gsva_matrix)

#run GSVA 
gsva_results <- gsva(gsva_matrix, gene_sets)
gsva_results <- as.data.frame(t(gsva_results))
gsva_results <- gsva_results%>%
        rownames_to_column("Patient_ID")
gsva_results <- gsva_results%>%
        dplyr::rename("GSVA_score" = "Gene_ID")

# merge GSVA and ICS risk for violinplots
gsva_results <- merge(gsva_results,gsva_clinical)
sum(is.na(gsva_results))
gsva_results <- gsva_results%>%
        drop_na()

WriteXLS(gsva_results, "GSVA_GSE65858_Wichmann_TLS.xlsx")

#--------------------------------------------------------------------
# VIOLINPLOTS & T-TEST

library(tidyverse)
library(ggpubr)
library(forcats)
library(psych)
library(WriteXLS)
library(car)
library(lsr)

#violinplot
gsva_violin <- gsva_results

ggplot(data = gsva_violin, aes(x = ICS_risk ,fill = ICS_risk ,
                             y = GSVA_score))+
        geom_violin()+
        geom_boxplot(width=0.3)+ 
        theme_classic()+
        ggtitle("TCGA_GSE65858_Wichmann_TLS")+
        #ylim(-2,2)+
        scale_fill_manual(values = c("red","blue", "yellow"))


#t.test
y <- gsva_violin$ICS_risk
x <- gsva_violin$GSVA_score

describeBy(x,y) # auf sd schauen und schauen ob in etwa gleich. 

#Levene Test:
leveneTest(x,y) # wenn p<0.05 -> Varianz verschieden -> var.equal = F
                # wenn p>0.05 -> Varianz gleich -> var.equal = T

# t-test:
t.test(x ~ y, var.equal = T)

# Cohens D für Effektstärke:
cohensD(x~y)    # cohensD ab 0.2 -> kleiner Effekt
# cohensD ab 0.5 -> mittlerer Effekt
# cohensD ab 0.8 -> großer Effekt

#------------------------------------------------
# spearman correlation

cor.test(df$MS4A1, df$GSVA_score, method = "spearman")

ggplot(df, aes(x=MS4A1_expression, y=GSVA_score))+ 
        geom_point()+
        xlim(0,15)+
        ylim(-1,1)+
        geom_smooth(method=lm)+
        labs(title = "TCGA HPVneg TLS Genesignature MS4A1 Expression")







