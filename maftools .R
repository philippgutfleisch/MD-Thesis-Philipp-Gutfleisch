# somatic mutation analysis with maftools

#if (!require("BiocManager", quietly = TRUE))
        #install.packages("BiocManager")
#BiocManager::install(version = "3.17")
#BiocManager::install("maftools")

library(maftools)
library(TCGAbiolinks)
library(tidyverse)
library(WriteXLS)

# download maf file with reference genome hg38
query <- GDCquery(
        project = "TCGA-HNSC", 
        data.category = "Simple Nucleotide Variation", 
        access = "open", 
        legacy = FALSE, 
        data.type = "Masked Somatic Mutation", 
        workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")

GDCdownload(query)
maf_TCGA <- GDCprepare(query)
maf_TCGA <- maftools::read.maf(maf_TCGA,
                               rmFlags = T, 
                               isTCGA = TRUE)

#-------------------------------------------------------------------------------
#subset in high and low risk patients
#create character variable of Patient_IDs of high and low risk group
clinical <- TCGA_clinical_Data 
high_risk <- clinical%>%
        filter(ICS_risk == "high")
high_risk_maf <- subsetMaf(maf = maf_TCGA, tsb = high_risk[["Patient_ID"]])

low_risk <- clinical%>%
        filter(ICS_risk == "low")
low_risk_maf <- subsetMaf(maf = maf_TCGA, tsb = low_risk[["Patient_ID"]])

#Analysis:
#define subgroup which should be analysed (fill in high or low risk as risk maf to do 
#seperate analysis of the subgroups)
risk_maf <- low_risk_maf

#visualisation
#general visualisation of maf dataset
plotmafSummary(maf = risk_maf, addStat = "median")
#oncoplot
oncoplot(maf = risk_maf, legend_height = 8, altered = F)
#transition & transversion plot
TCGA_titv <- titv(maf = risk_maf, plot = F, useSyn = T)
plotTiTv(res = TCGA_titv)
#compare to TCGA Cohorts
tcgaCompare(risk_maf, cohortName = "risk_MAF")
#variant allel frequency plot
plotVaf(maf = risk_maf)

#Somatic Analysis
#somatic interaction analysis (which genes mutate or don`t mutate with other genes)
somaticInteractions(maf = risk_maf,top = 20, 
                    pvalue = c(0.05,0.1),
                    fontSize = 0.4, 
                    nShiftSymbols = 5)
#oncodrive
maf.sig = oncodrive(maf = risk_maf, minMut = 5, pvalMethod = 'zscore')
head(maf.sig)
plotOncodrive(res = maf.sig, fdrCutOff = 0.1, useFraction = TRUE)

#comparison of high vs. low risk patients
mutation_comparison <- mafCompare(m1 = low_risk_maf,
                                  m1Name = "ICS_low_risk",
                                  m2 = high_risk_maf,
                                  m2Name = "ICS_high_risk",
                                  minMut = 5,
                                  useCNV = F)
mutation_comparison
WriteXLS(mutation_comparison$results, "somatic mutation comparison high vs low risk all cases.xlsx")

#forest plot of the differently mutated genes
forestPlot(mafCompareRes = mutation_comparison,
           pVal = 0.05,
           fdr = 0.1, # fdr = p value adjustment method, use fdr = 0.1 if no genes are significant
           geneFontSize = 0.6,
           titleSize = 1.2,
           lineWidth = 1)

#select differently mutated genes based on adj.p.value, use cutoff adjPval<0.1 if no genes are significant
differently_mutated_genes <- mutation_comparison$results$Hugo_Symbol[mutation_comparison$results$pval< 0.05]
differently_mutated_genes <- differently_mutated_genes[1:10]

coOncoplot(m1 = low_risk_maf,
           m1Name = "ICS_low_risk",
           m2 = high_risk_maf,
           m2Name = "ICS_high_risk",
           genes = differently_mutated_genes,
           geneNamefont = 0.5,
           legend_height = 6)


coBarplot(m1 = low_risk_maf,
          m1Name = "ICS_low_risk",
          m2 = high_risk_maf,
          m2Name = "ICS_high_risk",
          genes = differently_mutated_genes,
          showLegend = T, 
          geneSize = 0.5)

# drug gene interactions
dgi_high_risk <- drugInteractions(high_risk_maf)
dgi_low_risk <- drugInteractions(low_risk_maf)

# oncogenic pathways
pathways_high_risk <- OncogenicPathways(high_risk_maf) 
pathways_low_risk <- OncogenicPathways(low_risk_maf)

#plot oncogenic pathway to further look into the genes mutated in pathways of interest
PlotOncogenicPathways(maf = high_risk_maf, pathways = "NOTCH")

