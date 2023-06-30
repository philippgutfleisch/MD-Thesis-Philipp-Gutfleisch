#Subset Mutation Data

mutation <- HNSC.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg

#Herausfiltern von Primärtumoren
Patient_ID <- mutation$Sample
Patient_ID <- substr(Patient_ID,1,16)
mutation <- cbind(Patient_ID,mutation)

tumors <- mutation$Patient_ID
tumor_sample_A <- tumors %in% grep("01A",tumors,value = T)
tumor_sample_B <- tumors %in% grep("01B", tumors,value = T)

primary_tumors_A <- subset(mutation, tumor_sample_A)
primary_tumors_B <- subset(mutation, tumor_sample_B)

primary_tumors <- rbind(primary_tumors_A,primary_tumors_B)

# Patient_ID auf 12 spaces begrenzen, damit df mit anderen dfs gemerged werden kann
Patient_ID <- primary_tumors$Patient_ID
Patient_ID <- substr(Patient_ID,1,12)
primary_tumors <- primary_tumors%>%
        select(!Patient_ID)
primary_tumors <- cbind(Patient_ID,primary_tumors)


ICS_risk <- TCGA_clinical_Data%>%
        select(Patient_ID, ICS_risk)
# zsm.fügen von cnv data und risk, um später high von low risk Gruppen trennen zu können.
ICS_risk_mutation <- merge(ICS_risk, primary_tumors)

write.table(ICS_risk_mutation,"TCGA_HNSC_CNV.txt", sep = "\t", quote = F, col.names = T, row.names = F)


# set cutoffs ( gain: >0,2 ; loss: <-0,2)

high_risk <- CNV_Analysis_high_risk
as.numeric(high_risk$Segment_Mean)

gains <- ifelse(high_risk$Segment_Mean > 0.2,"GAIN",0)
losses <- ifelse(high_risk$Segment_Mean < (-0.2), "LOSS",0)

high_risk <- cbind(high_risk,gains,losses)
