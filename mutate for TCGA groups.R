df <- TCGA_clinical_Data
df <- df%>%
        select(Patient_ID, OS_5y_event,OS_5y_months, PFI_5y_event,PFI_5y_months,
               DSS_5y_event, DSS_5y_months, Gender, Age, Smoking, Alcohol, HPV,
               Subsite, pT, pN, Grading)

df_mutated <- df%>%
        mutate(Age_groups = case_when(Age < 60.91 ~ "Age<Mean",
                                      Age > 60.91 ~ "Age>Mean"))%>%
        mutate(Subsite_groups = case_when( Subsite == "Oropharynx" ~ "OPSCC",
                                           Subsite != "Oropharynx" ~ "non_OPSCC"))%>%
        mutate(T_Status_groups = case_when( pT == "T1" | pT == "T2" ~ "T1-2",
                                            pT == "T3"| pT == "T4a" | pT == "T4b" ~ "T3-4",
                                            pT == "TX" | pT == "NA" ~ "NA"))%>%
        mutate(N_Status_groups = case_when(pN == "N0" ~ "N0",
                                           pN == "NX" ~ "NA",
                                           pN != "N0" | pN != "NX" | pN != "NA" ~ "N+"))%>%
        mutate(Grading_groups = case_when (Grading == "G1" | Grading == "G2" ~ "G1-2",
                                           Grading == "G3" ~ "G3"))


df_mutated <- df%>%
        mutate(Resection = case_when( resection_margin_status == "Negative" | resection_margin_status == "Close"  ~ "R0",
                                           resection_margin_status == "Positive" ~ "R1"))
