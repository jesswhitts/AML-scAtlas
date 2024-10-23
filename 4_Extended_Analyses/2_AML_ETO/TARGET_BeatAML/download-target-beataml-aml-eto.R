library(readxl)
library(dplyr)
library(TCGAbiolinks)
library(SummarizedExperiment)

### Download TARGET Data ###
query <- GDCquery(
  project = c("TARGET-AML"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification")
GDCdownload(query = query)

# Remove duplicate cases
query.2=query
tmp=query.2$results[[1]]
tmp=tmp[which(!duplicated(tmp$cases)),]
query.2$results[[1]]=tmp
data <- GDCprepare(query.2)
# Save SE
saveRDS(data, "data/TARGET_se.rds")
# Save count mtx
rownames(data) <- data@rowRanges$gene_name
counts <- as.data.frame(assays(data)$unstranded)
write.csv(counts, "data/TARGET_counts.tsv", sep='\t')

# Select AML-ETO from TARGET cohort
TARGET_AML <- read_excel("data/TARGET_AML_ClinicalData_Discovery_20230720.xlsx")
AML_ETO_TARGET <- subset(TARGET_AML, `t(8;21)`=="Yes")
AML_ETO_TARGET$`Age at Diagnosis Calculated Years` <- AML_ETO_TARGET$`Age at Diagnosis in Days`/365.25
AML_ETO_TARGET$`Age at Diagnosis Years` <- floor(AML_ETO_TARGET$`Age at Diagnosis Calculated Years`) 
write.csv(AML_ETO_TARGET, "data/AML_ETO_TARGET.clinical.csv")

# Clean TARGET Metadata
patientIDs <- AML_ETO_TARGET$`TARGET USI`
TARGET_GEX <- readRDS("data/TARGET_se.rds")
patients <- unique(TARGET_GEX@colData@listData[["patient"]])
patient_intersect <- intersect(patientIDs, patients)
subset_TARGET_GEX <- subset(TARGET_GEX, select=TARGET_GEX[["patient"]] %in% patient_intersect)
rownames(subset_TARGET_GEX) <- subset_TARGET_GEX@rowRanges$gene_name
TARGET_counts <- as.data.frame(assays(subset_TARGET_GEX)$unstranded)
TARGET_AML_ETO_obs <- as.data.frame(subset_TARGET_GEX@colData@listData[c("sample", "specimen_type","race","ethnicity","gender",
                                                                         "barcode","patient","sample_submitter_id",
                                                                         "age_at_index","days_to_birth")])
write.csv(TARGET_AML_ETO_obs, "data/TARGET_AML_ETO_obs.csv")


### Download BeatAML Data ###
query <- GDCquery(
  project = c("BEATAML1.0-COHORT"),
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification")
GDCdownload(query = query)
data <- GDCprepare(query = query)
# Save SE
saveRDS(data, "data/BEATAML_se.rds")
# Save count mtx
rownames(data) <- data@rowRanges$gene_name
counts <- as.data.frame(assays(data)$unstranded)
write.csv(counts, "data/BEATAML_counts.tsv", sep='\t')

# Select AML-ETO from BeatAML cohort
BEATAML_GEX <- readRDS("data/BEATAML_se.rds")
BEATAML_GEX@colData@listData[["age_at_diagnosis"]] <- as.integer(BEATAML_GEX@colData@listData[["age_at_diagnosis"]])
BEATAML_GEX@colData@listData[["age_at_diagnosis_calculated_years"]]  <- BEATAML_GEX@colData@listData[["age_at_diagnosis"]]/365.25
BEATAML_GEX@colData@listData[["age_at_diagnosis_years"]]  <- floor(BEATAML_GEX@colData@listData[["age_at_diagnosis_calculated_years"]]) 
subset_BEATAML <- subset(BEATAML_GEX, select=
                           BEATAML_GEX@colData@listData[["primary_diagnosis"]]==
                           "Acute myeloid leukemia with t(8;21)(q22;q22); RUNX1-RUNX1T1")
BEATAML_obs <- as.data.frame(subset_BEATAML@colData@listData[c("gender","race","ethnicity","primary_diagnosis","eln_risk_classification",
                                                               "age_at_diagnosis","age_at_diagnosis_calculated_years","age_at_diagnosis_years",
                                                               "specimen_type","sample_type","tumor_descriptor","sample")])
write.csv(BEATAML_obs, "data/BEATAML_obs.csv")

# Merge BeatAML and TARGET counts 
rownames(subset_BEATAML) <- subset_BEATAML@rowRanges$gene_name
BEATAML_counts <- as.data.frame(assays(subset_BEATAML)$unstranded)
counts <- cbind(TARGET_counts, BEATAML_counts)
write.csv(counts, "data/AML_ETO_sample_counts.csv")