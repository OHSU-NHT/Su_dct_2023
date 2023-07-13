
library(Seurat)
library(dplyr)
library(here)
library(ggplot2)
library(RColorBrewer)
library(tibble)
library(geneName)
library(ggrepel)
library(ggpmisc)

# Download KPMP snRNAseq dataset (DOI:10.48698/yyvc-ak78) from https://atlas.kpmp.org/repository/

# Load the dataset and add participants' information as metadata --------------------------------------------------------------------
kpmp <- LoadH5Seurat(here("WashU-UCSD_HuBMAP_KPMP-Biopsy_10X-R_12032021.h5Seurat"))

#Tissue source
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(tissue_source = dplyr::case_when(
    specimen_id == '18-142' ~ "KPMP Pilot",
    specimen_id == '18-162' ~ "KPMP Pilot",
    specimen_id == '18-312' ~ "KPMP Pilot",
    specimen_id == '3490' ~ "Tissue Interrogation Site",
    specimen_id == '3499' ~ "Tissue Interrogation Site",
    specimen_id == '3504' ~ "Tissue Interrogation Site",
    specimen_id == '3535' ~ "Tissue Interrogation Site",
    specimen_id == '3593' ~ "Tissue Interrogation Site",
    specimen_id == '3613' ~ "Tissue Interrogation Site",
    specimen_id == 'KRP446' ~ "Tissue Interrogation Site",
    specimen_id == 'KRP460' ~ "Tissue Interrogation Site",
    specimen_id == 'KRP461' ~ "Tissue Interrogation Site",
    specimen_id == 'KRP462' ~ "Tissue Interrogation Site",
    specimen_id == '30-10034' ~ "KPMP Recruitment Site",
    specimen_id == '32-10003' ~ "KPMP Recruitment Site",
    specimen_id == '32-10034' ~ "KPMP Recruitment Site",
    specimen_id == '32-2' ~ "KPMP Recruitment Site",
    specimen_id == '33-10005' ~ "KPMP Recruitment Site",
    specimen_id == '33-10006' ~ "KPMP Recruitment Site",
    specimen_id == '29-10006' ~ "KPMP Recruitment Site",
    specimen_id == '29-10008' ~ "KPMP Recruitment Site",
    specimen_id == '29-10010' ~ "KPMP Recruitment Site",
    specimen_id == '29-10012' ~ "KPMP Recruitment Site",
    specimen_id == '29-10013' ~ "KPMP Recruitment Site",
    specimen_id == '31-10000' ~ "KPMP Recruitment Site",
    specimen_id == '31-10001' ~ "KPMP Recruitment Site",
    specimen_id == '31-10006' ~ "KPMP Recruitment Site",
    specimen_id == '31-10013' ~ "KPMP Recruitment Site",
    specimen_id == '31-10035' ~ "KPMP Recruitment Site"
  ))
#Protocol
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(protocol = dplyr::case_when(
    specimen_id == '18-142' ~ "KPMP Pilot 1 Protocol",
    specimen_id == '18-162' ~ "KPMP Pilot 1 Protocol",
    specimen_id == '18-312' ~ "KPMP Pilot 1 Protocol",
    specimen_id == '3490' ~ "Non-KPMP Protocol",
    specimen_id == '3499' ~ "Non-KPMP Protocol",
    specimen_id == '3504' ~ "Non-KPMP Protocol",
    specimen_id == '3535' ~ "Non-KPMP Protocol",
    specimen_id == '3593' ~ "Non-KPMP Protocol",
    specimen_id == '3613' ~ "Non-KPMP Protocol",
    specimen_id == 'KRP446' ~ "Non-KPMP Protocol",
    specimen_id == 'KRP460' ~ "Non-KPMP Protocol",
    specimen_id == 'KRP461' ~ "Non-KPMP Protocol",
    specimen_id == 'KRP462' ~ "Non-KPMP Protocol",
    specimen_id == '30-10034' ~ "KPMP Main Protocol",
    specimen_id == '32-10003' ~ "KPMP Main Protocol",
    specimen_id == '32-10034' ~ "KPMP Main Protocol",
    specimen_id == '32-2' ~ "KPMP Main Protocol",
    specimen_id == '33-10005' ~ "KPMP Main Protocol",
    specimen_id == '33-10006' ~ "KPMP Main Protocol",
    specimen_id == '29-10006' ~ "KPMP Main Protocol",
    specimen_id == '29-10008' ~ "KPMP Main Protocol",
    specimen_id == '29-10010' ~ "KPMP Main Protocol",
    specimen_id == '29-10012' ~ "KPMP Main Protocol",
    specimen_id == '29-10013' ~ "KPMP Main Protocol",
    specimen_id == '31-10000' ~ "KPMP Main Protocol",
    specimen_id == '31-10001' ~ "KPMP Main Protocol",
    specimen_id == '31-10006' ~ "KPMP Main Protocol",
    specimen_id == '31-10013' ~ "KPMP Main Protocol",
    specimen_id == '31-10035' ~ "KPMP Main Protocol"
  ))
#Sample type
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(sample_type = dplyr::case_when(
    specimen_id == '18-142' ~ "Tumor Nephrectomy",
    specimen_id == '18-162' ~ "Tumor Nephrectomy",
    specimen_id == '18-312' ~ "Tumor Nephrectomy",
    specimen_id == '3490' ~ "Percutaneous Needle Biopsy",
    specimen_id == '3499' ~ "Deceased Donor Nephrectomy",
    specimen_id == '3504' ~ "Partial Tumor Nephrectomy",
    specimen_id == '3535' ~ "Total Tumor Nephrectomy",
    specimen_id == '3593' ~ "Deceased Donor Nephrectomy",
    specimen_id == '3613' ~ "Deceased Donor Nephrectomy",
    specimen_id == 'KRP446' ~ "Intra-operative Biopsy",
    specimen_id == 'KRP460' ~ "Intra-operative Biopsy",
    specimen_id == 'KRP461' ~ "Intra-operative Biopsy",
    specimen_id == 'KRP462' ~ "Intra-operative Biopsy",
    specimen_id == '30-10034' ~ "Percutaneous Needle Biopsy",
    specimen_id == '32-10003' ~ "Percutaneous Needle Biopsy",
    specimen_id == '32-10034' ~ "Percutaneous Needle Biopsy",
    specimen_id == '32-2' ~ "Percutaneous Needle Biopsy",
    specimen_id == '33-10005' ~ "Percutaneous Needle Biopsy",
    specimen_id == '33-10006' ~ "Percutaneous Needle Biopsy",
    specimen_id == '29-10006' ~ "Percutaneous Needle Biopsy",
    specimen_id == '29-10008' ~ "Percutaneous Needle Biopsy",
    specimen_id == '29-10010' ~ "Percutaneous Needle Biopsy",
    specimen_id == '29-10012' ~ "Percutaneous Needle Biopsy",
    specimen_id == '29-10013' ~ "Percutaneous Needle Biopsy",
    specimen_id == '31-10000' ~ "Percutaneous Needle Biopsy",
    specimen_id == '31-10001' ~ "Percutaneous Needle Biopsy",
    specimen_id == '31-10006' ~ "Percutaneous Needle Biopsy",
    specimen_id == '31-10013' ~ "Percutaneous Needle Biopsy",
    specimen_id == '31-10035' ~ "Percutaneous Needle Biopsy"
  ))
#Sex
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(sex = dplyr::case_when(
    specimen_id == '18-142' ~ "Female",
    specimen_id == '18-162' ~ "Female",
    specimen_id == '18-312' ~ "Male",
    specimen_id == '3490' ~ "Female",
    specimen_id == '3499' ~ "Female",
    specimen_id == '3504' ~ "Female",
    specimen_id == '3535' ~ "Male",
    specimen_id == '3593' ~ "Male",
    specimen_id == '3613' ~ "Female",
    specimen_id == 'KRP446' ~ "Male",
    specimen_id == 'KRP460' ~ "Female",
    specimen_id == 'KRP461' ~ "Male",
    specimen_id == 'KRP462' ~ "Female",
    specimen_id == '30-10034' ~ "Female",
    specimen_id == '32-10003' ~ "Male",
    specimen_id == '32-10034' ~ "Male",
    specimen_id == '32-2' ~ "Male",
    specimen_id == '33-10005' ~ "Male",
    specimen_id == '33-10006' ~ "Male",
    specimen_id == '29-10006' ~ "Female",
    specimen_id == '29-10008' ~ "Female",
    specimen_id == '29-10010' ~ "Male",
    specimen_id == '29-10012' ~ "Female",
    specimen_id == '29-10013' ~ "Female",
    specimen_id == '31-10000' ~ "Male",
    specimen_id == '31-10001' ~ "Male",
    specimen_id == '31-10006' ~ "Male",
    specimen_id == '31-10013' ~ "Male",
    specimen_id == '31-10035' ~ "Female"
  ))
#Age
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(age = dplyr::case_when(
    specimen_id == '18-142' ~ "60-69",
    specimen_id == '18-162' ~ "60-69",
    specimen_id == '18-312' ~ "70-79",
    specimen_id == '3490' ~ "",
    specimen_id == '3499' ~ "",
    specimen_id == '3504' ~ "",
    specimen_id == '3535' ~ "",
    specimen_id == '3593' ~ "",
    specimen_id == '3613' ~ "",
    specimen_id == 'KRP446' ~ "",
    specimen_id == 'KRP460' ~ "",
    specimen_id == 'KRP461' ~ "",
    specimen_id == 'KRP462' ~ "",
    specimen_id == '30-10034' ~ "70-79",
    specimen_id == '32-10003' ~ "30-39",
    specimen_id == '32-10034' ~ "60-69",
    specimen_id == '32-2' ~ "60-69",
    specimen_id == '33-10005' ~ "60-69",
    specimen_id == '33-10006' ~ "20-29",
    specimen_id == '29-10006' ~ "60-69",
    specimen_id == '29-10008' ~ "60-69",
    specimen_id == '29-10010' ~ "70-79",
    specimen_id == '29-10012' ~ "50-59",
    specimen_id == '29-10013' ~ "60-69",
    specimen_id == '31-10000' ~ "50-59",
    specimen_id == '31-10001' ~ "70-79",
    specimen_id == '31-10006' ~ "60-69",
    specimen_id == '31-10013' ~ "60-69",
    specimen_id == '31-10035' ~ "70-79"
  ))
#Race
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(race = dplyr::case_when(
    specimen_id == '18-142' ~ "",
    specimen_id == '18-162' ~ "",
    specimen_id == '18-312' ~ "",
    specimen_id == '3490' ~ "White",
    specimen_id == '3499' ~ "White",
    specimen_id == '3504' ~ "White",
    specimen_id == '3535' ~ "White",
    specimen_id == '3593' ~ "White",
    specimen_id == '3613' ~ "White",
    specimen_id == 'KRP446' ~ "White",
    specimen_id == 'KRP460' ~ "White",
    specimen_id == 'KRP461' ~ "White",
    specimen_id == 'KRP462' ~ "White",
    specimen_id == '30-10034' ~ "Black or African-American",
    specimen_id == '32-10003' ~ "White",
    specimen_id == '32-10034' ~ "White",
    specimen_id == '32-2' ~ "White",
    specimen_id == '33-10005' ~ "White",
    specimen_id == '33-10006' ~ "White",
    specimen_id == '29-10006' ~ "Black or African-American",
    specimen_id == '29-10008' ~ "White",
    specimen_id == '29-10010' ~ "Asian",
    specimen_id == '29-10012' ~ "Black or African-American",
    specimen_id == '29-10013' ~ "Black or African-American",
    specimen_id == '31-10000' ~ "White",
    specimen_id == '31-10001' ~ "White",
    specimen_id == '31-10006' ~ "White",
    specimen_id == '31-10013' ~ "White",
    specimen_id == '31-10035' ~ "White"
  ))
#KDIGO Stage
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(KDIGO_stage = dplyr::case_when(
    specimen_id == '18-142' ~ "",
    specimen_id == '18-162' ~ "",
    specimen_id == '18-312' ~ "",
    specimen_id == '3490' ~ "",
    specimen_id == '3499' ~ "",
    specimen_id == '3504' ~ "",
    specimen_id == '3535' ~ "",
    specimen_id == '3593' ~ "",
    specimen_id == '3613' ~ "",
    specimen_id == 'KRP446' ~ "",
    specimen_id == 'KRP460' ~ "",
    specimen_id == 'KRP461' ~ "",
    specimen_id == 'KRP462' ~ "",
    specimen_id == '30-10034' ~ "Stage 2",
    specimen_id == '32-10003' ~ "Stage 2",
    specimen_id == '32-10034' ~ "Stage 2",
    specimen_id == '32-2' ~ "Stage 3",
    specimen_id == '33-10005' ~ "Stage 3",
    specimen_id == '33-10006' ~ "Stage 3",
    specimen_id == '29-10006' ~ "",
    specimen_id == '29-10008' ~ "",
    specimen_id == '29-10010' ~ "",
    specimen_id == '29-10012' ~ "",
    specimen_id == '29-10013' ~ "",
    specimen_id == '31-10000' ~ "",
    specimen_id == '31-10001' ~ "",
    specimen_id == '31-10006' ~ "",
    specimen_id == '31-10013' ~ "",
    specimen_id == '31-10035' ~ ""
  ))
#Baseline eGFR
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(baseline_eGFR = dplyr::case_when(
    specimen_id == '18-142' ~ "",
    specimen_id == '18-162' ~ "",
    specimen_id == '18-312' ~ "",
    specimen_id == '3490' ~ "",
    specimen_id == '3499' ~ "",
    specimen_id == '3504' ~ "",
    specimen_id == '3535' ~ "",
    specimen_id == '3593' ~ "",
    specimen_id == '3613' ~ "",
    specimen_id == 'KRP446' ~ "",
    specimen_id == 'KRP460' ~ "",
    specimen_id == 'KRP461' ~ "",
    specimen_id == 'KRP462' ~ "",
    specimen_id == '30-10034' ~ "30-39",
    specimen_id == '32-10003' ~ "60-69",
    specimen_id == '32-10034' ~ "50-59",
    specimen_id == '32-2' ~ "50-59",
    specimen_id == '33-10005' ~ "50-59",
    specimen_id == '33-10006' ~ "60-69",
    specimen_id == '29-10006' ~ "70-79",
    specimen_id == '29-10008' ~ "40-49",
    specimen_id == '29-10010' ~ "40-49",
    specimen_id == '29-10012' ~ "30-39",
    specimen_id == '29-10013' ~ "100-109",
    specimen_id == '31-10000' ~ "20-29",
    specimen_id == '31-10001' ~ "40-49",
    specimen_id == '31-10006' ~ "50-59",
    specimen_id == '31-10013' ~ "50-59",
    specimen_id == '31-10035' ~ "60-69"
  ))
#Diabetes history
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(diabetes_history = dplyr::case_when(
    specimen_id == '18-142' ~ "",
    specimen_id == '18-162' ~ "",
    specimen_id == '18-312' ~ "",
    specimen_id == '3490' ~ "",
    specimen_id == '3499' ~ "",
    specimen_id == '3504' ~ "",
    specimen_id == '3535' ~ "",
    specimen_id == '3593' ~ "",
    specimen_id == '3613' ~ "",
    specimen_id == 'KRP446' ~ "",
    specimen_id == 'KRP460' ~ "",
    specimen_id == 'KRP461' ~ "",
    specimen_id == 'KRP462' ~ "",
    specimen_id == '30-10034' ~ "Yes",
    specimen_id == '32-10003' ~ "No",
    specimen_id == '32-10034' ~ "Yes",
    specimen_id == '32-2' ~ "Yes",
    specimen_id == '33-10005' ~ "No",
    specimen_id == '33-10006' ~ "No",
    specimen_id == '29-10006' ~ "Yes",
    specimen_id == '29-10008' ~ "No",
    specimen_id == '29-10010' ~ "Yes",
    specimen_id == '29-10012' ~ "Yes",
    specimen_id == '29-10013' ~ "Yes",
    specimen_id == '31-10000' ~ "No",
    specimen_id == '31-10001' ~ "Yes",
    specimen_id == '31-10006' ~ "No",
    specimen_id == '31-10013' ~ "Yes",
    specimen_id == '31-10035' ~ "Yes"
  ))
#Diabetes duration
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(diabetes_duration_years = dplyr::case_when(
    specimen_id == '18-142' ~ "",
    specimen_id == '18-162' ~ "",
    specimen_id == '18-312' ~ "",
    specimen_id == '3490' ~ "",
    specimen_id == '3499' ~ "",
    specimen_id == '3504' ~ "",
    specimen_id == '3535' ~ "",
    specimen_id == '3593' ~ "",
    specimen_id == '3613' ~ "",
    specimen_id == 'KRP446' ~ "",
    specimen_id == 'KRP460' ~ "",
    specimen_id == 'KRP461' ~ "",
    specimen_id == 'KRP462' ~ "",
    specimen_id == '30-10034' ~ "",
    specimen_id == '32-10003' ~ "",
    specimen_id == '32-10034' ~ "10-14",
    specimen_id == '32-2' ~ "0-4",
    specimen_id == '33-10005' ~ "",
    specimen_id == '33-10006' ~ "",
    specimen_id == '29-10006' ~ "5-9",
    specimen_id == '29-10008' ~ "",
    specimen_id == '29-10010' ~ "10-14",
    specimen_id == '29-10012' ~ "10-14",
    specimen_id == '29-10013' ~ "5-9",
    specimen_id == '31-10000' ~ "",
    specimen_id == '31-10001' ~ "25-29",
    specimen_id == '31-10006' ~ "",
    specimen_id == '31-10013' ~ "0-4",
    specimen_id == '31-10035' ~ "10-14"
  ))
#Hypertension history
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(hypertension_history = dplyr::case_when(
    specimen_id == '18-142' ~ "",
    specimen_id == '18-162' ~ "",
    specimen_id == '18-312' ~ "",
    specimen_id == '3490' ~ "",
    specimen_id == '3499' ~ "",
    specimen_id == '3504' ~ "",
    specimen_id == '3535' ~ "",
    specimen_id == '3593' ~ "",
    specimen_id == '3613' ~ "",
    specimen_id == 'KRP446' ~ "",
    specimen_id == 'KRP460' ~ "",
    specimen_id == 'KRP461' ~ "",
    specimen_id == 'KRP462' ~ "",
    specimen_id == '30-10034' ~ "Yes",
    specimen_id == '32-10003' ~ "No",
    specimen_id == '32-10034' ~ "No",
    specimen_id == '32-2' ~ "Yes",
    specimen_id == '33-10005' ~ "Yes",
    specimen_id == '33-10006' ~ "Yes",
    specimen_id == '29-10006' ~ "Yes",
    specimen_id == '29-10008' ~ "Yes",
    specimen_id == '29-10010' ~ "Yes",
    specimen_id == '29-10012' ~ "Yes",
    specimen_id == '29-10013' ~ "Yes",
    specimen_id == '31-10000' ~ "Yes",
    specimen_id == '31-10001' ~ "Yes",
    specimen_id == '31-10006' ~ "Yes",
    specimen_id == '31-10013' ~ "Yes",
    specimen_id == '31-10035' ~ "Yes"
  ))
#Hypertension duration
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(hypertension_duration_years = dplyr::case_when(
    specimen_id == '18-142' ~ "",
    specimen_id == '18-162' ~ "",
    specimen_id == '18-312' ~ "",
    specimen_id == '3490' ~ "",
    specimen_id == '3499' ~ "",
    specimen_id == '3504' ~ "",
    specimen_id == '3535' ~ "",
    specimen_id == '3593' ~ "",
    specimen_id == '3613' ~ "",
    specimen_id == 'KRP446' ~ "",
    specimen_id == 'KRP460' ~ "",
    specimen_id == 'KRP461' ~ "",
    specimen_id == 'KRP462' ~ "",
    specimen_id == '30-10034' ~ "",
    specimen_id == '32-10003' ~ "",
    specimen_id == '32-10034' ~ "",
    specimen_id == '32-2' ~ "0-4",
    specimen_id == '33-10005' ~ "",
    specimen_id == '33-10006' ~ "",
    specimen_id == '29-10006' ~ "30-34",
    specimen_id == '29-10008' ~ "0-4",
    specimen_id == '29-10010' ~ "0-4",
    specimen_id == '29-10012' ~ "10-14",
    specimen_id == '29-10013' ~ "20-24",
    specimen_id == '31-10000' ~ "15-19",
    specimen_id == '31-10001' ~ "25-29",
    specimen_id == '31-10006' ~ "5-9",
    specimen_id == '31-10013' ~ "10-14",
    specimen_id == '31-10035' ~ "10-14"
  ))
#On RAS blockade
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(RAAS_blockade = dplyr::case_when(
    specimen_id == '18-142' ~ "",
    specimen_id == '18-162' ~ "",
    specimen_id == '18-312' ~ "",
    specimen_id == '3490' ~ "",
    specimen_id == '3499' ~ "",
    specimen_id == '3504' ~ "",
    specimen_id == '3535' ~ "",
    specimen_id == '3593' ~ "",
    specimen_id == '3613' ~ "",
    specimen_id == 'KRP446' ~ "",
    specimen_id == 'KRP460' ~ "",
    specimen_id == 'KRP461' ~ "",
    specimen_id == 'KRP462' ~ "",
    specimen_id == '30-10034' ~ "No",
    specimen_id == '32-10003' ~ "No",
    specimen_id == '32-10034' ~ "No",
    specimen_id == '32-2' ~ "No",
    specimen_id == '33-10005' ~ "Yes",
    specimen_id == '33-10006' ~ "Yes",
    specimen_id == '29-10006' ~ "Yes",
    specimen_id == '29-10008' ~ "No",
    specimen_id == '29-10010' ~ "No",
    specimen_id == '29-10012' ~ "No",
    specimen_id == '29-10013' ~ "Yes",
    specimen_id == '31-10000' ~ "Yes",
    specimen_id == '31-10001' ~ "Yes",
    specimen_id == '31-10006' ~ "No",
    specimen_id == '31-10013' ~ "Yes",
    specimen_id == '31-10035' ~ "Yes"
  ))
#Proteinuria
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(proteinuria = dplyr::case_when(
    specimen_id == '18-142' ~ "",
    specimen_id == '18-162' ~ "",
    specimen_id == '18-312' ~ "",
    specimen_id == '3490' ~ "",
    specimen_id == '3499' ~ "",
    specimen_id == '3504' ~ "",
    specimen_id == '3535' ~ "",
    specimen_id == '3593' ~ "",
    specimen_id == '3613' ~ "",
    specimen_id == 'KRP446' ~ "",
    specimen_id == 'KRP460' ~ "",
    specimen_id == 'KRP461' ~ "",
    specimen_id == 'KRP462' ~ "",
    specimen_id == '30-10034' ~ ">=1000 mg/g cr",
    specimen_id == '32-10003' ~ "150 to <500 mg/g cr",
    specimen_id == '32-10034' ~ "<150 mg/g cr",
    specimen_id == '32-2' ~ "",
    specimen_id == '33-10005' ~ "",
    specimen_id == '33-10006' ~ "",
    specimen_id == '29-10006' ~ "500 to <1000 mg/g cr",
    specimen_id == '29-10008' ~ "",
    specimen_id == '29-10010' ~ "",
    specimen_id == '29-10012' ~ ">=1000 mg/g cr",
    specimen_id == '29-10013' ~ "<150 mg/g cr",
    specimen_id == '31-10000' ~ "",
    specimen_id == '31-10001' ~ "150 to <500 mg/g cr",
    specimen_id == '31-10006' ~ "",
    specimen_id == '31-10013' ~ "<150 mg/g cr",
    specimen_id == '31-10035' ~ ">=1000 mg/g cr"
  ))
#A1c
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(A1c = dplyr::case_when(
    specimen_id == '18-142' ~ "",
    specimen_id == '18-162' ~ "",
    specimen_id == '18-312' ~ "",
    specimen_id == '3490' ~ "",
    specimen_id == '3499' ~ "",
    specimen_id == '3504' ~ "",
    specimen_id == '3535' ~ "",
    specimen_id == '3593' ~ "",
    specimen_id == '3613' ~ "",
    specimen_id == 'KRP446' ~ "",
    specimen_id == 'KRP460' ~ "",
    specimen_id == 'KRP461' ~ "",
    specimen_id == 'KRP462' ~ "",
    specimen_id == '30-10034' ~ "6.5 to <7.5%",
    specimen_id == '32-10003' ~ "",
    specimen_id == '32-10034' ~ ">=8.5%",
    specimen_id == '32-2' ~ "<6.5%",
    specimen_id == '33-10005' ~ "",
    specimen_id == '33-10006' ~ "",
    specimen_id == '29-10006' ~ "",
    specimen_id == '29-10008' ~ "",
    specimen_id == '29-10010' ~ "",
    specimen_id == '29-10012' ~ "",
    specimen_id == '29-10013' ~ "6.5 to <7.5%",
    specimen_id == '31-10000' ~ "<6.5%",
    specimen_id == '31-10001' ~ "6.5 to <7.5%",
    specimen_id == '31-10006' ~ "<6.5%",
    specimen_id == '31-10013' ~ "<6.5%",
    specimen_id == '31-10035' ~ "6.5 to <7.5%"
  ))
#Albuminuria
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(albuminuria = dplyr::case_when(
    specimen_id == '18-142' ~ "",
    specimen_id == '18-162' ~ "",
    specimen_id == '18-312' ~ "",
    specimen_id == '3490' ~ "",
    specimen_id == '3499' ~ "",
    specimen_id == '3504' ~ "",
    specimen_id == '3535' ~ "",
    specimen_id == '3593' ~ "",
    specimen_id == '3613' ~ "",
    specimen_id == 'KRP446' ~ "",
    specimen_id == 'KRP460' ~ "",
    specimen_id == 'KRP461' ~ "",
    specimen_id == 'KRP462' ~ "",
    specimen_id == '30-10034' ~ ">=1000 mg/g cr",
    specimen_id == '32-10003' ~ "",
    specimen_id == '32-10034' ~ "",
    specimen_id == '32-2' ~ "",
    specimen_id == '33-10005' ~ "",
    specimen_id == '33-10006' ~ "",
    specimen_id == '29-10006' ~ "",
    specimen_id == '29-10008' ~ "30 to <300 mg/g cr",
    specimen_id == '29-10010' ~ "",
    specimen_id == '29-10012' ~ "",
    specimen_id == '29-10013' ~ "30 to <300 mg/g cr",
    specimen_id == '31-10000' ~ "",
    specimen_id == '31-10001' ~ "<30 mg/g cr",
    specimen_id == '31-10006' ~ "",
    specimen_id == '31-10013' ~ "",
    specimen_id == '31-10035' ~ ""
  ))
#Disease
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(disease_type = dplyr::case_when(
    specimen_id == '18-142' ~ "Healthy Reference",
    specimen_id == '18-162' ~ "Healthy Reference",
    specimen_id == '18-312' ~ "Healthy Reference",
    specimen_id == '3490' ~ "Healthy Reference",
    specimen_id == '3499' ~ "Healthy Reference",
    specimen_id == '3504' ~ "Healthy Reference",
    specimen_id == '3535' ~ "Healthy Reference",
    specimen_id == '3593' ~ "Healthy Reference",
    specimen_id == '3613' ~ "Healthy Reference",
    specimen_id == 'KRP446' ~ "Healthy Reference",
    specimen_id == 'KRP460' ~ "Healthy Reference",
    specimen_id == 'KRP461' ~ "Healthy Reference",
    specimen_id == 'KRP462' ~ "Healthy Reference",
    specimen_id == '30-10034' ~ "AKI",
    specimen_id == '32-10003' ~ "AKI",
    specimen_id == '32-10034' ~ "AKI",
    specimen_id == '32-2' ~ "AKI",
    specimen_id == '33-10005' ~ "AKI",
    specimen_id == '33-10006' ~ "AKI",
    specimen_id == '29-10006' ~ "CKD",
    specimen_id == '29-10008' ~ "CKD",
    specimen_id == '29-10010' ~ "CKD",
    specimen_id == '29-10012' ~ "CKD",
    specimen_id == '29-10013' ~ "CKD",
    specimen_id == '31-10000' ~ "CKD",
    specimen_id == '31-10001' ~ "CKD",
    specimen_id == '31-10006' ~ "CKD",
    specimen_id == '31-10013' ~ "CKD",
    specimen_id == '31-10035' ~ "CKD"
  ))
#Replicates
kpmp@meta.data <- kpmp@meta.data %>%
  mutate(Rep = dplyr::case_when(
    specimen_id == '18-142' ~ "HR_1",
    specimen_id == '18-162' ~ "HR_2",
    specimen_id == '18-312' ~ "HR_3",
    specimen_id == '3490' ~ "HR_4",
    specimen_id == '3499' ~ "HR_5",
    specimen_id == '3504' ~ "HR_6",
    specimen_id == '3535' ~ "HR_7",
    specimen_id == '3593' ~ "HR_8",
    specimen_id == '3613' ~ "HR_9",
    specimen_id == 'KRP446' ~ "HR_10",
    specimen_id == 'KRP460' ~ "HR_11",
    specimen_id == 'KRP461' ~ "HR_12",
    specimen_id == 'KRP462' ~ "HR_13",
    specimen_id == '30-10034' ~ "AKI_1",
    specimen_id == '32-10003' ~ "AKI_2",
    specimen_id == '32-10034' ~ "AKI_3",
    specimen_id == '32-2' ~ "AKI_4",
    specimen_id == '33-10005' ~ "AKI_5",
    specimen_id == '33-10006' ~ "AKI_6",
    specimen_id == '29-10006' ~ "CKD_1",
    specimen_id == '29-10008' ~ "CKD_2",
    specimen_id == '29-10010' ~ "CKD_3",
    specimen_id == '29-10012' ~ "CKD_4",
    specimen_id == '29-10013' ~ "CKD_5",
    specimen_id == '31-10000' ~ "CKD_6",
    specimen_id == '31-10001' ~ "CKD_7",
    specimen_id == '31-10006' ~ "CKD_8",
    specimen_id == '31-10013' ~ "CKD_9",
    specimen_id == '31-10035' ~ "CKD_10"
  ))

# Fig.4A. Subset DCT nuclei ----------------------------------------------------------------------------------------------------
Idents(kpmp) <- "subclass.l2"
kpmp.dct <- subset(kpmp, idents = c("DCT1", "DCT2", "dDCT", "cycDCT"))
kpmp.dct <- subset(kpmp.dct, subset = SLC12A3 > 0)
DimPlot(kpmp.dct, 
        reduction = "umap", 
        label.size = 6, 
        cols = c("#29A498", "orange", "#F48FB1", "red"), 
        pt.size = 0) + 
  xlim(6, 10) + 
  ylim(1, 5)

ggsave(
  "F4A.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4.5,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.4B-E. AddModuleScore ------------------------------------------------------------------------------------------------------
DCT1_marker_gene_list <- list(c("ERBB4", "EGF", "TRPM7", "FGF13", "COL5A2", "UMOD", "PTGFR", "STK32B", "RTL4", "ABCA13"))
DCT2_marker_gene_list <- list(c("SLC8A1", "ARL15", "CALB1", "SLC2A9", "PHACTR1", "GLS", "S100G", "KL", "KLK1", "EGFEM1"))
Ca_cassette <- list(c("SLC8A1", "VDR", "TRPV5", "CALB1", "S100G", "RYR2"))
Mg_cassette <- list(c("TRPM6", "TRPM7", "CNNM2", "EGF", "SLC41A3", "FXYD2", "UMOD", "PROX1"))

kpmp.dct <- AddModuleScore(object = kpmp.dct, features = DCT1_marker_gene_list, name = "DCT1_score")
kpmp.dct <- AddModuleScore(object = kpmp.dct, features = DCT2_marker_gene_list, name = "DCT2_score")
kpmp.dct <- AddModuleScore(object = kpmp.dct, features = Ca_cassette, name = "Ca_score")
kpmp.dct <- AddModuleScore(object = kpmp.dct, features = Mg_cassette, name = "Mg_score")

# Fig.4B
FeaturePlot(kpmp.dct,
            features = c("DCT1_score1")) + 
  xlim(6, 10) + 
  ylim(1, 5) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdYlBu")))

ggsave(
  "F4B.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.4C
FeaturePlot(kpmp.dct,
            features = c("DCT2_score1")) + 
  xlim(6, 10) + 
  ylim(1, 5) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdYlBu")))

ggsave(
  "F4C.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.4D
FeaturePlot(kpmp.dct,
            features = c("Mg_score1")) + 
  xlim(6, 10) + 
  ylim(1, 5) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdYlBu")))

ggsave(
  "F4D.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.4E
FeaturePlot(kpmp.dct,
            features = c("Ca_score1")) + 
  xlim(6, 10) + 
  ylim(1, 5) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 9, name = "RdYlBu")))

ggsave(
  "F4E.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.4F-G. DEGs -----------------------------------------------------------------------------------------------------------------

DCT.markers <- FindMarkers(kpmp.dct, ident.1 = "DCT1", ident.2 = "DCT2", min.pct = 0.1, logfc.threshold = 0, test.use = "DESeq2")
saveRDS(DCT.markers, here("DCT1vsDCT2_kpmp.dct.rda"))

cyc.markers <- FindMarkers(kpmp.dct, ident.1 = "cycDCT", min.pct = 0.1, logfc.threshold = 0, test.use = "DESeq2")
saveRDS(cyc.markers, here("cycDCT_kpmp.dct.rda"))

# Import the 'ALL DEG files.RData' file (supplemental data from original publication) into the R Studio environment

# Fig.4F -------------------------------------------------------------------------------------------------------------------------
x <- DCT.markers
x2 <- "Average log2FC KPMP"
y <- DCT1vsDCT2
y2 <- "Average log2FC Mouse"

x_tb <- x %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

y_tb <- y %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# convert mouse gene to human gene
y_tb <- mousegnameConverter(y_tb, "gene")

# Tables combined
z <- inner_join(x_tb, y_tb, by = "gene")

#Coorelation plot
rsqrd = lm(avg_log2FC.x~avg_log2FC.y, data=z)
r2 <- summary(rsqrd)$r.squared

r2

genes.to.label <- c("ARL15", "SLC8A1", "CALB1", "CASR",  "TRPM7", "EGF", "UMOD", "ERBB4", "SLC2A9", "TRPM6", "TRPV5", "KLK1")
df.label <- subset(z, gene %in% genes.to.label)

ggplot(z, aes(x = avg_log2FC.x, y = avg_log2FC.y, label=gene)) +
  theme_bw() +
  geom_point(color=dplyr::case_when((z$p_val_adj.x < 0.01 & z$p_val_adj.y < 0.01) ~ "#21918c",
                                    TRUE ~ "gray")) +
  geom_text_repel(data=df.label,
                  segment.size  = 0.2,
                  segment.color = "grey50") +
  geom_smooth (method=lm) +
  labs(x = x2, y = y2) +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")),   # R square
               label.x.npc = "left", label.y.npc = 0.90,
               rr.digits = 3) +
  xlim(-3, 3) +
  ylim(-4, 4)

ggsave(
  "F4F.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)

# Fig.4G -------------------------------------------------------------------------------------------------------------------------
x <- cyc.markers
x2 <- "Average log2FC kpmp"
y <- Prolif
y2 <- "Average log2FC Mouse"

x_tb <- x %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

y_tb <- y %>%
  data.frame() %>%
  rownames_to_column(var="gene") %>% 
  as_tibble()

# convert mouse gene to human gene
y_tb <- mousegnameConverter(y_tb, "gene")

# Tables combined
z <- inner_join(x_tb, y_tb, by = "gene")

#Correlation plot
rsqrd = lm(avg_log2FC.x~avg_log2FC.y, data=z)
r2 <- summary(rsqrd)$r.squared

r2

#Set Range for Far Right Data Points
df.upper <- subset(z, avg_log2FC.x > 0.4 & avg_log2FC.y > 0.4)
#Set Range for Far Left Data Points
df.lower <- subset(z, avg_log2FC.x < -0.4 & avg_log2FC.y < -0.4)

genes.to.label <- c("TOP2A", "CENPP", "DIAPH3", "CENPF",  "EZH2", "NAV2", "FXYD2", "ATP1A1", "UMOD", "TRPM6", "SLC12A3", "EGF")
df.label <- subset(z, gene %in% genes.to.label)

ggplot(z, aes(x = avg_log2FC.x, y = avg_log2FC.y, label=gene)) +
  theme_bw() +
  geom_point(color=dplyr::case_when((z$p_val_adj.x < 0.01 & z$p_val_adj.y < 0.01) ~ "#21918c",
                                    TRUE ~ "gray")) +
  geom_text_repel(data=df.label,
                  segment.size  = 0.2,
                  segment.color = "grey50") +
  geom_smooth (method=lm) +
  labs(x = x2, y = y2) +
  stat_poly_eq(aes(label = paste(..rr.label.., sep = "~~~")),   # R square
               label.x.npc = "left", label.y.npc = 0.90,
               rr.digits = 3) +
  xlim(-2.5, 2.5) +
  ylim(-2.5, 2.5)

ggsave(
  "F4G.tiff",
  device = "tiff",
  plot = last_plot(),
  width = 4,
  height = 3,
  units = c("in"),
  dpi = 700,
  compression = 'lzw'
)
