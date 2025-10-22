library("tidyverse")
library("readxl")
library("ggbeeswarm")
library("ggplotify")
library("limma")
library("tidymodels")
library("ggrepel")
library("vip")
library("pROC")
library("patchwork")
library("pheatmap")
library("ggupset")
library("plotly")
library("ggiraph")

getwd()
data_soma <- read_tsv("../../Data/Somascan/data_phase2_somalogic_curated_20250922.tsv")
meta_soma <-  read_tsv("../../Data/Somascan/meta_phase2_somalogic_20250922.tsv")
data_ht_all <- read_tsv("../../Data/Olink/phase2_unbridged_data_overlapping_HT_pandisease_20250926.tsv")
manifest <- read_tsv("../../Data/Olink/metadata_olink_soma_phase2.tsv")
#manifest <- read_excel("../metadata/samples_2025-05-12.xlsx")
#data_ht_all <- read_tsv( "../olink_data/phase2_unbridged_data.tsv")

# dengue_id <- 
#   manifest |> 
#   filter(Disease == "Dengue")

# data_ht_all |> 
#   filter(DAid %in% dengue_id$DAid)


# a <- 
# data_ht_all |> 
#   distinct(DAid) |> 
#   left_join(manifest) |> 
#   distinct(Disease)


# data_soma |> 
#   distinct(DAid) |> 
#   left_join(manifest) |> 
#   distinct(Disease) |> 
#   filter(!Disease %in% a$Disease)


data_soma_all <- 
  data_soma |> 
    pivot_longer(names_to = "Assay", values_to = "RFU", cols = -DAid)

manifest_common <- 
  manifest |> 
  filter(!(Disease == "Lung cancer" & Cohort == "ANMU"),
Disease != "Longevity") 


# |> 

#     filter(DAid %in% data_soma_all$DAid) |> 
#     mutate(Disease = case_when(Disease == "Pancreatic ductal adenocarcinoma" ~ "Pancreatic cancer", 
#   Diagnose == "OVC" ~ "Ovarian cancer", 
# Diagnose == "BRC" ~ "Breast cancer",
# Diagnose == "LUNG" ~ "Lung cancer",
# Diagnose == "CRC" ~ "Colorectal cancer",
# Diagnose == "PRC" ~ "Prostate cancer",
# T ~ Disease),
# Disease = ifelse(Cohort == "EPIL", "Epilepsy", Disease),
# Disease = ifelse(Cohort == "PREG", "Pregnancy", Disease),
# Disease = ifelse(Cohort %in% c("FIBR", "PARD"), Diagnose, Disease),
# Disease = ifelse(Disease %in% c("Healthy control", "healthy", "Control"), "Healthy", Disease),
# Class = case_when(Disease == "PD" ~ "Neurologic",
# Disease %in% c("Breast cancer", "Lung cancer", "Prostate cancer", "Colorectal cancer", "Ovarian cancer") ~ "Cancer",
# Disease == "Pregnancy" ~ "Healthy",
# Disease == "Healthy" ~ "Healthy",
# T ~ Class))  |> 
#   filter(!is.na(Disease),
# Disease != "Longevity")



# manifest_olink <- 
#   manifest |> 
#     filter(!(Disease == "Lung cancer" & Cohort == "ANMU")) |> 
#     filter(DAid %in% data_ht_all$DAid) |> 
#     mutate(Disease = case_when(Disease == "Pancreatic ductal adenocarcinoma" ~ "Pancreatic cancer", 
#   Diagnose == "OVC" ~ "Ovarian cancer", 
# Diagnose == "BRC" ~ "Breast cancer",
# Diagnose == "LUNG" ~ "Lung cancer",
# Diagnose == "CRC" ~ "Colorectal cancer",
# Diagnose == "PRC" ~ "Prostate cancer",
# T ~ Disease),
# Disease = ifelse(Cohort == "EPIL", "Epilepsy", Disease),
# Disease = ifelse(Cohort %in% c("FIBR", "PARD"), Diagnose, Disease))  |> 
#   filter(!is.na(Disease))

