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

setwd("/Users/mariabueno/Library/CloudStorage/OneDrive-KTH/Repos/olink-somalogic-exploration")
data_soma <- read_tsv("../soma_data/data_phase2_somalogic_curated_20250922.tsv")
meta_soma <-  read_tsv("../soma_data/meta_phase2_somalogic_20250922.tsv")
#data_ht_all <- read_tsv( "../olink_data/phase2_unbridged_data.tsv")
data_ht_all <- read_tsv("../olink_data/phase2_unbridged_data_overlapping_HT_pandisease_20250926.tsv")
manifest <- read_excel("../metadata/samples_2025-05-12.xlsx")

data_soma_all <- 
  data_soma |> 
    pivot_longer(names_to = "Assay", values_to = "RFU", cols = -DAid)

manifest_common <- 
  manifest |> 
        filter(!(Disease == "Lung cancer" & Cohort == "ANMU")) |> 

    filter(DAid %in% data_soma_all$DAid) |> 
    mutate(Disease = case_when(Disease == "Pancreatic ductal adenocarcinoma" ~ "Pancreatic cancer", 
  Diagnose == "OVC" ~ "Ovarian cancer", 
Diagnose == "BRC" ~ "Breast cancer",
Diagnose == "LUNG" ~ "Lung cancer",
Diagnose == "CRC" ~ "Colorectal cancer",
Diagnose == "PRC" ~ "Prostate cancer",
T ~ Disease),
Disease = ifelse(Cohort == "EPIL", "Epilepsy", Disease),
Disease = ifelse(Cohort == "PREG", "Pregnancy", Disease),
Disease = ifelse(Cohort %in% c("FIBR", "PARD"), Diagnose, Disease),
Disease = ifelse(Disease %in% c("Healthy control", "healthy", "Control"), "Healthy", Disease),
Class = case_when(Disease == "PD" ~ "Neurologic",
Disease %in% c("Breast cancer", "Lung cancer", "Prostate cancer", "Colorectal cancer", "Ovarian cancer") ~ "Cancer",
Disease == "Pregnancy" ~ "Healthy",
Disease == "Healthy" ~ "Healthy",
T ~ Class))  |> 
  filter(!is.na(Disease),
Disease != "Longevity")



manifest_olink <- 
  manifest |> 
    filter(!(Disease == "Lung cancer" & Cohort == "ANMU")) |> 
    filter(DAid %in% data_ht_all$DAid) |> 
    mutate(Disease = case_when(Disease == "Pancreatic ductal adenocarcinoma" ~ "Pancreatic cancer", 
  Diagnose == "OVC" ~ "Ovarian cancer", 
Diagnose == "BRC" ~ "Breast cancer",
Diagnose == "LUNG" ~ "Lung cancer",
Diagnose == "CRC" ~ "Colorectal cancer",
Diagnose == "PRC" ~ "Prostate cancer",
T ~ Disease),
Disease = ifelse(Cohort == "EPIL", "Epilepsy", Disease),
Disease = ifelse(Cohort %in% c("FIBR", "PARD"), Diagnose, Disease))  |> 
  filter(!is.na(Disease))



# pal_type <- 
#   c("Disease" = "darkred", 
# "Related diseases" = "", 
# "Other diseases" = "", 
# "Healthy" = "grey")
pal_type <- c(
  "Disease"          = "#B22222",   # firebrick/darkred (strong, stands out)
  "Related diseases" = "#E06F1F",   # warm orange
  "Other diseases"   = "#E5C494",   # soft beige / muted tan
  "Healthy"          = "grey60"     # neutral grey
)

pal_de <-
  c(
    "not significant" = "#D3D3D3",
    "significant up" = "#FF7176",
    "significant down" = "#92C9DA"
  )

theme_hpa <-
  function(angled = F,
           axis_x = T,
           axis_y = T,
           facet_title = T) {
    t <-
      theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.2, "lines"),
        panel.background = element_rect(fill = "white"),
        panel.border = element_blank(),
        plot.title = element_text(
          face = "bold",
          size = rel(1),
          hjust = 0.5
        ),
        plot.subtitle = element_text(
          face = "bold",
          hjust = 0.5,
          size = rel(1),
          vjust = 1
        ),
        axis.title = element_text(face = "bold", size = rel(1)),
        axis.ticks.length = unit(.25, "cm"),
        axis.line = element_line(linewidth = 0.5),
        axis.text = element_text(size = rel(1), color = 'black'),
        legend.key = element_blank(),
        legend.position = "right",
        legend.text = element_text(size = rel(0.8)),
        legend.key.size = unit(0.7, "cm"),
        legend.title = element_text(size = rel(1)),
        plot.margin = unit(c(10, 5, 5, 5), "mm"),
        strip.background = element_rect(colour = "grey90", fill = "grey90"),
        strip.text = element_text(face = "bold")
      )
    
    if (angled) {
      t <-
        t + theme(axis.text.x = element_text(
          angle = 90,
          vjust = 0.5,
          hjust = 1
        ))
    }
    
    if (axis_x == F) {
      t <- t +
        theme(
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line.x = element_blank(),
          axis.title.x = element_blank()
        )
    }
    
    if (axis_y == F) {
      t <- t +
        theme(
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.line.y = element_blank(),
          axis.title.y = element_blank()
        )
    }
    if (facet_title == F) {
      t <- t + theme(strip.text = element_blank())
    }
    return(t)
  }

savepath <-
  function(savename) {
    result_folder <- paste0("results/", Sys.Date())
    dir.create(result_folder, showWarnings = FALSE)

    savename <-
      paste0(result_folder, "/", savename)


    return(savename)

  }
 
do_limma_disease <-
  function(data_wide,
           metadata,
           disease,
           controls,
           correct = T,
           cutoff = 0) {
    # Select current disease
    dat <-
      data_wide %>%
      inner_join(metadata %>%
                   select(DAid, Sex, Age, BMI, Disease), by = "DAid") %>%
      rename(Group = Disease) %>%
      mutate(Group = ifelse(Group == disease, "1_Case", "0_Control"))
    
    n_males <-
      metadata |>
      filter(Disease == disease,
             Sex == "M") |>
      nrow()
    
    n_females <-
      metadata |>
      filter(Disease == disease,
             Sex == "F") |>
      nrow()
    
    # Filter samples for sex-specific cancers
    if (n_males == 0) {
      
      metadata <- 
        metadata |> 
        filter(Sex == "F")
      
      dat <- 
        dat |> 
        filter(DAid %in% metadata$DAid)
    } else if (n_females == 0) {
      
      metadata <- 
        metadata |> 
        filter(Sex == "M")
      
      dat <- 
        dat |> 
        filter(DAid %in% metadata$DAid)
    } else {
      
      dat <- dat
      metadata <- metadata
      
    }
    
    # Filter missing age and sex
    if (correct == T) {
      dat <-
        dat |>
        filter(!is.na(Sex), !is.na(Age))
    } else {
      dat <- dat
    }
    
    # Design a model
    if (correct == T) {
      if (n_males == 0 | n_females == 0) {
        design <- model.matrix(~ 0 + as.factor(dat$Group) + dat$Age)
        colnames(design) <- c("control", "case", "Age")
      } else {
        design <-
          model.matrix(~ 0 + as.factor(dat$Group) + as.factor(dat$Sex) + dat$Age)
        colnames(design) <- c("control", "case",  "Sex", "Age")
      }
      
    } else {
      design <- model.matrix(~ 0 + as.factor(dat$Group))
      colnames(design) <- c("control", "case")
    }
    
    # Make contrast
    contrast <-
      makeContrasts(Diff = case - control, levels = design)
    
    # Fit linear model to each protein assay
    dat_fit <-
      dat %>%
      select(-Sex, -Age, -BMI, -Group)  %>%
      column_to_rownames("DAid") %>%
      t()
    
    fit <-
      lmFit(dat_fit, design = design, maxit = 100000) #method = "robust",
    
    # Apply contrast
    contrast_fit <- contrasts.fit(fit, contrast)
    
    # Apply empirical Bayes smoothing to the SE
    ebays_fit <- eBayes(contrast_fit, robust = T)
    
    # Extract DE results
    DE_results <-
      topTable(
        ebays_fit,
        n = nrow(ebays_fit$p.value),
        adjust.method = "fdr",
        confint = TRUE
      )
    
    DE_res <-
      DE_results %>%
      as_tibble(rownames = "Assay") %>%
      mutate(
        Disease = disease,
        sig = case_when(
          adj.P.Val < 0.05 & logFC < -cutoff ~ "significant down",
          adj.P.Val < 0.05 &
            logFC > cutoff ~ "significant up",
          T ~ "not significant"
        ),
        Control = controls
      )
    
    return(DE_res)
  }

  ## Generate volcano plot from differential expression results                                                                                                                                  
plot_volcano <- function(de_results, cutoff = 0, labels_balanced = F, multiple = F, nrow = 1) {
  
  if (labels_balanced) {
    labels <- 
      de_results |> 
      filter(sig == "significant up") |> 
      top_n(n = 10, wt = -log10(adj.P.Val))  |> 
      bind_rows(de_results |> 
                  filter(sig == "significant down") |> 
                  top_n(n = 10, wt = -log10(adj.P.Val)) 
      )
    
  } else {
    labels <- 
      de_results |> 
      top_n(n = 10, wt = -log10(adj.P.Val)) 
  }
  
  
  volcano_plot <- 
    de_results |> 
    ggplot(aes(x = logFC, y = -log10(adj.P.Val), color = sig, label = Assay)) +
    geom_point(size = 1, alpha = 0.4, show.legend = F) + 
    geom_text_repel(data = labels, aes(label = Assay), size = 2, show.legend = F) +
    geom_hline(yintercept = -log10(0.05), linetype = 'dashed', color = "darkgrey") +
    geom_vline(xintercept = -cutoff, linetype = 'dashed', color = "darkgrey") +
    geom_vline(xintercept = cutoff, linetype = 'dashed', color = "darkgrey") +
    scale_color_manual(values = pal_de) +
    theme_hpa() +   
    theme(axis.text = element_text(size = 8),
          legend.position = "top") 
  
  if(multiple == T) {
    
    volcano_plot <- 
      volcano_plot +
      facet_wrap(~Disease, scales = "free", nrow = nrow) 
  }
  
  return(volcano_plot)
}


plot_boxplot <- function(proteins,
                         data,
                         metadata,
                         platform = "Olink_HT",
                         title = "") {
  
  if (platform == "Olink_HT") {
    data_filtered <-
      data |>
      filter(Assay %in% proteins) |>
      select(DAid, Assay, NPX) |>
      left_join(metadata |>
                  select(DAid, Disease),
                by = "DAid") 
    

    boxplot <-
        data_filtered |>
        ggplot(aes(Disease,
                   NPX,
                   fill = Disease,
                   color = Disease)) +
        geom_quasirandom(alpha = 0.7) +
        geom_boxplot(
          alpha = 0.3,
          outlier.color = NA,
          color = "grey20"
        ) +
        scale_color_manual(values = pal_cancers) +
        scale_fill_manual(values = pal_cancers) +
        facet_wrap( ~ Assay, scales = "free_y") +
        theme_hpa(angled = T) +
        xlab("") +
        ggtitle(title)

    
  } else if (platform == "1.5K") {
    boxplot <- ggplot()
  } else {
    stop("Platform not recognized")
  }
  
  return(boxplot)
}





translate_soma <- function(soma_df) {
      soma_df |> 
  left_join(meta_soma |> 
              distinct(Assay = AptName, EntrezGeneSymbol), by = "Assay") |> 
  rename(AptName = Assay,
  Assay = EntrezGeneSymbol) |>
  relocate(Assay)

} 

# Cohort palette
# Class palette
library(RColorBrewer)
getPalette3 = colorRampPalette(brewer.pal(8, "Set2"))
pal_class <- getPalette3(8)
names(pal_class) <-
  c(
    "Neurologic",
    "Cardiovascular",
    "Cancer",
    "Autoimmune",
    "Pregnancy",
    "Infection",
    "Metabolic",
    "Healthy"
  )

disease_n <- 
    manifest_common |> 
    count(Class, Disease) |> 
    mutate(Class = factor(Class)) |> 
    arrange(Class) 

disease_levels <- unique(disease_n$Disease)
