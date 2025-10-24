library("tidyverse")
library("readxl")

# BAMSE analyses in Soma
data_soma <- read_tsv("../../Data/Somascan/data_phase2_somalogic_curated_20250922.tsv")
manifest <- read_excel("../../Data/metadata/samples_2025-05-12.xlsx")


manifest |> 
  filter(Cohort == "BAMS") |> 
  select(DAid, Age, Sex, BMI, `Extra data`)


# UMAP
umap_bamse <-
  do_umap(data = data_bamse_detected,
          wide = F,
          plots = F)

plot_umap <-
  umap_bamse |>
  left_join(meta_bamse, by = c("Sample" = "DAid")) |>
  mutate(Visit = factor(Visit)) |>
  ggplot(aes(UMAP1, UMAP2, color = Visit, fill = Visit)) +
  geom_point(alpha = 0.7) +
  stat_ellipse(geom = "polygon",
               alpha = 0.3,
               color = NA) +
  scale_color_manual(values = pal_bamse) +
  scale_fill_manual(values = pal_bamse) +
  theme_hpa() +
  ggtitle("UMAP")

ggsave(
  savepath_results("Manuscript-figures", "BAMSE_UMAP_v2.pdf"),
  h = 5,
  w = 6
)



# MM
data_mm <-
  data_bamse_detected |>
  select(DAid, Assay, NPX) %>%
  left_join(meta_bamse, by = "DAid") |>
  mutate(Sex = as.factor(Sex),
         Age = as.factor(Age),
         Subject = as.factor(Subject))

# Run mixed effect models on all proteins extracting the fixed effects
fixed_effects <-
  map_df(unique(data_mm$Assay), function(protein) {
    do_mixed_effect_model(df = data_mm,
                          type = "fixed_effects",
                          protein = protein)
    
  })

# Run mixed effect models on all proteins extracting the variance explained
variance_explained <-
  map_df(unique(data_mm$Assay), function(protein) {
    do_mixed_effect_model(df = data_mm,
                          type = "variance_explained",
                          protein = protein)
  })


# Variance explained
variance_explained |>
  filter(Component != "Residual") |>
  mutate(Variance = Variance * 100,
         Component = factor(
           Component,
           levels = c("Random effects (subject)", "Fixed effects (age & sex)")
         )) |>
  ggplot(aes(Variance, fill = Component)) +
  geom_histogram(show.legend = F) +
  facet_wrap( ~ Component, scales = "free") +
  scale_fill_manual(values = pal_mm) +
  xlab("Variance explained (%)") +
  ylab("Number of proteins") +
  theme_hpa()

ggsave(
  savepath_results("Manuscript-figures", "MM_variance_explained_BAMSE_v2.pdf"),
  h = 3,
  w = 7
)

# Supplementary
order <-
  variance_explained |>
  filter(Component == "Residual") |>
  arrange(Variance)

variance_explained |>
  filter(Component != "Residual") |>
  mutate(Protein = factor(Protein, levels = order$Protein)) |>
  mutate(Variance = Variance * 100) |>
  ggplot(aes(Protein, Variance, fill = Component)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = pal_mm) +
  theme_hpa() +
  theme(
    legend.position = "top",
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank()
  ) +
  ylab("Variance explained (%)")

ggsave(
  savepath_results("Manuscript-figures", "Fig2_MM-variance-explained_v2.pdf"),
  h = 5,
  w = 14
)