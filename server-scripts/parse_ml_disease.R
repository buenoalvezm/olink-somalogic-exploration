# Parse results
library(tidyverse)
library(ggbeeswarm)

source("scripts/functions/themes_levels_palettes.R")
meta <- read_tsv("../../Data/metadata/meta_phase2_hpav25_20251031.tsv")
classes <- 
  meta |> 
  distinct(Disease, Class) |> 
  mutate(Class = factor(Class, levels = disease_class_order)) |> 
  arrange(Class)

combined_df <- map_dfr(1:10, function(s) {
  file_path <- paste0("results_combined/important_proteins_multiclass_seed_", s, ".tsv")
  
  if (file.exists(file_path)) {
    read_tsv(file_path, show_col_types = FALSE) %>%
      mutate(seed = s)
  } else {
    NULL
  }
})

combined_df |> write_tsv("data/processed/ml_disease_seeds_importance.tsv")

combined_df |> 
  group_by(seed) |> 
  count(class) |> 
  left_join(classes, by = c("class" = "Disease")) |> 
  mutate(class = factor(class, levels = classes$Disease)) |> 
  ggplot(aes(class, n, color = Class, fill = Class)) +
  geom_quasirandom() +
  geom_boxplot(outlier.colour = NA, alpha = 0.5) +
  scale_color_manual(values = pal_class) +
  scale_fill_manual(values = pal_class) +
  theme_hpa(angled= T) +
  xlab("") +
  ylab("Number of features") +
  ggtitle("Stability across seeds")

ggsave(savepath("disease_n_features_seeds.png"), h = 5, w = 8)



plot_top_proteins_seed <- function(protein_importance,
                                  disease,
                                   n = 25) {
  top_proteins <-
    protein_importance |>
    filter(class == disease) |> 
    filter(abs(estimate) > 0) |>
    group_by(term) |>
    summarise(avg_importance = mean(abs(estimate)),
              n = n_distinct(seed)) |>
    arrange(-avg_importance) |>
    head(n)
  
  protein_importance |>
    filter(class == disease) |> 
    mutate(sign = ifelse(estimate > 0, "Positive", "Negative")) |> 
    filter(term %in% top_proteins$term) |>
    mutate(
      term = factor(term, levels = rev(top_proteins$term))
    ) |>
    ggplot(aes(fct_reorder(term, abs(estimate)), abs(estimate), fill = sign)) +
    geom_quasirandom(size = 0.5, aes(color = sign)) +
    geom_boxplot(alpha = 0.5, outlier.color = NA) +
    scale_color_manual(values = c(
      "Positive" = "#FF7176",
      "Negative" = "#92C9DA"
    )) +
    scale_fill_manual(values = c(
      "Positive" = "#FF7176",
      "Negative" = "#92C9DA"
    )) +
    coord_flip() +
    xlab("")  +
    ylab("Protein importance") +
    theme_hpa() +
    ggtitle(disease)
}


plot_top_proteins_seed(protein_importance = combined_df, disease = "Ovarian cancer")

pdf(savepath("top_25_disease_seeds.pdf"), h = 6, w = 6)
map(unique(combined_df$class), ~plot_top_proteins_seed(protein_importance = combined_df, disease = .x))
dev.off()


combined_df |> 
  group_by(class) |> 
  count(term) |> 
  ggplot(aes(n)) +
  geom_histogram() +
  facet_wrap(~class, scales = "free_y") +
  theme_hpa()

n_all_seeds <- 
  combined_df |> 
  group_by(class) |> 
  count(term) |> 
  arrange(-n) |> 
  filter(n == 9) |> 
  count(class) |> 
  left_join(classes, by = c("class" = "Disease")) 

n_all_seeds |> 
  mutate(class = factor(class, levels = classes$Disease)) |> 
  ggplot(aes(class, n, fill = Class)) +
  geom_col() +
  scale_fill_manual(values = pal_class) +
  theme_hpa(angled = T) +
  ggtitle("Number of features selected in all seeds")

ggsave(savepath("n_features_all_seeds.png"), h = 6, w = 7)

n_more_half_seeds <- 
  combined_df |> 
  group_by(class) |> 
  count(term) |> 
  arrange(-n) |> 
  filter(n > 5) |> 
  count(class) |> 
  left_join(classes, by = c("class" = "Disease")) |> 
  mutate(class = factor(class, levels = classes$Disease)) 

n_more_half_seeds |> 
  ggplot(aes(class, n, fill = Class)) +
  geom_col() +
  scale_fill_manual(values = pal_class) +
  theme_hpa(angled = T) +
  ggtitle("Number of features selected in > 5 seeds")

ggsave(savepath("n_features_half_seeds.png"), h = 6, w = 7)


combined_df |> 
  group_by(class) |> 
  count(term) |> 
  arrange(-n) |> #filter(term == "ASCL2")
  filter(n > 5) |> 
  ungroup() |> 
  count(term) |> 
  arrange(-n) |> 
  filter(n > 1) |> 
  ggplot(aes(term, n)) +
  geom_col() +
  theme_hpa(angled = T) +
  ggtitle("Proteins selected for more than one disease")

ggsave(savepath("features_half_seeds.png"), h = 6, w = 7)
