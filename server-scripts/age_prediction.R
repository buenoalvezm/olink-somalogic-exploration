###############################################################################
# Plasma proteomics age prediction pipeline
###############################################################################

# ===============================
# 1. Libraries
# ===============================
library(tidyverse)
library(tidymodels)
library(themis)
library(multiROC)
library(furrr)
library(future.callr)

source("server-scripts/functions_ml.R")

# ===============================
# 2. Global settings
# ===============================
seed <- 3
set.seed(seed)

run_mode <- c("olink", "soma", "combined")
# Choose ONE:
run_mode <- "combined"

dir.create(paste0("results_", run_mode ,"_age/model"), recursive = TRUE, showWarnings = FALSE)


# Parallelization
n_workers <- parallel::detectCores() - 1
plan(multicore, workers = n_workers)

# ===============================
# 3. Output directories
# ===============================
#dir.create("results_olink_age/model", recursive = TRUE, showWarnings = FALSE)
#dir.create("results_soma_age/model", recursive = TRUE, showWarnings = FALSE)
#dir.create("results_combined_age/model", recursive = TRUE, showWarnings = FALSE)

# ===============================
# 4. Data loading
# ===============================
data_soma_wide <- read_tsv(
  "../../Data/Somascan/processed/data_phase2_somalogic_curated_subset_20251031.tsv"
)

data_olink <- read_tsv(
  "../../Data/Olink/processed/data_phase2_HT_bridged_curated_subset_20251031.tsv"
)

meta <- read_tsv(
  "../../Data/metadata/meta_phase2_hpav25_20251031.tsv"
)

# ===============================
# 5. Metadata filtering
#    - Remove rare diseases (<45 samples)
#    - Remove Healthy class
# ===============================
exclude_diseases <- 
  meta |>
  filter(DAid %in% data_soma_wide$DAid,
        Disease %in% c("BAMSE - visit 16", "BAMSE - visit 24", "BAMSE - visit 4", "Healthy (Children)", "Pregnancy")) |>
  distinct(Disease)


resource_meta <- 
  meta |>
  filter(!Disease %in% exclude_diseases$Disease)

ages <- readRDS("data/AAV DAid Age.rds")

# meta |> 
#   filter(DAid %in% ages$DAid) |> 
#   left_join(ages) |> 
#   ggplot(aes(Age, Age.at.sampling)) +
#   geom_point() +
#   theme_hpa()

# ggsave("results/age_vasculitis.png", h = 6, w = 6)
# DA16721 (79 originally, 82 really)

resource_meta <- 
  resource_meta |> 
  left_join(ages) |> 
  mutate(Age = ifelse(!is.na(Age.at.sampling), Age.at.sampling, Age))


# ===============================
# 6. Prepare Olink data (wide)
# ===============================
olink_ml <- 
  data_olink |>
  filter(DAid %in% resource_meta$DAid) |>
  select(DAid, Assay, NPX) |>
  left_join(
    resource_meta |> select(DAid, Age),
    by = "DAid"
  ) |>
  pivot_wider(
    names_from = Assay,
    values_from = NPX
  )

# ===============================
# 7. Combine Soma + Olink
#    (only samples present in both)
# ===============================
data_combined <- 
  data_soma_wide |>
  filter(DAid %in% olink_ml$DAid) |>
  left_join(olink_ml, by = c("DAid"))

current_split <- generate_split(
  data = data_combined,
  proportion = 0.7,
  variable_stratify = "Age",
  seed = seed
)

saveRDS(current_split, paste("data/processed/current_split_age_", seed, ".rds"))


train_data <- subset_assays(current_split$data_train, run_mode)
test_data  <- subset_assays(current_split$data_test,  run_mode)

#print_split_summary(train_data, test_data, run_mode)


results_path <- paste0("results_", run_mode, "_age/")

age_ml_results_splits <-
  
  continuous_prediction(
    split_train = train_data,
    split_test = test_data,
    variable_predict = "Age",
    path_save =  results_path,
    seed = seed
  )

saveRDS(
  age_ml_results_splits,
  paste0("data/processed/ml_age_", run_mode, ".rds")
)





############
run_mode <- "olink"

dir.create("results_olink_age/model", recursive = TRUE, showWarnings = FALSE)


train_data <- subset_assays(current_split$data_train, run_mode)
test_data  <- subset_assays(current_split$data_test,  run_mode)

print_split_summary(train_data, test_data, run_mode)


results_path <- paste0("results_", run_mode, "_age/")

age_ml_results_splits <-
  
  continuous_prediction(
    split_train = train_data,
    split_test = test_data,
    variable_predict = "Age",
    path_save =  results_path,
    seed = seed
  )

saveRDS(
  age_ml_results_splits,
  paste0("data/processed/ml_age_", run_mode, ".rds")
)
run_mode <- "combined"

dir.create("results_combined_age/model", recursive = TRUE, showWarnings = FALSE)

train_data <- subset_assays(current_split$data_train, run_mode)
test_data  <- subset_assays(current_split$data_test,  run_mode)

print_split_summary(train_data, test_data, run_mode)


results_path <- paste0("results_", run_mode, "_age/")

age_ml_results_splits <-
  
  continuous_prediction(
    split_train = train_data,
    split_test = test_data,
    variable_predict = "Age",
    path_save =  results_path,
    seed = seed
  )

saveRDS(
  age_ml_results_splits,
  paste0("data/processed/ml_age_", run_mode, ".rds")
)

