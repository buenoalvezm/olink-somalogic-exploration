###############################################################################
# Plasma proteomics disease classification pipeline
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

source("scripts/functions_ml.R")

# ===============================
# 2. Global settings
# ===============================
seed <- 8
set.seed(seed)

run_mode <- c("olink", "soma", "combined")
# Choose ONE:
run_mode <- "combined"

# Parallelization
n_workers <- parallel::detectCores() - 1
plan(multicore, workers = n_workers)

# ===============================
# 3. Output directories
# ===============================
dir.create("results_olink/model", recursive = TRUE, showWarnings = FALSE)
dir.create("results_soma/model", recursive = TRUE, showWarnings = FALSE)
dir.create("results_combined/model", recursive = TRUE, showWarnings = FALSE)

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
  filter(DAid %in% data_soma_wide$DAid) |>
  count(Class, Disease) |>
  filter(n < 45 | Class == "Healthy")

resource_meta <- 
  meta |>
  filter(!Disease %in% exclude_diseases$Disease)

# ===============================
# 6. Prepare Olink data (wide)
# ===============================
olink_ml <- 
  data_olink |>
  filter(DAid %in% resource_meta$DAid) |>
  select(DAid, Assay, NPX) |>
  left_join(
    resource_meta |> select(DAid, Disease),
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

ml_daid <- 
  data_combined |> 
  distinct(DAid)

# ===============================
# 8. Train/test split (ONCE)
# ===============================
current_split <- generate_split(
  data = data_combined,
  proportion = 0.7,
  variable_stratify = "Disease",
  seed = seed
)

saveRDS(current_split, paste0("data/processed/current_split_", seed, ".rds"))


test_set_daid <- 
  current_split$data_test |> 
  distinct(DAid)

# Subset data according to run mode
train_data <- subset_assays(current_split$data_train , run_mode)
test_data  <- subset_assays(current_split$data_test ,  run_mode)

#print_split_summary(train_data, test_data, run_mode)


# Run disease-vs-all models
results_path <- paste0("results_", run_mode, "/")

multiclass_results <- disease_against_all(
  split_train = train_data,
  split_test  = test_data,
  path_save   = results_path,
  seed        = seed
)


# Save results
saveRDS(
  multiclass_results,
  paste0("data/processed/ml_", seed, "_", run_mode, ".rds")
)

