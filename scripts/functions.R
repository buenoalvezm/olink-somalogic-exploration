library(effectsize)

savepath <-
  function(savename) {
    result_folder <- paste0("results/", Sys.Date())
    dir.create(result_folder, showWarnings = FALSE)

    savename <-
      paste0(result_folder, "/", savename)


    return(savename)

  }

# Generate UMAP
do_umap <- function(data,
                    meta = NULL,
                    variable = NULL,
                    wide = T,
                    impute = T,
                    plots = F,
                    seed = 213,
                    n_neighbors = 15) {
  if (wide) {
    data_w <-
      data |>
      rename(Sample = DAid)
  } else {
    data_w <-
      data |>
      select(Sample = DAid, Assay, NPX) |>
      pivot_wider(values_from = NPX,
                  names_from = Assay)
  }
  
  if (impute) {
    umap_rec <-
      recipe( ~ ., data = data_w) %>%
      update_role(Sample, new_role = "id")  |>
      step_zv(all_predictors()) |> 
      step_normalize(all_predictors()) |>
      step_impute_knn(all_predictors()) |>
      step_umap(all_predictors(), neighbors = n_neighbors, seed = c(seed, seed))
    
    set.seed(seed)
    umap_prep <- prep(umap_rec)
    
  } else {
    umap_rec <-
      recipe( ~ ., data = data_w) %>%
      update_role(Sample, new_role = "id")  |>
      step_normalize(all_predictors()) |>
      step_umap(all_predictors(), neighbors = n_neighbors, seed = c(seed, seed))
    
    set.seed(seed)
    umap_prep <- prep(umap_rec)
    
  }
  
  umap_res <-  juice(umap_prep)
  
  if (plots) {
    # Loadings plot
    umap_plot <-
      umap_res |>
      left_join(meta |>
                  rename(Sample = DAid), by = "Sample") |>
      ggplot(aes(UMAP1, UMAP2, color = !!sym(variable))) +
      geom_point(alpha = 0.7, size = 2) +
      theme_hpa()
    
    return(list("umap_res" = umap_res,
                "umap_plot" = umap_plot))
  } else {
    return(umap_res)
  }
  
}



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


continuous_prediction <-
  function(split_train,
           split_test,
           variable_predict,
           seed) {
    cat(paste0("\nPreparing training and testing data for ", variable_predict))
    
    split_train <-
      split_train |>
      rename(Variable = !!sym(variable_predict))
    
    split_test <-
      split_test |>
      rename(Variable = !!sym(variable_predict))
    
    variable_split <- make_splits(split_train,
                                  split_test)
    
    cat(paste0("\nDefining ML specs for ", variable_predict))
    
    # Recipe with ML steps
    ml_recipe <-
      recipe(Variable ~ ., data = split_train) |>
      update_role(DAid, new_role = "id") |>
      step_normalize(all_numeric_predictors()) |>
      step_nzv(all_numeric_predictors()) |>
      #step_corr(all_numeric_predictors()) |>
      step_impute_knn(all_numeric_predictors())
    
    # LASSO model specifications
    glmnet_specs <-
      linear_reg() |>
      set_mode("regression") |>
      set_engine("glmnet") |>
      set_args(penalty = tune(),
               mixture = 1)
    
    # ML workflow
    glmnet_wflow <-
      workflow() |>
      add_recipe(ml_recipe) |>
      add_model(glmnet_specs)
    
    # Define glmnet grid
    set.seed(213)
    glmnet_grid <-
      glmnet_wflow |>
      extract_parameter_set_dials() |>
      grid_latin_hypercube(size = 20)
    
    # Define the resamples (CV)
    set.seed(213)
    ml_rs <- vfold_cv(split_train, v = 10, strata = Variable)
    
    # Define control_grid
    set.seed(213)
    ctrl <-
      control_grid(save_pred = TRUE, parallel_over = "everything")
    
    cat(paste0("\nFitting glmnet model for ", variable_predict))
    
    # Glmnet grid search
    set.seed(213)
    glmnet_res <-
      glmnet_wflow |>
      tune_grid(resamples = ml_rs,
                grid = glmnet_grid,
                control = ctrl)
    
    predictions_train <-
      glmnet_res |>
      collect_predictions()
    
    metrics_train <-
      glmnet_res |>
      collect_metrics()
    
    cat(paste0("\nSelecting best performing model for ", variable_predict))
    
    # Select best hyperparameter
    best_glmnet <-
      select_best(glmnet_res, metric = "rmse") |>
      select(-.config)
    
    #Finalize the workflow and fit the final model
    glmnet_wflow <-
      glmnet_wflow |>
      finalize_workflow(best_glmnet)
    
    final_glmnet_fit <- last_fit(glmnet_wflow, variable_split)
    
    # Extract model performance
    performance <-
      final_glmnet_fit |>
      collect_metrics() |>
      select(-.config, -.estimator)
    
    glmnet_metrics <-
      final_glmnet_fit |>
      collect_metrics() |>
      filter(.metric == "rmse")
    
    # Extract protein importance
    important_proteins <-
      final_glmnet_fit |>
      extract_fit_parsnip()  |>
      vip::vi(lambda = best_glmnet$penalty, event_level = "second")  |>
      mutate(Importance = abs(Importance),
             Variable = fct_reorder(Variable, Importance))
    
    # Extract model predictions
    predictions <-
      final_glmnet_fit |>
      collect_predictions(summarize = F)
    
    return(
      list(
        "penalty" = best_glmnet,
        "glmnet_model" = glmnet_res,
        "predictions_train" = predictions_train,
        "performance_train" = metrics_train,
        "final_workflow" = glmnet_wflow,
        "final_fit" = final_glmnet_fit,
        "predictions" = predictions,
        "performance" = performance,
        "important_proteins" = important_proteins
      )
    )
  }
generate_split <- function(data,
                           proportion = 0.7,
                           seed = 213,
                           variable_stratify) {
  set.seed(seed)
  
  data_split <-
    data |>
    initial_split(prop = proportion, strata = !!sym(variable_stratify))
  
  data_train <- training(data_split)
  data_test <- testing(data_split)
  
  return(list(
    "data_split" = data_split,
    "data_train" = data_train,
    "data_test" = data_test
  ))
}

do_anova_check <- function(df, protein) {
  
    df_protein <-
    df |>
    filter(Assay == protein)

  res <- sapply(df_protein[c("Age","Sex","BMI","Disease")], function(x) length(unique(x)))
  res |> 
    enframe() |> 
    mutate(Protein = protein) |> 
    relocate(Protein)
}


# Function to run ANOVA for a given protein
do_anova <- function(df, protein) {
  df_protein <-
    df |>
    filter(Assay == protein)
  
  # Linear model with all variables
  model <- lm(NPX ~ Age + Sex + BMI + Disease, data = df_protein)
  
  # Conduct ANOVA
  anova_res <- car::Anova(model, type = 3)
  
  # Calculate Eta-squared using effectsize package
  eta_squared_res <-
    eta_squared(anova_res, partial = TRUE) |>
    as.data.frame()
  
  # Get tidy model results
  tidy_anova <- broom::tidy(anova_res)
  
  # Add Eta-squared to the result
  tidy_anova <-
    tidy_anova |>
    left_join(eta_squared_res, by = c("term" = "Parameter")) |>
    mutate(Eta2_perc = Eta2_partial * 100) |>
    mutate(Protein = protein) |>
    relocate(Protein)
  
  return(tidy_anova)
  
}

plot_variance <- function(anova_res_adjusted, perc_variation = 45) {

  dat <-
  anova_res_adjusted |>
  group_by(Protein) |>
  filter(term != "Residuals") |>
  mutate(Total_variance = sum(Eta2_perc)) |>
  arrange(-Total_variance)

# Proteins with the most explained variance
top_variance <-
  dat |>
  filter(Total_variance > perc_variation) |>
  distinct(Protein) |>
  pull()

# Plot
dat |>
  ggplot(aes(
    x = reorder(Protein,-Total_variance),
    y = Eta2_perc,
    fill = term
  )) +
  geom_col() +
  theme_hpa(angled = T, axis_x = F) +
  scale_fill_manual(values = pal_anova) +
  ylab("Variance explained (%)") +
  xlab("Proteins")
  
}

pal_anova <- c(
  "Age" = "#75C8AE",
  "Sex" = "#EB7C6A",
  "BMI" = "#F7B84D",
  "Disease" = "#EEE2D1"
)
