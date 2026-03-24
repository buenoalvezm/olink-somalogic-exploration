# Functions

subset_assays <- function(data, mode) {
  core_cols <- c("DAid", "Age")

  if (mode == "olink") {
    data |>
      dplyr::select(
        -dplyr::starts_with("seq.")
      )

  } else if (mode == "soma") {
    data |>
      dplyr::select(
        dplyr::all_of(core_cols),
        dplyr::starts_with("seq.")
      )

  } else if (mode == "combined") {
    data

  } else {
    stop("Invalid run_mode")
  }
}


#Palette
library(RColorBrewer)
getPalette3 = colorRampPalette(brewer.pal(8, "Set2"))
pal_class<-getPalette3(8)
names(pal_class)<-c("Neurological","Cardiovascular","Cancer","Autoimmune","Pediatric","Infection","Metabolic","Healthy") 
class_order <- c("Healthy", "Cardiovascular","Metabolic","Cancer","Neurological","Autoimmune","Infection","Pediatric")

# HPA theme
theme_hpa <- 
  function(angled = F, axis_x = T, axis_y = T, facet_title = T) {
    t <- 
      theme(
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.spacing = unit(0.2, "lines"), 
        panel.background=element_rect(fill="white"),
        panel.border = element_blank(),
        plot.title = element_text(face = "bold",
                                  size = rel(1), hjust = 0.5),
        plot.subtitle=element_text(face = "bold",hjust = 0.5, size=rel(1),vjust=1),
        axis.title = element_text(face = "bold",size = rel(1)),
        axis.ticks.length = unit(.25, "cm"),
        axis.line = element_line(linewidth = 0.5),
        axis.text = element_text(size = rel(1), color = 'black'),
        legend.key = element_blank(),
        legend.position = "right",
        legend.text = element_text(size=rel(0.8)),
        legend.key.size= unit(0.7, "cm"),
        legend.title = element_text(size=rel(1)),
        plot.margin=unit(c(10,5,5,5),"mm"),
        strip.background=element_rect(colour="grey90",fill="grey90"),
        strip.text = element_text(face="bold")
      )
    
    if(angled) {
      t <- t + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) 
    }
    
    if(axis_x == F) {
      t <- t +
        theme(axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              axis.line.x = element_blank(),
              axis.title.x = element_blank())
    } 
    
    if(axis_y == F) {
      t <- t +
        theme(axis.text.y = element_blank(),
              axis.ticks.y = element_blank(),
              axis.line.y = element_blank(),
              axis.title.y = element_blank())
    }
    if(facet_title == F) {
      t <- t + theme(strip.text = element_blank())
    }
    return(t)
  }
# Function to generate splits
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
  
  return(list("data_split" = data_split, 
              "data_train" = data_train,
              "data_test" = data_test))
  
}



# Function to run multiclassification against all other diseases
disease_against_all<-  
  function(split_train, 
           split_test,
           seed,
          path_save) {                    
  
    
    disease_split <- make_splits(split_train, split_test)
    
    
    ### Define general recipe
    ml_recipe <- 
      recipe(Disease ~ ., data = split_train) |> 
      update_role(DAid, new_role = "id") |> 
      step_normalize(all_numeric()) |> 
      step_nzv(all_numeric()) |> 
      #step_corr(all_numeric()) |> 
      step_impute_knn(all_numeric()) 
    
    # Generate resamples
    set.seed(213)
    disease_rs <- vfold_cv(split_train, v = 10, strata = Disease)
    
    # Define evaluation metrics for all workflows
    eval_metrics <- metric_set(roc_auc)
    
    # Define control grid
    set.seed(213)
    ctrl <- control_grid(verbose = TRUE, 
                         allow_par = TRUE,
                         save_pred = TRUE, 
                         parallel_over = "everything") 
    
    # Tidymodels lasso multiclassification recipe
    glmnet_lasso_specs <-
      multinom_reg() |>
      set_mode("classification") |>
      set_engine("glmnet") |>
      set_args(penalty = tune(),
               mixture = 1)
    
    # Set up lasso workflow
    glmnet_wflow <-
      workflow() |> 
      add_recipe(ml_recipe) |> 
      add_model(glmnet_lasso_specs) 
    
    # Define hyperparameter tuning grid
    set.seed(213)
    glmnet_grid <-
      glmnet_wflow |>
      extract_parameter_set_dials() |>
      grid_latin_hypercube(size = 20)
    
    # Hyperparameter tuning
    set.seed(213)
    glmnet_res <-
      glmnet_wflow |>
      tune_grid(
        resamples = disease_rs,
        grid = glmnet_grid,
        control = ctrl,
        metrics = eval_metrics)
    
    # Explore results and select the best performing hyperparameter combination
    autoplot(glmnet_res)
    
    predictions_train <- 
      glmnet_res |> 
      collect_predictions()
    
    metrics_train <- 
      glmnet_res |> 
      collect_metrics()
    
    best_glmnet <- 
      glmnet_res |> 
      select_best(metric = "roc_auc")
    
    # Final fit
    final_glmnet <- 
      glmnet_wflow |> 
      finalize_workflow(best_glmnet)
    
    final_glmnet_fit <- 
      last_fit(final_glmnet, disease_split)
    
    # Explore performance
    final_predictions <- 
      final_glmnet_fit |> 
      collect_predictions() |> 
      mutate(DAid = split_test$DAid) |> 
      relocate(DAid) |> 
      select(-id)
    
  write_tsv(final_predictions, paste(path_save, "predictions_multiclass_seed_", seed, ".tsv", sep = ""))
    
    final_metrics <- 
      final_glmnet_fit |> 
      collect_metrics()
    
    # ROC
    dat <-
      split_test %>% 
      select(DAid, Disease) %>% 
      mutate(value = 1) %>% 
      spread(Disease,value, fill= 0) 
    
    true_dat <- 
      dat %>% 
      set_names(paste(names(dat), "_true", sep = "")) %>%
      rename(DAid = `DAid_true`)
    
    dat_prob <- 
      final_predictions %>% 
      rename_all(~stringr::str_replace_all(.,".pred_","")) |> 
      select(-".row", -"class", -"Disease", -".config") 
    
    prob_data <- 
      dat_prob %>% 
      set_names(paste(names(dat_prob), "_pred_glmnet", sep = "")) %>% 
      rename(DAid = DAid_pred_glmnet)
    
    final_df <- 
      true_dat %>% 
      left_join(prob_data, by = "DAid") %>% 
      select(-DAid) |> 
      as.data.frame()
    
    roc_res <- multi_roc(final_df, force_diag=T)
    plot_roc_df <- plot_roc_data(roc_res)
    
#    auc_ci <- roc_auc_with_ci(final_df)
    
    roc_dat <- 
      plot_roc_df %>%
      filter(!Group %in% c("Macro","Micro")) %>% 
      mutate(Performance = paste(Group, ": ", round(AUC, 4), sep = "")) %>% 
      arrange(-AUC)
    
 #   roc_dat_ext <- 
  #    roc_dat |> 
   #   left_join(auc_ci |> 
    #              mutate(Group = gsub("glmnet.", "", Var)) |> 
     #             rename(AUC_CI = AUC), by = "Group")
    
    write_tsv(roc_dat, paste(path_save, "roc_multiclass_seed_", seed, ".tsv", sep = "")) 
   # write_tsv(final_predictions, paste(path_save, "/predictions_multiclass_seed_", seed, ".tsv", sep = ""))
    
    roc_plot <- 
      roc_dat %>% 
      left_join(resource_meta |> 
                  distinct(Disease, Class), by = c("Group" = "Disease")) |> 
      # mutate(Group = factor(Group, levels = disease_class_order_roc)) |> 
      ggplot(aes(x = 1-Specificity, y=Sensitivity, color= Class)) +
      geom_path(size=1, show.legend = F) +
      geom_segment(aes(x = 0, y = 0, xend = 1, yend = 1), 
                   color = "grey",
                   linetype = 'dotdash',
                   show.legend = F) +
      geom_text(aes(label = round(AUC, 2), x = 0.75, y = 0.25), show.legend = F) +
      theme_hpa() +
      scale_y_continuous(breaks = c(0, 1)) +
      scale_x_continuous(breaks = c(0, 1)) +
      scale_color_manual(values = pal_class) +
      facet_wrap(~Group, nrow = 5) +
      coord_fixed()
    
    # Confusion matrix
    cm <- 
      final_glmnet_fit %>% 
      collect_predictions() %>%
      conf_mat(truth = Disease, estimate = .pred_class) 
    
    cm_plot <- 
      cm |>
      autoplot(type = "heatmap") +
      theme_hpa(angled = T)
 
    # Important proteins
    important_proteins <- 
      final_glmnet_fit  |> 
      extract_fit_parsnip() %>%
      tidy() |> 
      filter(term != "(Intercept)") |> 
      arrange(-abs(estimate)) |> 
      filter(abs(estimate) > 0)

     write_tsv(important_proteins, paste(path_save, "important_proteins_multiclass_seed_", seed, ".tsv", sep = ""))
    
    return(list("penalty" = best_glmnet,
                "glmnet_model" = glmnet_res,
                "predictions_train" = predictions_train, 
                "performance_train" = metrics_train,
                "final_workflow" = final_glmnet,
                "final_fit" = final_glmnet_fit,
                "predictions" = final_predictions,
                "final_metrics" = final_metrics,
                #"performance" = performance,
                "confusion_matrix" = cm,
                "confusion_matrix_plot" = cm_plot,
                "roc_curve" = roc_dat,
                "roc_plot" = roc_plot,
                "important_proteins" = important_proteins))
    
  }
