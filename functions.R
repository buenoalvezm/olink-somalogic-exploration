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
