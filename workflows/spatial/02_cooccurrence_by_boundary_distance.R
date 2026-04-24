# Script: 02_cooccurrence_by_boundary_distance.R
# Purpose: Cell-type abundance changes across boundary-distance bins
# Workflow: spatial
# Required inputs: cell_anno_insitutype.h5ad; cosmx_synora.csv
# Main outputs: combined_data_Tex.csv and distance-stratified plots.
# Prerequisites: Run after Synora boundary detection.

library(tidyverse)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(reticulate)

ad <- import("anndata")
adata <- ad$read_h5ad("cell_anno_insitutype.h5ad")
meta_with_synora <- adata$obs
meta_with_synora$X_index <- rownames(meta_with_synora)
synora <- read.csv("cosmx_synora.csv", row.names = 1)
meta_with_synora <- merge(meta_with_synora, synora, by = "X_index")

global_min <- -400
global_max <- 400
n_intervals <- 50
cell_types <- "c07_CD8_Tex_PDCD1"

assign_distance_intervals <- function(data) {
  data %>%
    mutate(
      Distance_Group = cut(
        Distance2Boundary,
        breaks = breaks,
        include.lowest = TRUE,
        labels = interval_labels
      ),
      Distance_Center = interval_centers[as.integer(Distance_Group)]
    )
}
calculate_celltype_proportions <- function(cell_type_name) {
  cell_type_data <- meta_data %>% 
    filter(cluster_assignment == cell_type_name) %>%
    assign_distance_intervals()
  all_cells_data <- meta_data %>%
    assign_distance_intervals()
  total_counts <- all_cells_data %>%
    group_by(sample_id, Distance_Group, Distance_Center) %>%
    summarise(total_cells = n(), .groups = "drop")
  cell_counts <- cell_type_data %>%
    group_by(sample_id, Distance_Group, Distance_Center) %>%
    summarise(cell_count = n(), .groups = "drop")
  combined <- total_counts %>%
    left_join(cell_counts, by = c("sample_id", "Distance_Group", "Distance_Center")) %>%
    mutate(
      cell_count = replace_na(cell_count, 0),
      Proportion = log2((cell_count + 0.1) / total_cells),
      cell_type = cell_type_name
    ) %>%
    select(sample_id, Distance_Center, Proportion, cell_type)
  return(combined)
}
calculate_bacteria_proportions <- function() {
  bacteria_data <- meta_data %>%
    mutate(has_bacteria = ifelse(Fusobacterium.genus > 0, "Bacteria_positive", "Bacteria_negative")) %>%
    assign_distance_intervals()
  all_cells_data <- meta_data %>%
    assign_distance_intervals()
  total_counts <- all_cells_data %>%
    group_by(sample_id, Distance_Group, Distance_Center) %>%
    summarise(total_cells = n(), .groups = "drop")
  bacteria_counts <- bacteria_data %>%
    filter(has_bacteria == "Bacteria_positive") %>%
    group_by(sample_id, Distance_Group, Distance_Center) %>%
    summarise(bacteria_count = n(), .groups = "drop")
  combined <- total_counts %>%
    left_join(bacteria_counts, by = c("sample_id", "Distance_Group", "Distance_Center")) %>%
    mutate(
      bacteria_count = replace_na(bacteria_count, 0),
      Proportion = log2((bacteria_count + 0.1) / total_cells),
      cell_type = "Fusobacterium_positive"
    ) %>%
    select(sample_id, Distance_Center, Proportion, cell_type)
  return(combined)
}
calculate_zscore_data <- function(proportion_data, cell_type_name) {
  wide_data <- proportion_data %>%
    mutate(
      sample_id = as.character(sample_id),
      interval = as.numeric(cut(Distance_Center, breaks = breaks, include.lowest = TRUE, labels = interval_labels))
    ) %>%
    select(sample_id, interval, Proportion) %>%
    pivot_wider(names_from = interval, values_from = Proportion)
  mat <- as.matrix(wide_data[, -1])
  rownames(mat) <- wide_data$sample_id
  mat_zscore <- t(apply(mat, 1, function(x) {
    if (all(is.na(x)) || sd(x, na.rm = TRUE) == 0) {
      return(rep(0, length(x)))
    }
    (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  }))
  for (j in 1:ncol(mat_zscore)) {
    if (all(is.na(mat_zscore[, j]))) {
      mat_zscore[, j] <- 0
    } else {
      col_mean <- mean(mat_zscore[, j], na.rm = TRUE)
      mat_zscore[is.na(mat_zscore[, j]), j] <- col_mean
    }
  }
  zscore_data <- as.data.frame(mat_zscore)
  zscore_data$sample_id <- rownames(zscore_data)
  zscore_long <- zscore_data %>%
    pivot_longer(
      cols = -sample_id,
      names_to = "interval",
      values_to = "z_score"
    ) %>%
    mutate(
      interval = as.numeric(interval),
      sample_id = as.character(sample_id)
    )
  interval_means <- zscore_long %>%
    group_by(interval) %>%
    summarise(
      mean_z = mean(z_score, na.rm = TRUE),
      se_z = sd(z_score, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(
      distance_center = interval_centers[interval],
      cell_type = cell_type_name
    )
  return(interval_means)
}

breaks <- seq(global_min, global_max, length.out = n_intervals + 1)
interval_centers <- (head(breaks, -1) + tail(breaks, -1)) / 2
interval_labels <- 1:n_intervals
meta_data <- meta_with_synora %>%
  select(
    sample_id,
    Distance2Boundary,
    cluster_assignment,
    Fusobacterium.genus
  ) %>%
  filter(!is.na(Distance2Boundary) & !is.na(cluster_assignment))

all_celltype_results <- map_dfr(cell_types, ~{
  print(paste("celltype:", .x))
  tryCatch({
    calculate_celltype_proportions(.x)
  }, error = function(e) {
    print(paste("celltype", .x, "error:", e$message))
    return(NULL)
  })
})
bacteria_results <- calculate_bacteria_proportions()
all_proportion_data <- bind_rows(all_celltype_results, bacteria_results)
target_cell_types <- c("c07_CD8_Tex_PDCD1", "c41_mregDC_LAMP3")
zscore_results <- list()
for (cell_type in target_cell_types) {
  cell_data <- all_proportion_data %>% filter(cell_type == !!cell_type)
  if (nrow(cell_data) > 0) {
    zscore_results[[cell_type]] <- calculate_zscore_data(cell_data, cell_type)
  }
}
bacteria_zscore <- calculate_zscore_data(
  all_proportion_data %>% filter(cell_type == "Fusobacterium_positive"),
  "Fusobacterium_positive"
)
combined_zscore_data <- bind_rows(
  zscore_results$c07_CD8_Tex_PDCD1,
  zscore_results$c41_mregDC_LAMP3,
  bacteria_zscore
)
write.csv(combined_zscore_data, "combined_data_Tex.csv")

