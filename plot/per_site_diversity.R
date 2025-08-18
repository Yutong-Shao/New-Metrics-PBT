
library(readr)
library(dplyr)
library(tidyr)
library(purrr)
library(tools)
library(ggplot2)
library(patchwork)

# Define general-purpose functions

# 1. process_dataset_diversity(): Traverse all model folders under base_path and count the frequency of diversity values in each seq_i
process_dataset_diversity <- function(base_path) {
  # Get first-level subdirectories, each folder name represents a model name
  model_folders <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)

  map_dfr(model_folders, function(model_path) {
    model_name <- basename(model_path)
    # Construct full paths for 100 files under the current model folder
    file_paths <- file.path(model_path, paste0("seq_", 1:100, ".fa_gap.sitewise_diversity.txt"))

    map_dfr(file_paths, function(file) {
      if (!file.exists(file)) return(NULL)
      df <- read_tsv(file, col_types = cols())
      # Count occurrences of diversity values from 1 to 20
      df %>%
        count(diversity) %>%
        complete(diversity = 1:20, fill = list(n = 0)) %>%
        mutate(
          model   = model_name,
          file    = file_path_sans_ext(basename(file)),
          .before = 1
        ) %>%
        rename(count = n)
    })
  })
}

# 2. plot_model_with_empirical(): General plotting function
plot_model_with_empirical <- function(model_name, model_table, empirical_path, dataset_label) {
  # Filter the data of the specified model
  model_data <- model_table %>%
    filter(model == model_name)

  # Construct a list of smoothed spline curves
  spline_list <- model_data %>%
    group_by(file) %>%
    arrange(diversity) %>%
    group_split() %>%
    map_dfr(function(df) {
      spline_df <- as.data.frame(spline(x = df$diversity, y = df$count, n = 200))
      spline_df$file <- unique(df$file)
      return(spline_df)
    })

  # Read empirical data
  empirical_data <- read_tsv(empirical_path, col_types = cols())

  empirical_freq <- empirical_data %>%
    count(diversity) %>%
    complete(diversity = 1:20, fill = list(n = 0)) %>%
    arrange(diversity)

  empirical_spline <- as.data.frame(
    spline(x = empirical_freq$diversity, y = empirical_freq$n, n = 300)
  )

  # Plotting
  p <- ggplot() +
    geom_density(data = model_entropy, aes(x = entropy, group = file), color = "black", alpha = 0.15, size = 0.4) +
    geom_density(data = empirical_df, aes(x = entropy), color = "red", size = 0.4) +
    labs(title = model_name, x = "Entropy", y = "Density") +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 9),
      axis.text.x  = element_text(size = 6),
      axis.text.y  = element_text(size = 6),
      axis.title.y = element_text(size = 8),
      axis.title.x = element_text(size = 8)
    )

  return(p)
}

# Apply the function to each dataset
example_div_table  <- process_dataset_diversity("D:/8701/August/data/model_simulation")

# View sample output
# head(example_div_table)

# Plotting

p1 <- plot_model_with_empirical(
  model_name = "LG+R3",
  model_table = example_div_table,
  empirical_path = "D:/8701/August/data/model_simulation/original.sitewise_diversity.txt"
)

p2 <- plot_model_with_empirical(
  model_name = "LG+C60+R3",
  model_table = example_div_table,
  empirical_path = "D:/8701/August/data/model_simulation/original.sitewise_diversity.txt"
)

p3 <- plot_model_with_empirical(
  model_name = "LG+C60fix+R3-PMSF",
  model_table = example_div_table,
  empirical_path = "D:/8701/August/data/model_simulation/original.sitewise_diversity.txt"
)

p4 <- plot_model_with_empirical(
  model_name = "LG+C60opt+R3-PMSF",
  model_table = example_div_table,
  empirical_path = "D:/8701/August/data/model_simulation/original.sitewise_diversity.txt"
)

p5 <- plot_model_with_empirical(
  model_name = "GTR20+C60fix+R3-PMSF",
  model_table = example_div_table,
  empirical_path = "D:/8701/August/data/model_simulation/original.sitewise_diversity.txt"
)

p6 <- plot_model_with_empirical(
  model_name = "GTR20+C60opt+R3-PMSF",
  model_table = example_div_table,
  empirical_path = "D:/8701/August/data/model_simulation/original.sitewise_diversity.txt"
)

final_plot <- (p1 + p2 + p3) / (p4 + p5 + p6)

print(final_plot)

#ggsave("D:/8701/New-Metrics-PBT/plot/example_per_site_diversity.png", plot = final_plot, width = 8, height = 6, dpi = 900)
