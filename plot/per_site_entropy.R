
library(readr)
library(dplyr)
library(purrr)
library(tools)
library(ggplot2)
library(patchwork)

# Define general-purpose functions

read_entropy_data <- function(base_path) {
  # Get first-level subdirectories, each folder name represents a model name）
  model_folders <- list.dirs(base_path, full.names = TRUE, recursive = FALSE)

  map_dfr(model_folders, function(model_path) {
    model_name <- basename(model_path)

    # Construct full paths for 100 files under the current model folder
    file_paths <- file.path(model_path, paste0("seq_", 1:100, ".fa_gap.sitewise_entropy.txt"))

    map_dfr(file_paths, function(file) {
      if (!file.exists(file)) return(NULL)

      df <- read_tsv(file, col_types = cols())
      df %>%
        select(entropy) %>%
        mutate(
          model = model_name,
          file  = file_path_sans_ext(basename(file)),
          .before = 1
        )
    })
  })
}

plot_entropy_density <- function(entropy_table, empirical_path, model_name) {

  # Read empirical data
  empirical_df <- read_tsv(empirical_path, col_types = cols())

  # Filter the target model data from the full table
  model_entropy <- entropy_table %>% filter(model == model_name)

  # Plotting
  p <- ggplot() +
    geom_density(
      data = model_entropy,
      aes(x = entropy, group = file),
      color = "black", alpha = 0.15, size = 0.4
    ) +
    geom_density(
      data = empirical_df,
      aes(x = entropy),
      color = "red", size = 0.4
    ) +
    labs(
      title = model_name,
      x = "Entropy",
      y = "Density"
    ) +
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

example_ent_table  <- read_entropy_data("D:/8701/August/data/model_simulation")

# head(example_ent_table)


p1 <- plot_entropy_density(
  model_name     = "LG+R3",
  entropy_table  = example_ent_table,
  empirical_path = "D:/8701/August/data/model_simulation/original.sitewise_entropy.txt"
)

p1 <- plot_entropy_density(
  model_name     = "LG+R3",
  entropy_table  = example_ent_table,
  empirical_path = "D:/8701/August/data/model_simulation/original.sitewise_entropy.txt"
)

p2 <- plot_entropy_density(
  model_name     = "LG+C60+R3",
  entropy_table  = example_ent_table,
  empirical_path = "D:/8701/August/data/model_simulation/original.sitewise_entropy.txt"
)

p3 <- plot_entropy_density(
  model_name     = "LG+C60fix+R3-PMSF",
  entropy_table  = example_ent_table,
  empirical_path = "D:/8701/August/data/model_simulation/original.sitewise_entropy.txt"
)

p4 <- plot_entropy_density(
  model_name     = "LG+C60opt+R3-PMSF",
  entropy_table  = example_ent_table,
  empirical_path = "D:/8701/August/data/model_simulation/original.sitewise_entropy.txt"
)

p5 <- plot_entropy_density(
  model_name     = "GTR20+C60fix+R3-PMSF",
  entropy_table  = example_ent_table,
  empirical_path = "D:/8701/August/data/model_simulation/original.sitewise_entropy.txt"
)

p6 <- plot_entropy_density(
  model_name     = "GTR20+C60opt+R3-PMSF",
  entropy_table  = example_ent_table,
  empirical_path = "D:/8701/August/data/model_simulation/original.sitewise_entropy.txt"
)


final_plot <- (p1 + p2 + p3) / (p4 + p5 + p6)

print(final_plot)

#ggsave("D:/8701/New-Metrics-PBT/plot/example_per_site_entropy.png", plot = final_plot, width = 8, height = 6, dpi = 900)
