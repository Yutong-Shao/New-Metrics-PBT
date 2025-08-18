
library(readr)
library(dplyr)
library(ggplot2)
library(ggridges)
library(patchwork)

# Step 1: Read the txt file
cvm_results <- readr::read_tsv("D:/8701/New-Metrics-PBT/data/model_results/cvm_entropy_results.txt") %>%
  mutate(
    group = case_when(
      model == "LG+R3" ~ "LG",
      model == "LG+C60+R3" ~ "C60",
      model %in% c("LG+C60fix+R3-PMSF", "LG+C60opt+R3-PMSF", "GTR20+C60fix+R3-PMSF", "GTR20+C60opt+R3-PMSF") ~ "PMSF",
      TRUE ~ NA_character_
    ),
    weights = case_when(
      grepl("fix", model, ignore.case = TRUE) ~ "fix",
      grepl("opt", model, ignore.case = TRUE) ~ "opt",
      TRUE ~ "default"
    ),
    exchangeability = case_when(
      grepl("Poisson", model) ~ "Poisson",
      grepl("LG", model) ~ "LG",
      grepl("GTR20", model) ~ "GTR20",
      TRUE ~ NA_character_
    ),
    label_raw = model %>%
      gsub("fix", "", .) %>%
      gsub("opt", "", .)
  ) %>%
  mutate(
    sort_key = case_when(
      model == "LG+R3" ~ 1,
      model == "LG+C60+R3" ~ 2,
      model %in% c("LG+C60fix+R3-PMSF", "LG+C60opt+R3-PMSF", "GTR20+C60fix+R3-PMSF", "GTR20+C60opt+R3-PMSF") ~ 3,
      TRUE ~ 99
    )
  ) %>%
  arrange(group, exchangeability, sort_key, label_raw)

# Set label levels
cvm_results <- cvm_results %>%
  mutate(label = factor(label_raw, levels = rev(unique(label_raw))))

# Step 2: Prepare label_positions
labels <- levels(cvm_results$label)
label_positions <- data.frame(
  label = labels,
  ypos  = seq_along(labels)
)

# Step 3: Ridge Plot
p_ridges <- ggplot(cvm_results, aes(x = W2_statistic, y = label, fill = weights)) +
  geom_density_ridges(
    scale = 1.3, rel_min_height = 0, alpha = 0.7,
    color = "black", linewidth = 0.37,
    panel_scaling = FALSE, position = "identity"
  ) +
  geom_segment(
    data = label_positions,
    aes(x = min(cvm_results$W2_statistic), xend = max(cvm_results$W2_statistic),
        y = ypos, yend = ypos),
    inherit.aes = FALSE, color = "black", linewidth = 0.4
  ) +
  scale_x_continuous(expand = c(0,0)) +
  scale_y_discrete(expand = expansion(add = c(0,0))) +
  labs(x = expression(W^2 ~ "(Cramér-von Mises statistic) on sitewise entropy")) +
  scale_fill_manual(values = c("fix" = "#9AC4DB", "opt" = "#CB8D8B", "default" = "grey70")) +
  theme_void() +
  theme(
    axis.title.x = element_text(size = 9),
    axis.text.x = element_text(size = 8),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_line(color = "black", size = 0.5),
    axis.line.x = element_line(color = "black", size = 0.5),
    axis.ticks.x = element_line(color = "black", size = 0.5),
    axis.ticks.y.left = element_line(color = "black", size = 0.5),
    axis.ticks.length = unit(0.1, "cm"),
    panel.grid.major.x = element_blank(),
    legend.position = "top",
    legend.margin = margin(b = 5),
    legend.title = element_blank(),
    legend.text = element_text(size = 7)
  )

# Step 4: Labels plot
label_background <- cvm_results %>%
  select(label, group) %>%
  distinct() %>%
  mutate(ypos = seq(from = 1, by = 2, length.out = n()))

p_labels <- ggplot(label_background, aes(x = 2, y = ypos, fill = group)) +
  geom_tile(width = 3.7, height = 2.0) +
  scale_fill_manual(values = c(
    "LG" = "#E9CDDF",
    "C60" = "#C8BFD9",
    "PMSF" = "#B6C9C0"
  )) +
  geom_text(aes(label = label, y = ypos),
            hjust = 1, size = 2.5, nudge_x = 1.8) +
  scale_y_reverse(expand = expansion(add = c(-1,1))) +
  theme_void() +
  theme(legend.position = "none") +
  coord_cartesian(clip = "off", xlim = c(0.3,3.58)) +
  theme(plot.margin = margin(0.5,0.14,0,-0))

# Step 5: Combine the plots
final_plot <- p_labels + plot_spacer() + p_ridges +
  plot_layout(widths = c(0.8, 0.000001, 4))

print(final_plot)

ggsave("D:/8701/New-Metrics-PBT/plot/example_cvm_entropy.png", plot = final_plot, width = 8, height = 5, dpi = 900)
