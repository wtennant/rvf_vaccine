# plot_compare_opt_island_balance.R: plotting to compare the balance of strategy effectiveness of
# different optimal vaccine strategies on each island in the archipelago.
# Pre-requisites: optimal parameters already obtained for each vaccination rate and tagging strategy,
#                 optimal parameters fed through simulation under corresponding scenarios,
#                 simulations have already been summarised using data_summarise_scenario.R.

# Clear the workspace.
rm(list = ls())

# Load in the data visualisation and manipulation libraries.
library(sf)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(data.table)
library(colorspace)
library(ggokabeito)
library(grid)

# Load in our favourite ggplot theme.
source("gg_theme.R")

# Define the optimal vaccine strategies to compare.
# First character (identifiability):
#   (U)ntagged animals, or
#   (T)agged animals.
# Second character (vaccine distribution):
#   (P)roportional to population size,
#   (1) Optimised over infections averted across the archipelago, or
#   (2) Optimised over infections averted on Worst-performing island.
strat_ids <- c("U1", "U2", "T1", "T2")

# Define the folder locations of all optimisations/scenarios as a named vector.
in_folders <- c(
  paste0(
    rep("../cpp/out/", length(strat_ids)),
    strat_ids
  )
)
names(in_folders) <- strat_ids

# Define the identifiers and labels of strategy options.
tag <- data.frame(
  TAG_ID = factor(c("U", "T"), levels = c("U", "T")),
  TAG_LABEL = c("Untagged", "Tagged"),
  LIGHTEN = c(0.4, -0.4)
)
dist <- data.frame(
  DIST_ID = factor(c("P", "1", "2"), levels = c("P", "1", "2")),
  DIST_LABEL = c(
    "Proportional to island population size",
    "Optimised over infections averted across the archipelago",
    "Optimised over infections averted on the worst-peforming island"
  ),
  COLOUR = c("#777777", palette_okabe_ito(order = c(1, 2)))
)

# Define the identifiers and labels for all combinations of strategies.
all_strats <- bind_rows(tag, dist) %>%
  tidyr::expand(
    nesting(DIST_ID, DIST_LABEL, COLOUR),
    nesting(TAG_ID, TAG_LABEL, LIGHTEN)
  ) %>%
  drop_na(everything()) %>%
  unite(col = "STRAT_ID", TAG_ID, DIST_ID, sep = "", remove = FALSE)

# Calculate the total population of the Comoros archipelago.
total_pop <- fread(paste0(in_folders[1], "/simulations/default_pars_000.csv"), header = TRUE) %>%
  filter(PAR_NAME == "n_pop") %>%
  .$PAR_VALUE %>%
  sum()

# Define the settings for the main parameter of interest in common between all optimisations/scenarios.
sens_par_name_1 <- "vac_rate"
trans_sens_par_1 <- \(x) {
  round(x / total_pop, digits = 2)
}
sens_par_title_1 <- "Percentage of livestock vaccinated annually across\nthe Comoros archipelago"
sens_par_label_1 <- scales::percent
sens_par_1_filter <- c("5%", "30%")

# Data frame to store all summary information.
summary_df <- data.frame()

# Initialise a counter for ID numbers.
i <- j <- 0

# For each opt/scenario...
for (in_folder in in_folders) {
  # Get the number of unique summary files.
  file_ids <- gsub("[A-z.]", "", list.files(path = paste0(in_folder, "/simulations"), pattern = "summary_[0-9]+.csv"))

  # For each file...
  for (file_id in file_ids) {
    # Load in the scenario parameters.
    scen_par_df <- fread(paste0(in_folder, "/simulations/scenario_pars_", file_id, ".csv"), header = TRUE)

    # Calculate the number of indices per scenario parameter, combine the parameter
    # name and index columns, and then expand the data frame to a wide format.
    scen_par_df <- scen_par_df %>%
      mutate(PAR_NAME = ifelse(PAR_NAME == sens_par_name_1, "SENS_PAR_1", PAR_NAME)) %>%
      mutate(PAR_VALUE = ifelse(PAR_NAME == "SENS_PAR_1", trans_sens_par_1(PAR_VALUE), PAR_VALUE)) %>%
      filter(PAR_NAME == "SENS_PAR_1") %>%
      pivot_wider(id_cols = SAMPLE_ID, names_from = PAR_NAME, values_from = PAR_VALUE)

    # Only proceed if there are valid values for the sensitivity parameter.
    if (nrow(scen_par_df) > 0) {
      # Load in the summary data.
      one_summary_df <- fread(paste0(in_folder, "/simulations/summary_", file_id, ".csv"), header = TRUE)

      # Join into summary data frame.
      one_summary_df <- right_join(one_summary_df, scen_par_df, by = "SAMPLE_ID")

      # Note the input scenario / optimisation name.
      one_summary_df$STRAT_ID <- names(in_folders)[which(in_folder == in_folders)]
      one_summary_df$GROUP_ID <- i

      # Move onto the next id number.
      i <- i + 1

      # Append to the overall data frame.
      summary_df <- bind_rows(summary_df, one_summary_df)
    }
  }
}

# Set the indices of the optimised parameters to factors.
summary_df <- summary_df %>%
  mutate(
    ISLAND_ID = factor(ISLAND_ID, levels = c(1, 3, 0, 2, -1)),
    SENS_PAR_1 = factor(sens_par_label_1(SENS_PAR_1))
  ) %>%
  left_join(
    all_strats,
    by = "STRAT_ID"
  )
summary_df$SENS_PAR_1 <- factor(summary_df$SENS_PAR_1,
  levels = paste0(sort(as.numeric(gsub("%", "", levels(summary_df$SENS_PAR_1)))), "%")
)

# Per sample, calculate the mean effectiveness across all islands.
summary_df <- summary_df %>%
  filter(ISLAND_ID != -1) %>%
  mutate(
    DIFF_EFF_MEAN = EFFECTIVENESS - mean(EFFECTIVENESS),
    .by = c(SAMPLE_ID, REP_ID, SENS_PAR_1:LIGHTEN)
  )

# Summarise each scenario in two ways:
# the ef
stat_summary_df <- summary_df %>%
  summarise(
    EFF_Q025 = quantile(EFFECTIVENESS, probs = 0.025),
    EFF_Q500 = quantile(EFFECTIVENESS, probs = 0.5),
    EFF_Q975 = quantile(EFFECTIVENESS, probs = 0.975),
    DIFF_EFF_Q025 = quantile(DIFF_EFF_MEAN, probs = 0.025),
    DIFF_EFF_Q500 = quantile(DIFF_EFF_MEAN, probs = 0.5),
    DIFF_EFF_Q975 = quantile(DIFF_EFF_MEAN, probs = 0.975),
    EFF_SD = sd(EFFECTIVENESS),
    DIFF_EFF_SD = sd(DIFF_EFF_MEAN),
    .by = c(ISLAND_ID, SENS_PAR_1:LIGHTEN)
  ) %>%
  mutate(
    EFF_LABEL = paste0(
      100 * round(EFF_Q500, digits = 2), "%\n[", 100 * round(EFF_Q025, digits = 2),
      "%, ", 100 * round(EFF_Q975, digits = 2), "%]"
    ),
    DIFF_EFF_LABEL = paste0(
      100 * round(DIFF_EFF_Q500, digits = 2), "%\n[", 100 * round(DIFF_EFF_Q025, digits = 2),
      "%, ", 100 * round(DIFF_EFF_Q975, digits = 2), "%]"
    )
  )

# Map island indices to names.
island_labels <- c(
  "0" = "Anjouan",
  "1" = "Grande Comore",
  "2" = "Mayotte",
  "3" = "Mohéli",
  "-1" = "Comoros archipelago"
)

# Load in the shape files and manipulate the shapefiles into
# a single sf object.
comoros <- read_sf("../data/shapefiles/com_admbnda_adm1_cosep_ocha_20191205.shp") %>%
  select(ADM1_EN, geometry) %>%
  mutate(ADM1_EN = gsub(" \\([A-z]+\\)", "", ADM1_EN))
mayotte <- read_sf("../data/shapefiles/MYT_adm1.shp")
mayotte <- st_as_sf(st_union(mayotte)) %>%
  rename(geometry = x) %>%
  mutate(ADM1_EN = "Mayotte")
archipelago <- rbind(comoros, mayotte) %>%
  rowwise() %>%
  mutate(
    ISLAND_ID = factor(
      setNames(names(island_labels), island_labels)[ADM1_EN],
      levels = c(1, 3, 0, 2)
    )
  )

# Join in the summary data.
archi_data <- left_join(archipelago, stat_summary_df, by = "ISLAND_ID")

# Define functions to re-position the statistics on the map.
label_nudge <- function(geometry, island) {
  x <- suppressWarnings(st_point_on_surface(st_zm(geometry)))
  if (island == 0) {
    x <- x + c(0.5, 0)
  } else if (island == 1) {
    x <- x + c(0.5, 0.1)
  } else if (island == 2) {
    x <- x - c(0.5, 0)
  } else if (island == 3) {
    x <- x - c(0, 0.3)
  }
  x <- st_set_crs(x, st_crs(geometry))
  return(x)
}

# Use the function to define points of where the stats
# should be on the map.
archi_data <- archi_data %>%
  rowwise() %>%
  mutate(
    geometry_label = label_nudge(geometry, ISLAND_ID)
  )

# Define the way of labelling the sensitivity parameter panels.
sens_par_labels <- paste0(sens_par_title_1, ": ", unique(stat_summary_df$SENS_PAR_1))
names(sens_par_labels) <- unique(stat_summary_df$SENS_PAR_1)

# Plot the effectiveness per island.
gg_eff <- archi_data %>%
  filter(SENS_PAR_1 %in% sens_par_1_filter) %>%
  ggplot() +
  geom_sf(
    mapping = aes(
      geometry = geometry,
      fill = EFF_Q500,
    ),
    colour = "black",
    linewidth = 0.75
  ) +
  geom_sf_text(
    mapping = aes(
      geometry = geometry_label,
      label = EFF_LABEL
    ),
    colour = "black",
    fun.geometry = identity,
    size = 4
  ) +
  scale_fill_distiller(
    name = "Median percentage of infections averted on the island",
    direction = 1,
    palette = "YlOrRd",
    breaks = seq(0, 1, by = 0.1),
    limits = c(0, 1),
    labels = scales::percent
  ) +
  facet_grid(
    DIST_LABEL + TAG_LABEL ~ SENS_PAR_1,
    labeller = labeller(
      SENS_PAR_1 = as_labeller(sens_par_labels, default = label_wrap_gen(width = 30)),
      DIST_LABEL = label_wrap_gen(width = 30)
    )
  ) +
  gg_theme +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.key.width = unit(25, "mm")
  ) +
  guides(
    fill = guide_colorbar(
      show.limits = TRUE,
      frame.colour = "black",
      ticks.colour = "black",
      title.position = "top"
    )
  )

# Build the figure and edit the legend to hide legend lines / figures for summary statistics.
gg_eff <- ggplot_gtable(ggplot_build(gg_eff))

# Also edit the colour of the strips.
strip_ids <- grep("strip-r", gg_eff$layout$name)
for (strip_id in strip_ids){
  # In each strip, the first strip text is associated with
  # the tagged status of livestock.
  tag_label <- gg_eff$grobs[[strip_id]]$grobs[[1]]$children[[2]]$children[[1]]$label
  lighten <- tag$LIGHTEN[tag$TAG_LABEL == tag_label]
  gg_eff$grobs[[strip_id]]$grobs[[1]]$children[[1]]$gp$fill <- colorspace::lighten(
    "grey85", amount = lighten / 2
  )

  # Now change the colour for the vaccine distribution.
  dist_label <- gg_eff$grobs[[strip_id]]$grobs[[2]]$children[[2]]$children[[1]]$label
  colour <- dist$COLOUR[dist$DIST_LABEL == str_wrap(dist_label, width = 100)]
  gg_eff$grobs[[strip_id]]$grobs[[2]]$children[[1]]$gp$fill <- colorspace::lighten(
    colour,
    amount = 0.5
  )
}

# Save the figure: effectiveness on each island.
file_name <- paste(
  "fig_compare_opt_effect_islandalt",
  paste(strat_ids, collapse = "_"),
  paste(gsub("%", "", sens_par_1_filter), collapse = "_"),
  sep = "_"
)
ggsave(
  filename = paste0(file_name, ".svg"), plot = gg_eff, width = 9, height = 11, scale = 1.1,
  path = "figures"
)
ggsave(
  filename = paste0(file_name, ".pdf"), plot = gg_eff, width = 9, height = 11, scale = 1.1,
  path = "figures"
)


# Plot the effectiveness per island.
gg_diff_eff <- archi_data %>%
  filter(SENS_PAR_1 %in% sens_par_1_filter) %>%
  ggplot() +
  geom_sf(
    mapping = aes(
      geometry = geometry,
      fill = DIFF_EFF_Q500,
    ),
    colour = "black",
    linewidth = 0.75
  ) +
  geom_sf_text(
    mapping = aes(
      geometry = geometry_label,
      label = DIFF_EFF_LABEL
    ),
    colour = "black",
    fun.geometry = identity,
    size = 4
  ) +
  scale_fill_fermenter(
    name = paste(
      "Difference in the percentage of infections averted",
      "between each island and the mean across the archipelago",
      sep = "\n"
    ),
    palette = "RdBu",
    breaks = seq(-0.45, 0.45, by = 0.1),
    limits = c(-0.45, 0.45),
    labels = scales::percent
  ) +
  facet_grid(
    DIST_LABEL + TAG_LABEL ~ SENS_PAR_1,
    labeller = labeller(
      SENS_PAR_1 = as_labeller(sens_par_labels, default = label_wrap_gen(width = 30)),
      DIST_LABEL = label_wrap_gen(width = 30)
    )
  ) +
  gg_theme +
  theme(panel.border = element_rect(colour = "black", fill = NA)) +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    legend.key.width = unit(25, "mm")
  ) +
  guides(
    fill = guide_colorsteps(
      show.limits = TRUE,
      frame.colour = "black",
      title.position = "top"
    )
  )

# Build the figure and edit the legend to hide legend lines / figures for summary statistics.
gg_diff_eff <- ggplot_gtable(ggplot_build(gg_diff_eff))

# Also edit the colour of the strips.
strip_ids <- grep("strip-r", gg_diff_eff$layout$name)
for (strip_id in strip_ids){
  # In each strip, the first strip text is associated with
  # the tagged status of livestock.
  tag_label <- gg_diff_eff$grobs[[strip_id]]$grobs[[1]]$children[[2]]$children[[1]]$label
  lighten <- tag$LIGHTEN[tag$TAG_LABEL == tag_label]
  gg_diff_eff$grobs[[strip_id]]$grobs[[1]]$children[[1]]$gp$fill <- colorspace::lighten(
    "grey85", amount = lighten / 2
  )

  # Now change the colour for the vaccine distribution.
  dist_label <- gg_diff_eff$grobs[[strip_id]]$grobs[[2]]$children[[2]]$children[[1]]$label
  colour <- dist$COLOUR[dist$DIST_LABEL == str_wrap(dist_label, width = 100)]
  gg_diff_eff$grobs[[strip_id]]$grobs[[2]]$children[[1]]$gp$fill <- colorspace::lighten(
    colour,
    amount = 0.5
  )
}

# Save the figure: effectiveness on each island.
file_name <- paste(
  "fig_compare_opt_island_balance",
  paste(strat_ids, collapse = "_"),
  paste(gsub("%", "", sens_par_1_filter), collapse = "_"),
  sep = "_"
)
ggsave(
  filename = paste0(file_name, ".svg"), plot = gg_diff_eff, width = 9, height = 11, scale = 1.1,
  path = "figures"
)
ggsave(
  filename = paste0(file_name, ".pdf"), plot = gg_diff_eff, width = 9, height = 11, scale = 1.1,
  path = "figures"
)

# Generate a table of statistics for the main text.
stats_df <- archi_data %>%
  as.data.frame() %>%
  group_by(SENS_PAR_1, TAG_ID, DIST_ID, ISLAND_ID) %>%
  select(SENS_PAR_1, TAG_ID, DIST_ID, ISLAND_ID, DIFF_EFF_Q500, DIFF_EFF_Q025, DIFF_EFF_Q975) %>%
  mutate(across(.cols = contains("DIFF_EFF"), .fns = \(x) round(100 * x, digits = 2))) %>%
  as.data.frame()

# Format for a LaTeX table.
tex_tabular <- stats_df %>%
  arrange(SENS_PAR_1, TAG_ID, DIST_ID, ISLAND_ID) %>%
  mutate(
    LABEL = paste(
      paste0("\\texttt{", TAG_ID, DIST_ID, "}"),
      SENS_PAR_1,
      island_labels[as.character(ISLAND_ID)],
      paste0(DIFF_EFF_Q500, "% [", DIFF_EFF_Q025, "%, ", DIFF_EFF_Q975, "%] \\\\"),
      sep = " & "
    ),
    LABEL = gsub("\\%", "\\\\%", LABEL),
    .keep = "none"
  )
tex_tabular <- paste(tex_tabular$LABEL, collapse = "\n")