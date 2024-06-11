# plot_compare_opt_inf_time.R: plotting to compare the number of infections
# over time for a range of vaccination strategies.
# Pre-requisites: optimal parameters already obtained for each vaccination rate and tagging strategy,
#                 optimal parameters fed through simulation under corresponding scenarios,
#                 simulations have already been summarised using data_summarise_scenario_time.R.

# Clear the workspace.
rm(list = ls())

# Load in the data visualisation and manipulation libraries.
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
  expand(
    nesting(DIST_ID, DIST_LABEL, COLOUR),
    nesting(TAG_ID, TAG_LABEL, LIGHTEN)
  ) %>%
  drop_na(everything()) %>%
  unite(col = "STRAT_ID", TAG_ID, DIST_ID, sep = "", remove = FALSE)

# Calculate the total population of the Comoros archipelago.
total_pop <- fread(paste0(in_folders[1], "/weekly_simulations/default_pars_000.csv"), header = TRUE) %>%
  filter(PAR_NAME == "n_pop") %>%
  .$PAR_VALUE %>%
  sum()

# Get the population size per island.
island_pop <- fread(paste0(in_folders[1], "/weekly_simulations/default_pars_000.csv"), header = TRUE) %>%
  filter(PAR_NAME == "n_pop") %>%
  mutate(ISLAND_ID = factor(PAR_INDEX, levels = c(1, 3, 0, 2))) %>%
  select(ISLAND_ID, PAR_VALUE) %>%
  rename(POP_SIZE = PAR_VALUE)

# Define the settings for the main parameter of interest in common between all optimisations/scenarios.
sens_par_name_1 <- "vac_rate"
trans_sens_par_1 <- \(x) {
  round(x / total_pop, digits = 2)
}
sens_par_title_1 <- "Percentage of livestock vaccinated annually across\nthe Comoros archipelago"
sens_par_label_1 <- scales::percent
sens_par_filter_1 <- c("20%", "25%")

# Data frame to store all summary information.
time_df <- data.frame()

# Initialise a counter for ID numbers.
i <- 0

# For each opt/scenario...
for (in_folder in in_folders) {
  # Get the number of unique summary files.
  file_ids <- gsub("[A-z.]", "", list.files(
    path = paste0(in_folder, "/weekly_simulations"),
    pattern = "rnaught_[0-9]+.csv"
  ))

  # For each file...
  for (file_id in file_ids) {
    # Load in the scenario parameters.
    scen_par_df <- fread(paste0(in_folder, "/weekly_simulations/scenario_pars_", file_id, ".csv"), header = TRUE)

    # Calculate the number of indices per scenario parameter, combine the parameter
    # name and index columns, and then expand the data frame to a wide format.
    sens_par <- scen_par_df %>%
      mutate(PAR_NAME = ifelse(PAR_NAME == sens_par_name_1, "SENS_PAR_1", PAR_NAME)) %>%
      mutate(PAR_VALUE = ifelse(PAR_NAME == "SENS_PAR_1", trans_sens_par_1(PAR_VALUE), PAR_VALUE)) %>%
      filter(PAR_NAME == "SENS_PAR_1") %>%
      pivot_wider(id_cols = SAMPLE_ID, names_from = PAR_NAME, values_from = PAR_VALUE) %>%
      slice(1) %>%
      pull(SENS_PAR_1)

    # Only proceed if there are valid values for the sensitivity parameter.
    if (nrow(scen_par_df) > 0) {

      # Load in the summary data.
      one_time_df <- fread(paste0(in_folder, "/weekly_simulations/summary_time_", file_id, ".csv"), header = TRUE)

      # Add in the sensitivity parameter.
      one_time_df <- one_time_df %>%
        mutate(SENS_PAR_1 = sens_par)

      # Note the input scenario / optimisation name.
      one_time_df$STRAT_ID <- names(in_folders)[which(in_folder == in_folders)]
      one_time_df$GROUP_ID <- i

      # Move onto the next id number.
      i <- i + 1

      # Append to the overall data frame.
      time_df <- rbind(time_df, one_time_df)
    }
  }
}

# Do the same for the optimisation parameters.
time_df <- time_df %>%
  mutate(
    ISLAND_ID = factor(ISLAND_ID, levels = c(1, 3, 0, 2)),
    FACET_ID = factor(ISLAND_ID, levels = c(1, 0, 3, 2)),
    SENS_PAR_1 = factor(sens_par_label_1(SENS_PAR_1))
  ) %>%
  left_join(
    all_strats,
    by = "STRAT_ID"
  )

# Convert the levels of the sensitivity parameter.
time_df$SENS_PAR_1 <- factor(time_df$SENS_PAR_1,
  levels = paste0(sort(as.numeric(gsub("%", "", levels(time_df$SENS_PAR_1)))), "%")
)

# Map island indices to names.
island_labels <- c(
  "0" = "Anjouan",
  "1" = "Grande Comore",
  "2" = "Mayotte",
  "3" = "Mohéli"
)

# Define the way of labelling the sensitivity parameter panels.
sens_par_labels <- paste0(sens_par_title_1, ": ", c("0%", as.character(unique(time_df$SENS_PAR_1))))
names(sens_par_labels) <- c("0%", as.character(unique(time_df$SENS_PAR_1)))

# The number of infections over time for different strategies.
gg_inf <- time_df %>%
  filter(SENS_PAR_1 %in% sens_par_filter_1) %>%
  ggplot() +
  geom_line(
    data = ~ filter(., INT == "Median") %>%
      left_join(island_pop, by = "ISLAND_ID") %>%
      mutate(PROP_INF = 100 * INF / POP_SIZE),
    mapping = aes(
      x = TIME / 48,
      y = PROP_INF,
      colour = ISLAND_ID,
      alpha = INT
    )
  ) +
  geom_ribbon(
    data = ~ filter(., INT == "CrI95") %>%
      left_join(island_pop, by = "ISLAND_ID") %>%
      mutate(PROP_INF = 100 * INF / POP_SIZE) %>%
      pivot_wider(
        id_cols = c(DATE:ISLAND_ID, INT, SENS_PAR_1:LIGHTEN),
        names_from = SIDE, names_sep = "_",
        values_from = PROP_INF
      ),
    mapping = aes(
      x = TIME / 48,
      ymin = Lower,
      ymax = Upper,
      fill = ISLAND_ID,
      alpha = INT
    )
  ) +
  geom_vline(
    xintercept = 0,
    alpha = 0.5,
    linetype = "dashed"
  ) +
  scale_x_continuous(
    name = "Time since the first vaccine was administered (years)"
  ) +
  scale_y_continuous(
    paste0(
      "Number of livestock infected with Rift Valley fever virus per 100 animals"
    )
  ) +
  scale_colour_brewer(
    name = "Island",
    palette = "Set1",
    breaks = levels(time_df$FACET_ID),
    limits = levels(time_df$ISLAND_ID),
    labels = island_labels,
    guide = "none"
  ) +
  scale_fill_brewer(
    name = "Island",
    palette = "Set1",
    breaks = levels(time_df$FACET_ID),
    limits = levels(time_df$ISLAND_ID),
    labels = island_labels
  ) +
  scale_alpha_ordinal(
    name = "Summary statistic",
    limits = \(x) c("Median", setdiff(x, "Median")),
    labels = function(x) {
      y <- ifelse(
        x == "Median",
        "Median",
        paste0(gsub("CrI", "", x), "% credible interval")
      )
      return(y)
    },
    range = c(1, 0.2)
  ) +
  facet_grid(
    DIST_LABEL + TAG_LABEL ~ SENS_PAR_1,
    labeller = labeller(
      SENS_PAR_1 = as_labeller(sens_par_labels, default = label_wrap_gen(width = 30)),
      DIST_LABEL = label_wrap_gen(width = 30)
    )
  ) +
  gg_theme +
  theme(
    panel.border = element_rect(fill = NA, colour = "black")
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      order = 2,
      byrow = TRUE,
      nrow = 2,
      override.aes = list(
        colour = "black"
      )
    ),
    alpha = guide_legend(
      title.position = "top",
      order = 1,
      byrow = TRUE,
      nrow = 2,
      override.aes = list(
        fill = c(NA, "black"),
        colour = "black"
      )
    )
  )

# Build the figure and edit the legend to hide legend lines / figures for summary statistics.
gg_inf <- ggplot_gtable(ggplot_build(gg_inf))
gg_inf$grobs[[32]]$grobs[[1]]$grobs[[5]] <- zeroGrob()
gg_inf$grobs[[32]]$grobs[[1]]$grobs[[7]] <- zeroGrob()

# Also edit the colour of the strips.
strip_ids <- grep("strip-r", gg_inf$layout$name)
for (strip_id in strip_ids){
  # In each strip, the first strip text is associated with
  # the tagged status of livestock.
  tag_label <- gg_inf$grobs[[strip_id]]$grobs[[1]]$children[[2]]$children[[1]]$label
  lighten <- tag$LIGHTEN[tag$TAG_LABEL == tag_label]
  gg_inf$grobs[[strip_id]]$grobs[[1]]$children[[1]]$gp$fill <- colorspace::lighten(
    "grey85", amount = lighten / 2
  )

  # Now change the colour for the vaccine distribution.
  dist_label <- gg_inf$grobs[[strip_id]]$grobs[[2]]$children[[2]]$children[[1]]$label
  colour <- dist$COLOUR[dist$DIST_LABEL == str_wrap(dist_label, width = 100)]
  gg_inf$grobs[[strip_id]]$grobs[[2]]$children[[1]]$gp$fill <- colorspace::lighten(
    colour,
    amount = 0.5
  )
}

# Save the figure.
# Optimal parameters on each island.
file_name <- paste(
  "fig_compare_opt_inf_time",
  paste(strat_ids, collapse = "_"),
  paste(gsub("%", "", sens_par_filter_1), collapse = "_"),
  sep = "_"
)
ggsave(
  filename = paste0(file_name, ".svg"), plot = gg_inf, width = 9, height = 11, scale = 1.1,
  path = "figures"
)
ggsave(
  filename = paste0(file_name, ".pdf"), plot = gg_inf, width = 9, height = 11, scale = 1.1,
  path = "figures"
)

# The number of infections averted over time for different strategies.
gg_inf_avert <- time_df %>%
  filter(SENS_PAR_1 %in% sens_par_filter_1) %>%
  ggplot() +
  geom_line(
    data = ~ filter(., INT == "Median") %>%
      left_join(island_pop, by = "ISLAND_ID") %>%
      mutate(INF_AVERTED = 100 * INF_AVERTED / POP_SIZE),
    mapping = aes(
      x = TIME / 48,
      y = INF_AVERTED,
      colour = ISLAND_ID,
      alpha = INT
    )
  ) +
  geom_ribbon(
    data = ~ filter(., INT == "CrI95") %>%
      left_join(island_pop, by = "ISLAND_ID") %>%
      mutate(INF_AVERTED = 100 * INF_AVERTED / POP_SIZE) %>%
      pivot_wider(
        id_cols = c(DATE:ISLAND_ID, INT, SENS_PAR_1:LIGHTEN),
        names_from = SIDE, names_sep = "_",
        values_from = INF_AVERTED
      ),
    mapping = aes(
      x = TIME / 48,
      ymin = Lower,
      ymax = Upper,
      fill = ISLAND_ID,
      alpha = INT
    )
  ) +
  geom_vline(
    xintercept = 0,
    alpha = 0.5,
    linetype = "dashed"
  ) +
  scale_x_continuous(
    name = "Time since the first vaccine was administered (years)"
  ) +
  scale_y_continuous(
    paste0(
      "Number of Rift Valley fever virus infections averted per 100 animals"
    ),
    limits = \(x) c(-max(abs(x)), max(abs(x)))
  ) +
  scale_colour_brewer(
    name = "Island",
    palette = "Set1",
    breaks = levels(time_df$FACET_ID),
    limits = levels(time_df$ISLAND_ID),
    labels = island_labels,
    guide = "none"
  ) +
  scale_fill_brewer(
    name = "Island",
    palette = "Set1",
    breaks = levels(time_df$FACET_ID),
    limits = levels(time_df$ISLAND_ID),
    labels = island_labels
  ) +
  scale_alpha_ordinal(
    name = "Summary statistic",
    limits = \(x) c("Median", setdiff(x, "Median")),
    labels = function(x) {
      y <- ifelse(
        x == "Median",
        "Median",
        paste0(gsub("CrI", "", x), "% credible interval")
      )
      return(y)
    },
    range = c(1, 0.2)
  ) +
  facet_grid(
    DIST_LABEL + TAG_LABEL ~ SENS_PAR_1,
    labeller = labeller(
      SENS_PAR_1 = as_labeller(sens_par_labels, default = label_wrap_gen(width = 30)),
      DIST_LABEL = label_wrap_gen(width = 30)
    )
  ) +
  gg_theme +
  theme(
    panel.border = element_rect(fill = NA, colour = "black")
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      order = 2,
      byrow = TRUE,
      nrow = 2,
      override.aes = list(
        colour = "black"
      )
    ),
    alpha = guide_legend(
      title.position = "top",
      order = 1,
      byrow = TRUE,
      nrow = 2,
      override.aes = list(
        fill = c(NA, "black"),
        colour = "black"
      )
    )
  )

# Build the figure and edit the legend to hide legend lines / figures for summary statistics.
gg_inf_avert <- ggplot_gtable(ggplot_build(gg_inf_avert))
gg_inf_avert$grobs[[32]]$grobs[[1]]$grobs[[5]] <- zeroGrob()
gg_inf_avert$grobs[[32]]$grobs[[1]]$grobs[[7]] <- zeroGrob()

# Also edit the colour of the strips.
strip_ids <- grep("strip-r", gg_inf_avert$layout$name)
for (strip_id in strip_ids){
  # In each strip, the first strip text is associated with
  # the tagged status of livestock.
  tag_label <- gg_inf_avert$grobs[[strip_id]]$grobs[[1]]$children[[2]]$children[[1]]$label
  lighten <- tag$LIGHTEN[tag$TAG_LABEL == tag_label]
  gg_inf_avert$grobs[[strip_id]]$grobs[[1]]$children[[1]]$gp$fill <- colorspace::lighten(
    "grey85", amount = lighten / 2
  )

  # Now change the colour for the vaccine distribution.
  dist_label <- gg_inf_avert$grobs[[strip_id]]$grobs[[2]]$children[[2]]$children[[1]]$label
  colour <- dist$COLOUR[dist$DIST_LABEL == str_wrap(dist_label, width = 100)]
  gg_inf_avert$grobs[[strip_id]]$grobs[[2]]$children[[1]]$gp$fill <- colorspace::lighten(
    colour,
    amount = 0.5
  )
}

# Save the figure.
# Optimal parameters on each island.
file_name <- paste(
  "fig_compare_opt_infavert_time",
  paste(strat_ids, collapse = "_"),
  paste(gsub("%", "", sens_par_filter_1), collapse = "_"),
  sep = "_"
)
ggsave(
  filename = paste0(file_name, ".svg"), plot = gg_inf_avert, width = 9, height = 11, scale = 1.1,
  path = "figures"
)
ggsave(
  filename = paste0(file_name, ".pdf"), plot = gg_inf_avert, width = 9, height = 11, scale = 1.1,
  path = "figures"
)

# Plot the number of infections over time with no vaccination.
gg_inf_no_vac <- time_df %>%
  filter(SENS_PAR_1 %in% "5%", STRAT_ID == strat_ids[1]) %>%
  mutate(SENS_PAR_1 = "0%") %>%
  ggplot() +
  geom_line(
    data = ~ filter(., INT == "Median") %>%
      left_join(island_pop, by = "ISLAND_ID") %>%
      mutate(PROP_INF = 100 * INF_NO_VAC / POP_SIZE),
    mapping = aes(
      x = TIME / 48,
      y = PROP_INF,
      colour = ISLAND_ID,
      alpha = INT
    )
  ) +
  geom_ribbon(
    data = ~ filter(., INT == "CrI95") %>%
      left_join(island_pop, by = "ISLAND_ID") %>%
      mutate(PROP_INF = 100 * INF_NO_VAC / POP_SIZE) %>%
      pivot_wider(
        id_cols = c(DATE:ISLAND_ID, INT, SENS_PAR_1:LIGHTEN),
        names_from = SIDE, names_sep = "_",
        values_from = PROP_INF
      ),
    mapping = aes(
      x = TIME / 48,
      ymin = Lower,
      ymax = Upper,
      fill = ISLAND_ID,
      alpha = INT
    )
  ) +
  geom_vline(
    xintercept = 0,
    alpha = 0.5,
    linetype = "dashed"
  ) +
  scale_x_continuous(
    name = "Time since the first vaccine would have been administered (years)"
  ) +
  scale_y_continuous(
    paste0(
      "Number of livestock infected with Rift Valley\nfever virus per 100 animals"
    )
  ) +
  scale_colour_brewer(
    name = "Island",
    palette = "Set1",
    breaks = levels(time_df$FACET_ID),
    limits = levels(time_df$ISLAND_ID),
    labels = island_labels,
    guide = "none"
  ) +
  scale_fill_brewer(
    name = "Island",
    palette = "Set1",
    breaks = levels(time_df$FACET_ID),
    limits = levels(time_df$ISLAND_ID),
    labels = island_labels
  ) +
  scale_alpha_ordinal(
    name = "Summary statistic",
    limits = \(x) c("Median", setdiff(x, "Median")),
    labels = function(x) {
      y <- ifelse(
        x == "Median",
        "Median",
        paste0(gsub("CrI", "", x), "% credible interval")
      )
      return(y)
    },
    range = c(1, 0.2)
  ) +
  facet_grid(
     ~ SENS_PAR_1,
    labeller = labeller(
      SENS_PAR_1 = as_labeller(sens_par_labels, default = label_wrap_gen(width = 60))
    )
  ) +
  gg_theme +
  theme(
    panel.border = element_rect(fill = NA, colour = "black")
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      order = 2,
      byrow = TRUE,
      nrow = 2,
      override.aes = list(
        colour = "black"
      )
    ),
    alpha = guide_legend(
      title.position = "top",
      order = 1,
      byrow = TRUE,
      nrow = 2,
      override.aes = list(
        fill = c(NA, "black"),
        colour = "black"
      )
    )
  )

# Build the figure and edit the legend to hide legend lines / figures for summary statistics.
gg_inf_no_vac <- ggplot_gtable(ggplot_build(gg_inf_no_vac))
gg_inf_no_vac$grobs[[12]]$grobs[[1]]$grobs[[5]] <- zeroGrob()
gg_inf_no_vac$grobs[[12]]$grobs[[1]]$grobs[[7]] <- zeroGrob()

# Save the figure.
# Optimal parameters on each island.
file_name <- paste(
  "fig_compare_opt_inf_no_vac_time",
  paste(strat_ids, collapse = "_"),
  sep = "_"
)
ggsave(
  filename = paste0(file_name, ".svg"), plot = gg_inf_no_vac, width = 9, height = 6, scale = 1.1,
  path = "figures"
)
ggsave(
  filename = paste0(file_name, ".pdf"), plot = gg_inf_no_vac, width = 9, height = 6, scale = 1.1,
  path = "figures"
)
