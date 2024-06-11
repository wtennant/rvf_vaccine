# plot_compare_opt_freq_time.R: Visualise how optimised parameter sets depend on
# different vaccination rates, tagging strategies, frequency and timing of vaccination.
# Pre-requisites: optimal parameters already obtained for each vaccination rate and tagging strategy,
#                 optimal parameters fed through simulation under corresponding scenarios,
#                 simulations have already been summarised using data_summarise_scenario.R.

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
strat_ids <- c("U1", "T1")

# Define the folder locations of all optimisations/scenarios as a named vector.
in_folders <- c(
  paste0(
    rep("../cpp/out/", length(strat_ids)),
    strat_ids,
    "_freq_time"
  )
)
names(in_folders) <- strat_ids

# Define the identifiers and labels of strategy options.
tag <- data.frame(
  TAG_ID = factor(c("U", "T"), levels = c("U", "T")),
  TAG_LABEL = c("Untagged", "Tagged"),
  LIGHTEN = c(0.4, -0.4)
) %>%
  mutate(
    TAG_LABEL = factor(TAG_LABEL, levels = TAG_LABEL)
  )
dist <- data.frame(
  DIST_ID = factor(c("P", "1", "2"), levels = c("P", "1", "2")),
  DIST_LABEL = c(
    "Proportional to island population size",
    "Optimised over infections averted across the archipelago",
    "Optimised over infections averted on the worst-peforming island"
  ),
  COLOUR = c("#777777", palette_okabe_ito(order = c(1, 2)))
) %>%
  mutate(
    DIST_LABEL = factor(DIST_LABEL, levels = DIST_LABEL)
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
total_pop <- fread(paste0(in_folders[1], "/simulations/default_pars_000.csv"), header = TRUE) %>%
  filter(PAR_NAME == "n_pop") %>%
  .$PAR_VALUE %>%
  sum()

# Define the settings for the parameters of interest.
# Parameter definition of the x-axis.
sens_par_name_1 <- "vac_t_year_length"
trans_sens_par_1 <- \(x) {
  x
}
sens_par_title_1 <- "Vaccination timing each epidemiological year"
sens_par_label_1 <- \(x) {
  ifelse(x == 48, "Throughout", paste("First month"))
}

# Parameter definition for the y-axis.
sens_par_name_2 <- "vac_t_freq"
trans_sens_par_2 <- \(x) {
  1 / x
}
sens_par_title_2 <- "Frequency of vaccination campaigns (per year)"
sens_par_label_2 <- \(x) {
  as.character(MASS::fractions(x))
}

# Define the settings for the columns of the panel arrangement.
sens_par_name_3 <- "vac_rate"
trans_sens_par_3 <- \(x) {
  round(x / total_pop, digits = 2)
}
sens_par_title_3 <- "Percentage vaccinated annually"
sens_par_label_3 <- scales::percent

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
      mutate(PAR_NAME = ifelse(PAR_NAME == sens_par_name_2, "SENS_PAR_2", PAR_NAME)) %>%
      mutate(PAR_VALUE = ifelse(PAR_NAME == "SENS_PAR_2", trans_sens_par_2(PAR_VALUE), PAR_VALUE)) %>%
      mutate(PAR_NAME = ifelse(PAR_NAME == sens_par_name_3, "SENS_PAR_3", PAR_NAME)) %>%
      mutate(PAR_VALUE = ifelse(PAR_NAME == "SENS_PAR_3", trans_sens_par_3(PAR_VALUE), PAR_VALUE)) %>%
      filter(PAR_NAME %in% c("SENS_PAR_1", "SENS_PAR_2", "SENS_PAR_3")) %>%
      pivot_wider(id_cols = SAMPLE_ID, names_from = PAR_NAME, values_from = PAR_VALUE)

    # Only proceed if there are valid values for the sensitivity parameter.
    if (nrow(scen_par_df) > 0) {

      # Load in the summary data.
      one_summary_df <- fread(paste0(in_folder, "/simulations/summary_", file_id, ".csv"), header = TRUE)

      # Join into summary data frame.
      one_summary_df <- right_join(one_summary_df, scen_par_df, by = "SAMPLE_ID")

      # Summarise the summary data(!). Extract the median, 97.5% and 2.5% quantiles.
      one_summary_df <- one_summary_df %>%
        group_by(SENS_PAR_1, SENS_PAR_2, SENS_PAR_3, ISLAND_ID) %>%
        summarise(
          EFF_Q025 = quantile(EFFECTIVENESS, probs = 0.025),
          EFF_Q500 = quantile(EFFECTIVENESS, probs = 0.5),
          EFF_Q975 = quantile(EFFECTIVENESS, probs = 0.975),
          EFF_SD = sd(EFFECTIVENESS),
          .groups = "drop"
        ) %>%
        mutate(EFF_LABEL = paste0(
          100 * round(EFF_Q500, digits = 2), "%\n[", 100 * round(EFF_Q025, digits = 2),
          "%, ", 100 * round(EFF_Q975, digits = 2), "%]"
        ))

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
    FACET_ID = factor(ISLAND_ID, levels = c(-1, 1, 3, 0, 2)),
    SENS_PAR_1 = factor(sens_par_label_1(SENS_PAR_1),
      levels = sens_par_label_1(sort(unique(summary_df$SENS_PAR_1)))
    ),
    SENS_PAR_2 = factor(sens_par_label_2(SENS_PAR_2),
      levels = sens_par_label_2(sort(unique(summary_df$SENS_PAR_2)))
    ),
    SENS_PAR_3 = factor(sens_par_label_3(SENS_PAR_3),
      levels = sens_par_label_3(rev(sort(unique(summary_df$SENS_PAR_3))))
    )
  ) %>%
  left_join(
    all_strats,
    by = "STRAT_ID"
  )

# Map island indices to names.
island_labels <- c(
  "0" = "Anjouan",
  "1" = "Grande Comore",
  "2" = "Mayotte",
  "3" = "Mohéli",
  "-1" = "Comoros archipelago"
)

# Construct a look-up table for the facets.
facet_labels <- paste0(sens_par_title_3, ": ", unique(summary_df$SENS_PAR_3))
names(facet_labels) <- unique(summary_df$SENS_PAR_3)

# And for the entire archipelago...
gg_cost <- summary_df %>%
  filter(ISLAND_ID == -1) %>%
  ggplot() +
  geom_tile(
    mapping = aes(
      x = SENS_PAR_1,
      y = SENS_PAR_2,
      fill = stage(
        start = EFF_Q500,
        after_scale = after_scale(colorspace::lighten(fill, amount = 0.1, space = "HLS"))
      )
    ),
    colour = "black"
  ) +
  geom_text(mapping = aes(x = SENS_PAR_1, y = SENS_PAR_2, label = EFF_LABEL), colour = "black") +
  scale_x_discrete(sens_par_title_1, expand = c(0, 0)) +
  scale_y_discrete(sens_par_title_2, expand = c(0, 0)) +
  scale_fill_distiller(
    name = "Median percentage\nof infections averted\nacross the archipelago",
    palette = "YlOrRd",
    label = scales::percent,
    limits = c(0, 1),
    n.breaks = 6,
    direction = 1
  ) +
  facet_grid(
    SENS_PAR_3 ~ DIST_LABEL + TAG_LABEL,
    labeller = labeller(
      SENS_PAR_3 = as_labeller(facet_labels, default = label_wrap_gen(width = 30)),
      DIST_LABEL = label_wrap_gen(width = 30)
    )
  ) +
  gg_theme +
  theme(
    panel.border = element_rect(fill = NA, colour = "black"),
    legend.position = "right",
    plot.background = element_blank(),
    legend.background = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
    legend.spacing.y = unit(5, "mm")
  ) +
  guides(
    fill = guide_colourbar(
      frame.colour = "black",
      ticks.colour = "black",
      barheight = unit(100, "mm")
    )
  )

# Also edit the colour of the strips.
gg_cost <- ggplot_gtable(ggplot_build(gg_cost))
strip_ids <- grep("strip-t", gg_cost$layout$name)
for (strip_id in strip_ids){
  # In each strip, the first strip text is associated with
  # the tagged status of livestock.
  tag_label <- gg_cost$grobs[[strip_id]]$grobs[[2]]$children[[2]]$children[[1]]$label
  lighten <- tag$LIGHTEN[tag$TAG_LABEL == tag_label]
  gg_cost$grobs[[strip_id]]$grobs[[2]]$children[[1]]$gp$fill <- colorspace::lighten(
    "grey85", amount = lighten / 2
  )

  # Now change the colour for the vaccine distribution.
  dist_label <- gg_cost$grobs[[strip_id]]$grobs[[1]]$children[[2]]$children[[1]]$label
  colour <- dist$COLOUR[dist$DIST_LABEL == str_wrap(dist_label, width = 100)]
  gg_cost$grobs[[strip_id]]$grobs[[1]]$children[[1]]$gp$fill <- colorspace::lighten(
    colour,
    amount = 0.5
  )
}

# Save the figure.
file_name <- paste("fig_compare_opt_freq_timing", paste(strat_ids, collapse = "_"), sep = "_")
ggsave(
  filename = paste0(file_name, ".pdf"), plot = gg_cost, width = 8, height = 9, dpi = 600, scale = 1.1,
  path = "figures"
)
ggsave(
  filename = paste0(file_name, ".svg"), plot = gg_cost, width = 8, height = 9, dpi = 600, scale = 1.1,
  path = "figures"
)
