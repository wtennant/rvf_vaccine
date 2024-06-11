# plot_compare_opt_dist: plotting to compare the optimal distribution
# of vaccines between islands under a variety of different vaccine strategies.
# Pre-requisites: optimal parameters already obtained for each vaccination rate and tagging strategy,
#                 optimal parameters prepared and summarised with data_format_optimal.R.

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
strat_ids <- c("U1", "T1", "U2", "T2")

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
) %>%
  mutate(
    TAG_LABEL = factor(TAG_LABEL, levels = TAG_LABEL)
  )
dist <- data.frame(
  DIST_ID = factor(c("P", "1", "2"), levels = c("P", "1", "2")),
  DIST_LABEL = c(
    "Proportional\nto island\npopulation\nsize",
    "Optimised over infections averted\nacross the archipelago",
    "Optimised over infections averted\non the worst-peforming island"
  ),
  COLOUR = c("#777777", palette_okabe_ito(order = c(1, 2)))
) %>%
  mutate(
    DIST_LABEL = factor(DIST_LABEL, levels = DIST_LABEL)
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

# Calculate the distribution of vaccines according to population size.
prop_dist <- fread(paste0(in_folders[1], "/simulations/default_pars_000.csv"), header = TRUE) %>%
  filter(PAR_NAME == "n_pop") %>%
  mutate(PAR_VALUE = PAR_VALUE / total_pop) %>%
  rename(ISLAND_ID = PAR_INDEX, PROP_DIST = PAR_VALUE) %>%
  mutate(
    ISLAND_ID = factor(ISLAND_ID, levels = c(1, 3, 0, 2, -1)),
    FACET_ID = factor(ISLAND_ID, levels = c(1, 0, 3, 2, -1))
  ) %>%
  mutate("DIST_ID" = "P") %>%
  left_join(
    all_strats,
    by = "DIST_ID",
    relationship = "many-to-many"
  )

# Define the settings for the main parameter of interest in common between all optimisations/scenarios.
sens_par_name_1 <- "vac_rate"
trans_sens_par_1 <- \(x) {
  round(x / total_pop, digits = 2)
}
sens_par_title_1 <- "Percentage of livestock vaccinated annually across\nthe Comoros archipelago"
sens_par_label_1 <- scales::percent

# Data frame to store all summary information.
summary_df <- opt_df <- data.frame()

# Initialise a counter for ID numbers.
i <- j <- 0

# For each opt/scenario...
for (in_folder in in_folders) {
  # Get the names of the optimised parameter files.
  opt_files <- list.files(path = in_folder, pattern = "optimal_", full.names = TRUE)

  # For each optimisation file...
  for (opt_file in opt_files) {
    # Load in the set of optimal parameters.
    one_opt_df <- fread(opt_file, header = TRUE)

    # Re-format the sensitivity parameter, and extract the parameter being optimised
    # (i.e. vaccine distribution).
    one_opt_df <- one_opt_df %>%
      filter(PAR_NAME == sens_par_name_1) %>%
      mutate(PAR_NAME = ifelse(PAR_NAME == sens_par_name_1, "SENS_PAR_1", PAR_NAME)) %>%
      mutate(PAR_VALUE = ifelse(PAR_NAME == "SENS_PAR_1", trans_sens_par_1(PAR_VALUE), PAR_VALUE)) %>%
      pivot_wider(id_cols = SAMPLE_ID, names_from = PAR_NAME, values_from = PAR_VALUE) %>%
      left_join(one_opt_df, by = "SAMPLE_ID", multiple = "all") %>%
      filter(PAR_NAME == "vac_dist") %>%
      select(-PAR_NAME) %>%
      rename(ISLAND_ID = PAR_INDEX, OPT_DIST = PAR_VALUE)

    # Note the input scenario / optimisation name.
    one_opt_df$STRAT_ID <- names(in_folders)[which(in_folder == in_folders)]
    one_opt_df$GROUP_ID <- j
    j <- j + 1

    # Combine into the overall data frame of optimised parameters.
    opt_df <- rbind(opt_df, one_opt_df)
  }
}

# Do the same for the optimisation parameters.
opt_df <- opt_df %>%
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
opt_df$SENS_PAR_1 <- factor(opt_df$SENS_PAR_1,
  levels = paste0(sort(as.numeric(gsub("%", "", levels(opt_df$SENS_PAR_1)))), "%")
)

# Map island indices to names.
island_labels <- c(
  "0" = "Anjouan",
  "1" = "Grande Comore",
  "2" = "Mayotte",
  "3" = "Mohéli",
  "-1" = "Comoros archipelago"
)

# Add the proportional to population size allocation
opt_df <- opt_df %>%
  expand(SAMPLE_ID, ISLAND_ID) %>%
  mutate(SENS_PAR_1 = factor("Any%", levels = c("Any%", levels(opt_df$SENS_PAR_1)))) %>%
  right_join(prop_dist, by = c("ISLAND_ID"), relationship = "many-to-many") %>%
  select(-PAR_NAME) %>%
  rename(OPT_DIST = PROP_DIST) %>%
  mutate(
    GROUP_ID = -as.numeric(SENS_PAR_1)
  ) %>%
  bind_rows(opt_df)

# First visualise this as a 'cutting' problem.
# Define the X-axis as the sample number crossed with the
# vaccination rate.
prob_df <- opt_df %>%
  select(-FACET_ID) %>%
  pivot_wider(names_from = ISLAND_ID, values_from = OPT_DIST) %>%
  mutate(SENS_PAR_ID = as.numeric(SENS_PAR_1)) %>%
  group_by(STRAT_ID, SENS_PAR_1) %>%
  arrange(SAMPLE_ID, .by_group = TRUE) %>%
  mutate(X = seq_len(n()) + SENS_PAR_ID * n() * 1.25) %>%
  pivot_longer(cols = `1`:`2`, names_to = "ISLAND_ID", values_to = "OPT_DIST") %>%
  mutate(ISLAND_ID = factor(ISLAND_ID, levels = c(1, 3, 0, 2)))

# Define the labels for each vaccination rate.
x_axis_info <- prob_df %>%
  group_by(SENS_PAR_1) %>%
  reframe(BREAKS = mean(X)) %>%
  rename(LABELS = SENS_PAR_1)

# Comparing the distributions of vaccines across the islands.
gg_opt <- prob_df %>%
  ggplot() +
  geom_col(
    mapping = aes(
      x = X,
      y = OPT_DIST,
      fill = ISLAND_ID
    ),
    width = 1,
    colour = NA,
    position = "stack"
  ) +
  scale_x_continuous(
    name = sens_par_title_1,
    breaks = x_axis_info$BREAKS,
    labels = x_axis_info$LABELS,
    expand = c(0, 250)
  ) +
  scale_y_continuous(
    name = "Percentage share of vaccines between islands in the Comoros archipelago",
    label = scales::percent,
    expand = c(0, 0)
  ) +
  scale_fill_brewer(
    "Island",
    palette = "Set1",
    labels = island_labels
  ) +
  facet_grid(
    TAG_LABEL ~ DIST_LABEL,
    scales = "free_x",
    space = "free_x"
  ) +
  gg_theme +
  theme(
    panel.border = element_rect(fill = NA, colour = "black"),
    panel.spacing = unit(1, "lines"),
    strip.clip = "off"
  ) +
  guides(
    fill = guide_legend(
      title.position = "top",
      override.aes = list(
        colour = "black"
      )
    )
  )

# Convert to a gtable and scale down the strip text for the proportional category.
# This conversion is hard coded for a set input of strategies to plot.
if (all(strat_ids == c("U1", "T1", "U2", "T2"))) {
  gg_opt <- ggplot_gtable(ggplot_build(gg_opt))
  gg_opt$grobs[[18]]$grobs[[1]]$children[[2]]$children[[1]]$gp$fontsize <- 9

  # Also re-colour the facet strips to correspond with the strategy colouring.
  # The first strip text is associated with the tagged status of livestock.
  strip_ids <- grep("strip-r", gg_opt$layout$name)
  for (strip_id in strip_ids){
    tag_label <- gg_opt$grobs[[strip_id]]$grobs[[1]]$children[[2]]$children[[1]]$label
    lighten <- tag$LIGHTEN[tag$TAG_LABEL == tag_label]
    gg_opt$grobs[[strip_id]]$grobs[[1]]$children[[1]]$gp$fill <- colorspace::lighten(
      "grey85", amount = lighten / 2
    )
  }

  # Now change the colour for the vaccine distribution.
  strip_ids <- grep("strip-t", gg_opt$layout$name)
  for (strip_id in strip_ids){
    dist_label <- gg_opt$grobs[[strip_id]]$grobs[[1]]$children[[2]]$children[[1]]$label
    colour <- dist$COLOUR[str_wrap(dist$DIST_LABEL, width = 100) == str_wrap(dist_label, width = 100)]
    gg_opt$grobs[[strip_id]]$grobs[[1]]$children[[1]]$gp$fill <- colorspace::lighten(
      colour,
      amount = 0.5
    )
  }
}

# Save the figure.
# Optimal parameters on each island.
file_name <- paste("fig_compare_opt_dist", paste(strat_ids, collapse = "_"), sep = "_")
ggsave(
  filename = paste0(file_name, ".svg"), plot = gg_opt, width = 9, height = 9, scale = 1.1,
  path = "figures"
)
ggsave(
  filename = paste0(file_name, ".pdf"), plot = gg_opt, width = 9, height = 9, scale = 1.1,
  path = "figures"
)

# As an alternative, visualise the percentage of livestock vaccinated on each island per year.
prop_island_vac_df <- fread(paste0(in_folders[1], "/simulations/default_pars_000.csv"), header = TRUE) %>%
  filter(PAR_NAME == "n_pop") %>%
  mutate(TOTAL_POP = sum(PAR_VALUE)) %>%
  mutate(PAR_INDEX = factor(PAR_INDEX, levels = c(1, 3, 0, 2))) %>%
  right_join(opt_df, by = c("PAR_INDEX" = "ISLAND_ID")) %>%
  filter(SENS_PAR_1 != "Any%") %>%
  mutate(SENS_PAR_1 = droplevels(SENS_PAR_1, "Any%")) %>%
  mutate(VAC_RATE = as.numeric(gsub("%", "", SENS_PAR_1)) / 100) %>%
  mutate(PROP_ANIMAL_VAC = VAC_RATE * OPT_DIST * TOTAL_POP / PAR_VALUE)

# Visualise the percentage of livestock vaccinated on each island per year.
gg_island_vac <- prop_island_vac_df %>%
  ggplot() +
  geom_violin(
    mapping = aes(
      x = SENS_PAR_1,
      y = PROP_ANIMAL_VAC,
      fill = DIST_ID,
      group = GROUP_ID
    ),
    colour = NA,
    width = 0.9,
    scale = "width",
    position = "dodge"
  ) +
  stat_summary(
    geom = "pointrange",
    mapping = aes(
      x = SENS_PAR_1,
      y = PROP_ANIMAL_VAC,
      group = GROUP_ID,
      shape = TAG_ID,
      colour = stage(
        start = DIST_ID,
        after_scale = after_scale(darken(colour, amount = 0.5))
      )
    ),
    fill = "white",
    fun.min = \(x) quantile(x, probs = 0.25),
    fun = median,
    fun.max = \(x) quantile(x, probs = 0.75),
    position = position_dodge(width = 0.9),
    size = 0.5,
    stroke = 0.5
  ) +
  geom_vline(
    data = ~ data.frame(X = seq_len(nlevels(.$SENS_PAR_1) - 1)),
    mapping = aes(xintercept = X + 0.5),
    colour = "black",
    alpha = 0.25,
    linetype = "11"
  ) +
  geom_segment(
    data = ~ data.frame(
      X = seq_len(nlevels(.$SENS_PAR_1)),
      Y = as.numeric(gsub("%", "", levels(.$SENS_PAR_1))) / 100,
      DIST_ID = "P"
    ),
    mapping = aes(
      x = X - 0.5,
      xend = X + 0.5,
      y = Y,
      yend = Y,
      colour = DIST_ID
    ),
    linewidth = 1,
    linetype = "22",
  ) +
  scale_x_discrete(sens_par_title_1) +
  scale_y_continuous("Percentage of livestock vaccinated annually on the island",
    label = scales::percent
  ) +
  scale_fill_manual(
    name = "Vaccine distribution",
    limits = dist$DIST_ID,
    breaks = dist$DIST_ID,
    labels = str_wrap(dist$DIST_LABEL),
    values = dist$COLOUR
  ) +
  scale_colour_manual(
    name = "Vaccine distribution",
    limits = dist$DIST_ID,
    breaks = dist$DIST_ID,
    labels = dist$DIST_LABEL,
    values = dist$COLOUR,
    guide = "none"
  ) +
  scale_shape_manual(
    name = "Livestock",
    values = c("circle filled", "circle"),
    limits = tag$TAG_ID,
    breaks = tag$TAG_ID,
    labels = tag$TAG_LABEL
  ) +
  facet_wrap(~FACET_ID, labeller = labeller(FACET_ID = island_labels), ncol = 2) +
  gg_theme +
  theme(
    panel.border = element_rect(fill = NA, colour = "black"),
    legend.direction = "vertical",
    legend.box.just = "top"
  ) +
  guides(
    shape = guide_legend(
      title.position = "top",
      byrow = TRUE,
      order = 1,
      override.aes = list(colour = "black")
    ),
    fill = guide_legend(
      title.position = "top",
      ncol = 1,
      byrow = TRUE,
      override.aes = list(
        colour = "black"
      ),
      order = 2
    )
  )

# Convert to a gtable and change the legend for the proportional to a dashed line.
# This conversion is hard coded for a set input of strategies to plot.
if (all(strat_ids == c("U1", "T1", "U2", "T2"))) {
  gg_island_vac <- ggplot_gtable(ggplot_build(gg_island_vac))
  gg_island_vac$grobs[[30]]$grobs[[2]]$grobs[[4]] <- linesGrob(
    x = unit.c(
      gg_island_vac$grobs[[30]]$grobs[[2]]$grobs[[4]]$x - gg_island_vac$grobs[[30]]$grobs[[2]]$grobs[[4]]$width / 2,
      gg_island_vac$grobs[[30]]$grobs[[2]]$grobs[[4]]$x + gg_island_vac$grobs[[30]]$grobs[[2]]$grobs[[4]]$width / 2
    ),
    y = unit(c(0.5, 0.5), "npc"),
    gp = gpar(
      col = dist$COLOUR[dist$DIST_ID == "P"],
      lty = "22",
      lwd = 4,
      lineend = "square"
    )
  )
}

# Save the figure.
file_name <- paste("fig_compare_opt_islandvac", paste(strat_ids, collapse = "_"), sep = "_")
ggsave(
  filename = paste0(file_name, ".svg"), plot = gg_island_vac, width = 9, height = 9, scale = 1.1,
  path = "figures"
)
ggsave(
  filename = paste0(file_name, ".pdf"), plot = gg_island_vac, width = 9, height = 9, scale = 1.1,
  path = "figures"
)

# Calculate some summary statistics for the main text.
stats_df <- opt_df %>%
  group_by(SENS_PAR_1, TAG_ID, DIST_ID, ISLAND_ID) %>%
  reframe(PROBS = c(0.5, 0.025, 0.975),
          PROB_NAME = sprintf("Q%04.1f", (PROBS * 100)),
          DIST = sprintf("%.2f", quantile(100 * OPT_DIST, probs = PROBS))) %>%
  select(-PROBS) %>%
  pivot_wider(names_from = PROB_NAME, values_from = DIST)
stats2_df <- prop_island_vac_df %>%
  group_by(SENS_PAR_1, TAG_ID, DIST_ID, FACET_ID) %>%
  reframe(PROBS = c(0.5, 0.025, 0.975),
          PROB_NAME = sprintf("Q%04.1f", (PROBS * 100)),
          PROP_ANIMAL_VAC = sprintf("%.2f", quantile(100 * PROP_ANIMAL_VAC, probs = PROBS))) %>%
  select(-PROBS) %>%
  pivot_wider(names_from = PROB_NAME, values_from = PROP_ANIMAL_VAC)

# Format ready for a tabular environment.
tex_tabular <- stats_df %>%
  mutate(
    LABEL = paste(
      paste0("\\texttt{", TAG_ID, DIST_ID, "}"),
      SENS_PAR_1,
      island_labels[as.character(ISLAND_ID)],
      paste0(Q50.0, "% [", Q02.5, "%, ", Q97.5, "%] \\\\"),
      sep = " & "
    ),
    LABEL = gsub("\\%", "\\\\%", LABEL),
    .keep = "none"
  )
tex_tabular <- paste(tex_tabular$LABEL, collapse = "\n")
tex_tabular2 <- stats2_df %>%
  mutate(
    ISLAND_ID = factor(FACET_ID, levels = c(1, 3, 0, 2))
  ) %>%
  arrange(SENS_PAR_1, TAG_ID, DIST_ID, ISLAND_ID) %>%
  mutate(
    LABEL = paste(
      paste0("\\texttt{", TAG_ID, DIST_ID, "}"),
      SENS_PAR_1,
      island_labels[as.character(ISLAND_ID)],
      paste0(Q50.0, "% [", Q02.5, "%, ", Q97.5, "%] \\\\"),
      sep = " & "
    ),
    LABEL = gsub("\\%", "\\\\%", LABEL),
    .keep = "none"
  )
tex_tabular2 <- paste(tex_tabular2$LABEL, collapse = "\n")
