# plot_compare_opt_effect: plotting to compare the effectiveness of
# different optimal vaccine strategies across the archipelago and on each
# individual island.
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
library(egg)

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
strat_ids <- c("UP", "U1", "U2", "TP", "T1", "T2")

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
    "Optimised over infections averted across\nthe archipelago",
    "Optimised over infections averted on the\nworst-peforming island"
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
    FACET_ID = factor(ISLAND_ID, levels = c(1, 0, 3, 2, -1)),
    SENS_PAR_1 = factor(sens_par_label_1(SENS_PAR_1))
  ) %>%
  left_join(
    all_strats,
    by = "STRAT_ID"
  )
summary_df$SENS_PAR_1 <- factor(summary_df$SENS_PAR_1,
  levels = paste0(sort(as.numeric(gsub("%", "", levels(summary_df$SENS_PAR_1)))), "%")
)

# Map island indices to names.
island_labels <- c(
  "0" = "Anjouan",
  "1" = "Grande Comore",
  "2" = "Mayotte",
  "3" = "Mohéli",
  "-1" = "Comoros archipelago"
)

# Calculate the worst affected island on average (by mode of worst performer
# across all simulations in the group).
worst_island_id_df <- summary_df %>%
  filter(ISLAND_ID != -1) %>%
  group_by(GROUP_ID, SAMPLE_ID, REP_ID) %>%
  mutate(MIN_EFFECTIVENESS = min(EFFECTIVENESS)) %>%
  filter(near(MIN_EFFECTIVENESS, EFFECTIVENESS)) %>%
  group_by(GROUP_ID) %>%
  summarise(AVG_WORST_ISLAND_ID = factor(c(1, 3, 0, 2, -1)[which.max(tabulate(ISLAND_ID))],
    levels = c(1, 3, 0, 2, -1)
  ))

# What is the effectiveness on the worst affected islands?
worst_island_df <- worst_island_id_df %>%
  right_join(summary_df, by = "GROUP_ID", multiple = "all") %>%
  filter(ISLAND_ID == AVG_WORST_ISLAND_ID)

# What are the maximum and minimum effectiveness of every group and island?
highlight_worst_island_df <- worst_island_id_df %>%
  right_join(summary_df, by = "GROUP_ID", multiple = "all") %>%
  filter(ISLAND_ID != -1) %>%
  group_by(GROUP_ID, ISLAND_ID, FACET_ID, SENS_PAR_1, STRAT_ID, AVG_WORST_ISLAND_ID) %>%
  summarise(
    MIN_EFF = min(EFFECTIVENESS),
    MAX_EFF = max(EFFECTIVENESS), .groups = "drop"
  ) %>%
  group_by(GROUP_ID) %>%
  mutate(IS_WORST_ISLAND = FACET_ID == AVG_WORST_ISLAND_ID)

# Visualise the infections averted on each island over different vaccination rates.
gg_cost <- summary_df %>%
  filter(ISLAND_ID != -1) %>%
  ggplot() +
  geom_violin(
    mapping = aes(
      x = SENS_PAR_1,
      y = EFFECTIVENESS,
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
      y = EFFECTIVENESS,
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
  geom_rect(
    data = highlight_worst_island_df,
    mapping = aes(
      xmin = as.numeric(SENS_PAR_1) - 0.85 / 2,
      xmax = as.numeric(SENS_PAR_1) + 0.85 / 2,
      ymin = MIN_EFF - 0.01,
      ymax = MAX_EFF + 0.01,
      group = GROUP_ID,
      linetype = IS_WORST_ISLAND
    ),
    fill = NA,
    colour = "black",
    position = position_dodge(width = 0.9)
  ) +
  geom_vline(
    data = ~ data.frame(X = seq_len(nlevels(.$SENS_PAR_1) - 1)),
    mapping = aes(xintercept = X + 0.5),
    colour = "black",
    alpha = 0.25,
    linetype = "11"
  ) +
  scale_x_discrete(sens_par_title_1) +
  scale_y_continuous("Percentage of infections averted",
    label = scales::percent,
    breaks = seq(0, 1, by = 0.25),
    limits = c(NA, 1.05)
  ) +
  scale_linetype_manual("Summary statistic",
    breaks = c("TRUE", "FALSE"),
    labels = c("Island with the\nworst performance", ""),
    values = c("TRUE" = "22", "FALSE" = NA)
  ) +
  scale_fill_manual(
    name = "Vaccine distribution",
    limits = dist$DIST_ID,
    breaks = dist$DIST_ID,
    labels = dist$DIST_LABEL,
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
  facet_wrap(~ISLAND_ID, labeller = labeller(ISLAND_ID = island_labels), ncol = 1) +
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
    ),
    linetype = guide_legend(title.position = "top", nrow = 2, order = 3)
  )

# Save the figure: effectiveness on each island.
file_name <- paste("fig_compare_opt_effect_island", paste(strat_ids, collapse = "_"), sep = "_")
ggsave(
  filename = paste0(file_name, ".svg"), plot = gg_cost, width = 9, height = 12, scale = 1.1,
  path = "figures"
)
ggsave(
  filename = paste0(file_name, ".pdf"), plot = gg_cost, width = 9, height = 12, scale = 1.1,
  path = "figures"
)

# And for the entire archipelago...
gg_cost_ca <- summary_df %>%
  filter(ISLAND_ID == -1) %>%
  ggplot() +
  geom_violin(
    mapping = aes(
      x = SENS_PAR_1,
      y = EFFECTIVENESS,
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
      y = EFFECTIVENESS,
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
    linewidth = 1
  ) +
  geom_vline(
    data = ~ data.frame(X = seq_len(nlevels(.$SENS_PAR_1) - 1)),
    mapping = aes(xintercept = X + 0.5),
    colour = "black",
    alpha = 0.25,
    linetype = "11"
  ) +
  scale_x_discrete(sens_par_title_1) +
  scale_y_continuous("Percentage of infections averted\nacross the Comoros archipelago",
    label = scales::percent,
    breaks = seq(0, 1, by = 0.25),
    limits = c(NA, 1.05)
  ) +
  scale_fill_manual(
    name = "Vaccine distribution",
    limits = dist$DIST_ID,
    breaks = dist$DIST_ID,
    labels = dist$DIST_LABEL,
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
  gg_theme +
  theme(
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

# Save the figure: effectiveness on each island.
file_name <- paste("fig_compare_opt_effect_archipelago", paste(strat_ids, collapse = "_"), sep = "_")
ggsave(
  filename = paste0(file_name, ".svg"), plot = gg_cost_ca, width = 9, height = 8, scale = 1.1,
  path = "figures"
)
ggsave(
  filename = paste0(file_name, ".pdf"), plot = gg_cost_ca, width = 9, height = 8, scale = 1.1,
  path = "figures"
)

# Calculate the summary statistics for each strategy.
stats_df <- summary_df %>%
  mutate(TOTAL_INF_NO_VAC = TOTAL_INF_AVERTED / EFFECTIVENESS) %>%
  mutate(TOTAL_INF = TOTAL_INF_NO_VAC - TOTAL_INF_AVERTED) %>%
  group_by(SENS_PAR_1, TAG_ID, DIST_ID, ISLAND_ID) %>%
  reframe(PROBS = c(0.5, 0.025, 0.975),
          PROB_NAME = sprintf("Q%04.1f", (PROBS * 100)),
          INF_NOVAC = sprintf("%.0f", quantile(TOTAL_INF_NO_VAC, probs = PROBS)),
          INF_VAC = sprintf("%.0f", quantile(TOTAL_INF, probs = PROBS)),
          EFF = sprintf("%.2f", quantile(100 * EFFECTIVENESS, probs = PROBS))) %>%
  select(-PROBS) %>%
  pivot_wider(names_from = PROB_NAME, values_from = INF_NOVAC:EFF)

# Format ready for a tabular environment.
tex_tabular_ca <- stats_df %>%
  filter(ISLAND_ID == -1) %>%
  mutate(
    LABEL = paste(
      paste0("\\texttt{", TAG_ID, DIST_ID, "}"),
      SENS_PAR_1,
      paste0(INF_VAC_Q50.0, " [", INF_VAC_Q02.5, ", ", INF_VAC_Q97.5, "]"),
      paste0(EFF_Q50.0, "% [", EFF_Q02.5, "%, ", EFF_Q97.5, "%] \\\\"),
      sep = " & "
    ),
    LABEL = gsub("\\%", "\\\\%", LABEL),
    .keep = "none"
  )
tex_tabular_ca <- paste(tex_tabular_ca$LABEL, collapse = "\n")
tex_tabular_inf <- stats_df %>%
  filter(ISLAND_ID != -1) %>%
  mutate(
    LABEL = paste(
      paste0("\\texttt{", TAG_ID, DIST_ID, "}"),
      SENS_PAR_1,
      island_labels[as.character(ISLAND_ID)],
      paste0(INF_VAC_Q50.0, " [", INF_VAC_Q02.5, ", ", INF_VAC_Q97.5, "] \\\\"),
      sep = " & "
    ),
    LABEL = gsub("\\%", "\\\\%", LABEL),
    .keep = "none"
  )
tex_tabular_inf <- paste(tex_tabular_inf$LABEL, collapse = "\n")
tex_tabular_eff <- stats_df %>%
  filter(ISLAND_ID != -1) %>%
  mutate(
    LABEL = paste(
      paste0("\\texttt{", TAG_ID, DIST_ID, "}"),
      SENS_PAR_1,
      island_labels[as.character(ISLAND_ID)],
      paste0(EFF_Q50.0, "% [", EFF_Q02.5, "%, ", EFF_Q97.5, "%] \\\\"),
      sep = " & "
    ),
    LABEL = gsub("\\%", "\\\\%", LABEL),
    .keep = "none"
  )
tex_tabular_eff <- paste(tex_tabular_eff$LABEL, collapse = "\n")
tex_tabular_no_vac <- summary_df %>%
  mutate(TOTAL_INF_NO_VAC = TOTAL_INF_AVERTED / EFFECTIVENESS) %>%
  mutate(TOTAL_INF = TOTAL_INF_NO_VAC - TOTAL_INF_AVERTED) %>%
  group_by(ISLAND_ID) %>%
  reframe(PROBS = c(0.5, 0.025, 0.975),
          PROB_NAME = sprintf("Q%04.1f", (PROBS * 100)),
          INF_NOVAC = sprintf("%.0f", quantile(TOTAL_INF_NO_VAC, probs = PROBS)),
          INF_VAC = sprintf("%.0f", quantile(TOTAL_INF, probs = PROBS)),
          EFF = sprintf("%.2f", quantile(100 * EFFECTIVENESS, probs = PROBS))) %>%
  select(-PROBS) %>%
  pivot_wider(names_from = PROB_NAME, values_from = INF_NOVAC:EFF) %>%
  mutate(
    LABEL = paste(
      island_labels[as.character(ISLAND_ID)],
      paste0(INF_NOVAC_Q50.0, " [", INF_NOVAC_Q02.5, ", ", INF_NOVAC_Q97.5, "] \\\\"),
      sep = " & "
    ),
    LABEL = gsub("\\%", "\\\\%", LABEL),
    .keep = "none"
  )
tex_tabular_no_vac <- paste(tex_tabular_no_vac$LABEL, collapse = "\n")
