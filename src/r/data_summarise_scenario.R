# data_summarise_scenario: Reads in a group of scenario data and produces
# some summary results of each scenario that was processed. For example,
# this script can be used to calculate infections averted
# across all time points for different sets of simulations executed.
# Pre-requisites: optimal parameters already obtained for each vaccination rate and tagging strategy,
#                 optimal parameters fed through simulation under corresponding scenarios.

# Clear the workspace.
rm(list = ls())

# Load in data manipulation libraries.
library(tidyr)
library(dplyr)
library(data.table)

# Define the input folder.
in_folder <- "../cpp/out/U1"

# Define the directory where simulations are stored.
dir <- paste0(in_folder, "/simulations")

# Get the total number of scenarios and simulation files per scenario to summarise.
max_scenarios <- length(list.files(path = dir, pattern = "scenario_pars", full.names = TRUE))
max_sims <- length(list.files(path = dir, pattern = "no_vac_summary", full.names = TRUE)) / max_scenarios

# Create a progress bar and counter
cat("Summarising simulation results...\n")
pb <- txtProgressBar(style = 3, min = 0, max = max_sims * max_scenarios)
pb_id <- 0

# Get the ID numbers of all sets of scenario parameters.
scenario_ids <- gsub("[A-za-z[:punct:]]", "", list.files(path = dir, pattern = "scenario_pars"))

# Go through every scenario...
for (scenario_id in scenario_ids) {

  # Get a list of simulation files corresponding to the currency scenario.
  files <- list.files(path = dir, pattern = paste0("no_vac_summary_", scenario_id), full.names = TRUE)

  # Create an empty data frame to store the complete vaccination and no-vaccination info.
  vac_df <- no_vac_df <- data.frame()

  # For each simulation file in the directory...
  for (file in files) {

    # Read in the simulation file (where simulations did not include vaccination)
    no_vac_df <- rbind(no_vac_df, fread(file, header = TRUE, showProgress = FALSE))

    # Read in the simulation file (where simulations included vaccination)
    vac_df <- rbind(vac_df, fread(gsub("no_vac_", "", file), header = TRUE, showProgress = FALSE))

    # Update how many simulation files have been processed.
    pb_id <- pb_id + 1
    setTxtProgressBar(pb, pb_id)
  }

  # Get the first time point where vaccines are administered.
  t_vac_start <- filter(vac_df, DOSES > 0) %>%
    head(n = 1) %>%
    .$MIN_TIME

  # Change the origin of the time scale to when vaccines were first implemented.
  no_vac_df <- no_vac_df %>%
    mutate(
      MAX_TIME = MAX_TIME - t_vac_start,
      MIN_TIME = MIN_TIME - t_vac_start,
      POST_VACCINE_TIME = (MIN_TIME >= 0)
    )
  vac_df <- vac_df %>%
    mutate(
      MAX_TIME = MAX_TIME - t_vac_start,
      MIN_TIME = MIN_TIME - t_vac_start,
      POST_VACCINE_TIME = (MIN_TIME >= 0)
    )

  # Calculate the total number of infections for the no-vaccination and vaccination scenarios.
  no_vac_df <- no_vac_df %>%
    group_by(SAMPLE_ID, REP_ID, ISLAND_ID) %>%
    reframe(
      TOTAL_INF = sum(POST_VACCINE_TIME * I_U)
    )
  vac_df <- vac_df %>%
    group_by(SAMPLE_ID, REP_ID, ISLAND_ID) %>%
    reframe(
      TOTAL_VAC_TIME = sum(POST_VACCINE_TIME * (MAX_TIME - MIN_TIME + 1)),
      TOTAL_INF = sum(POST_VACCINE_TIME * (I_U + I_V1 + I_V2 + I_W))
    )

  # Add in the total infections across the entire archipelago.
  no_vac_df <- no_vac_df %>%
    group_by(SAMPLE_ID, REP_ID) %>%
    reframe(TOTAL_INF = sum(TOTAL_INF)) %>%
    mutate(ISLAND_ID = -1) %>%
    bind_rows(no_vac_df)
  vac_df <- vac_df %>%
    group_by(SAMPLE_ID, REP_ID) %>%
    reframe(
      TOTAL_VAC_TIME = mean(TOTAL_VAC_TIME),
      TOTAL_INF = sum(TOTAL_INF)
    ) %>%
    mutate(ISLAND_ID = -1) %>%
    bind_rows(vac_df)

  # Combine the vaccination and no-vaccination scenarios.
  stats_df <- vac_df %>%
    left_join(no_vac_df,
      by = c("SAMPLE_ID", "REP_ID", "ISLAND_ID"),
      suffix = c("", "_NO_VAC")
    )

  # Create some stats, e.g. effectiveness, across the entire time period for each replicate.
  stats_df <- stats_df %>%
    group_by(SAMPLE_ID, REP_ID, ISLAND_ID) %>%
    reframe(
      TOTAL_VAC_TIME = mean(TOTAL_VAC_TIME),
      TOTAL_INF_AVERTED = TOTAL_INF_NO_VAC - TOTAL_INF,
      EFFECTIVENESS = (TOTAL_INF_AVERTED) / TOTAL_INF_NO_VAC
    )

  # Write the stats data frames to file.
  write.table(stats_df,
    sep = ",", file = paste0(dir, "/summary_", scenario_id, ".csv"),
    row.names = FALSE, col.names = TRUE
  )
}
