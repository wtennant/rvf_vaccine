# data_summarise_scenario_time: Reads in a group of scenario data and produces
# some summary results of each scenario that was processed. For example,
# this script can be used to calculate the vaccine efficiency over time.
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
dir <- paste0(in_folder, "/weekly_simulations")

# Get the total number of scenarios and simulation files per scenario to summarise.
max_scenarios <- length(list.files(path = dir, pattern = "scenario_pars", full.names = TRUE))
max_sims <- length(list.files(path = dir, pattern = "no_vac_summary", full.names = TRUE)) / max_scenarios

# Read in the NDVI data and get the minimum NDVI for each island.
ndvi <- fread("../data/ndvi_extended.csv") %>%
  group_by(ISLAND_ID) %>%
  mutate(
    MIN_NDVI = min(NDVI),
    TIME = seq_len(n()) - 1
  )

# Create a progress bar and counter
cat("Summarising simulation results...\n")
pb <- txtProgressBar(style = 3, min = 0, max = max_sims * max_scenarios)
pb_id <- 0

# Get the ID numbers of all sets of scenario parameters.
scenario_ids <- gsub("[A-za-z[:punct:]]", "", list.files(path = dir, pattern = "scenario_pars"))

# Go through every scenario...
for (scenario_id in scenario_ids) {

  # Read in the posterior information.
  post <- fread(paste0(dir, "/posterior_", scenario_id, ".csv")) %>%
    pivot_wider(id_cols = c(SAMPLE_ID, REP_ID, POST_ID), names_from = PAR_NAME, values_from = PAR_VALUE) %>%
    pivot_longer(
      cols = trans_scale_anj:trans_scale_moh,
      names_to = "TRANS_SCALE_NAME", values_to = "TRANS_SCALE_VALUE"
    ) %>%
    group_by(SAMPLE_ID, REP_ID) %>%
    mutate(ISLAND_ID = 0:3) %>%
    rename(
      NDVI_RATE = ndvi_rate,
      TRANS_SCALE = TRANS_SCALE_VALUE
    ) %>%
    select(SAMPLE_ID, REP_ID, ISLAND_ID, NDVI_RATE, TRANS_SCALE)

  # Load in the sceanrio parameters.
  scen_pars <- fread(paste0(dir, "/scenario_pars_", scenario_id, ".csv")) %>%
    unite(col = "PAR_NAME", PAR_NAME, INDEX) %>%
    pivot_wider(id_cols = SAMPLE_ID, names_from = PAR_NAME, values_from = PAR_VALUE) %>%
    rename(VAC_EFF = vac_efficacy_0, VAC_IDENT = vac_identifiable_0) %>%
    select(SAMPLE_ID, VAC_EFF, VAC_IDENT)

  # Get a list of simulation files corresponding to the currency scenario.
  files <- list.files(path = dir, pattern = paste0("no_vac_summary_", scenario_id), full.names = TRUE)

  # Create an empty data frame to store the vaccine and vaccine-free scenarios.
  vac_df <- no_vac_df <- data.frame()

  # For each simulation file in the directory...
  for (file in files) {
    # Read in the simulation file (where simulations did not include vaccination)
    no_vac_df <- rbind(no_vac_df, fread(file, header = TRUE, showProgress = FALSE, nThread = 1))

    # Read in the simulation file (where simulations included vaccination)
    vac_df <- rbind(vac_df, fread(gsub("no_vac_", "", file), header = TRUE, showProgress = FALSE, nThread = 1))

    # Update how many simulation files have been processed.
    pb_id <- pb_id + 1
    setTxtProgressBar(pb, pb_id)
  }

  # Retrieve when doses were first adminsitered.
  first_vac_t <- vac_df %>%
    filter(DOSES > 0) %>%
    filter(MIN_TIME == min(MIN_TIME)) %>%
    head(n = 1) %>%
    .$MIN_TIME

  # Join in the scenario parameters, NDVI data and posterior values.
  no_vac_df <- no_vac_df %>%
    left_join(scen_pars, by = c("SAMPLE_ID"), relationship = "many-to-one") %>%
    left_join(post, by = c("SAMPLE_ID", "REP_ID", "ISLAND_ID"), relationship = "many-to-one") %>%
    left_join(ndvi, by = c("MIN_TIME" = "TIME", "ISLAND_ID"), relationship = "many-to-one")
  vac_df <- vac_df %>%
    left_join(scen_pars, by = c("SAMPLE_ID"), relationship = "many-to-one") %>%
    left_join(post, by = c("SAMPLE_ID", "REP_ID", "ISLAND_ID"), relationship = "many-to-one") %>%
    left_join(ndvi, by = c("MIN_TIME" = "TIME", "ISLAND_ID"), relationship = "many-to-one")

  # Filter out the final time point and any time points prior to vaccination.
  no_vac_df <- no_vac_df %>%
    filter(MIN_TIME <= max(ndvi$TIME)) %>%
    mutate(TIME = MIN_TIME - first_vac_t)
  vac_df <- vac_df %>%
    filter(MIN_TIME <= max(ndvi$TIME)) %>%
    mutate(TIME = MIN_TIME - first_vac_t)

  # Calculate the basic and effective reproduction number.
  no_vac_df <- no_vac_df %>%
    mutate(
      RNAUGHT = exp(NDVI_RATE * (NDVI - MIN_NDVI) + TRANS_SCALE),
      N = rowSums(across(.cols = S_U:R_U)),
      PROP_S = S_U / N,
      REFF = RNAUGHT * PROP_S,
      INF = I_U
    )
  vac_df <- vac_df %>%
    mutate(
      RNAUGHT = exp(NDVI_RATE * (NDVI - MIN_NDVI) + TRANS_SCALE),
      N = rowSums(across(.cols = S_U:R_W)),
      PROP_S = (S_U + S_V1 + S_V2 + (1 - VAC_EFF) * S_W) / N,
      REFF = RNAUGHT * PROP_S,
      VACCINATABLE = ifelse(VAC_IDENT < 0.5, N, rowSums(across(.cols = S_U:R_U))),
      VAC_EFFICIENCY = S_U / VACCINATABLE,
      INF = I_U + I_V1 + I_V2 + I_W
    )

  # Combine the two data frames.
  r_df <- no_vac_df %>%
    left_join(
      vac_df,
      by = c("SAMPLE_ID", "REP_ID", "ISLAND_ID", "DATE"),
      suffix = c("_NO_VAC", "")
    ) %>%
    mutate(INF_AVERTED = INF_NO_VAC - INF)

  # Calculate the time-average (geometric mean) of
  # each statistic of interest.
  sum_r_df <- r_df %>%
    filter(TIME >= 0) %>%
    group_by(SAMPLE_ID, REP_ID, ISLAND_ID) %>%
    reframe(
      across(
        .cols = c(
          contains("RNAUGHT"),
          contains("REFF")
        ),
        .fns = \(x) exp(mean(log(x))),
        .names = "{.col}_GEOMEAN"
      ),
      across(
        .cols = c(
          contains("PROP_S"),
          contains("VAC_EFFICIENCY")
        ),
        .fns = \(x) mean(x),
        .names = "{.col}_MEAN"
      )
    )

  # Summarise the statistics for each time point and island.
  r_df <- r_df %>%
    group_by(DATE, TIME, ISLAND_ID) %>%
    reframe(
      QUAN = c(0.5, 0.25, 0.75, 0.025, 0.975),
      INT = paste0("CrI", sprintf("%02.0f", 100 * abs(QUAN - 0.5) * 2)),
      INT = ifelse(INT == "CrI00", "Median", INT),
      SIDE = ifelse(INT == "Median", "Median", ifelse(QUAN < 0.5, "Lower", "Upper")),
      across(
        .cols = c(
          contains("RNAUGHT"),
          contains("PROP_S"),
          contains("REFF"),
          contains("VAC_EFFICIENCY"),
          contains("INF")
        ),
        .fns = \(x) quantile(x, probs = QUAN)
      )
    )

  # Write the stats data frames to file.
  fwrite(r_df,
    sep = ",", file = paste0(dir, "/summary_time_", scenario_id, ".csv"),
    row.names = FALSE, col.names = TRUE
  )
  fwrite(sum_r_df,
    sep = ",", file = paste0(dir, "/summary_", scenario_id, ".csv"),
    row.names = FALSE, col.names = TRUE
  )
}
