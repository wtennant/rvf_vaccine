# data_format_optimal.R: Format the parameter sets obtained from optimisation
# of a particular scenario. These optimal parameter sets are then arranged in a
# data frame of four columns: Sample ID, parameter name, parameter index and
# parameter value.
# Pre-requisites: optimisation algorithm executed for a given vaccination rate and tagging strategy.

# Clear the workspace.
rm(list = ls())

# Load in data manipulation libraries.
library(dplyr)
library(tidyr)
library(data.table)

# Define the input folder.
in_folder <- "in/"

# List of default parameters to add.
default_par_list <- c("vac_rate", "vac_identifiable", "vac_protect_duration")

# Define a folder to output the data summaries to.
out_folder <- "out/"

# Get the total number of files in the folder.
sens_par_opts <- as.numeric(gsub("[A-z.]", "", list.files(path = in_folder, pattern = "cost_"))) + 1

# For each optimal distribution...
for (i in sens_par_opts) {
  # Load in the optimal parameter set.
  opt_df <- fread(paste0(in_folder, "optim_pars_", sprintf("%.3d", i - 1), ".csv"),
    header = TRUE, data.table = FALSE
  )

  # Load in the associated cost data.
  cost_df <- fread(paste0(in_folder, "cost_", sprintf("%.3d", i - 1), ".csv"),
    header = TRUE, data.table = FALSE
  )

  # For each optimisation, calculate the mean cost per particle, and then
  # select only the particle with the maximum mean cost.
  max_cost_df <- group_by(cost_df, OPTIM_ID, PARTICLE_ID) %>%
    filter(STEP_ID == max(STEP_ID)) %>%
    summarise(MEAN_COST = mean(COST_VALUE), .groups = "drop") %>%
    group_by(OPTIM_ID) %>%
    filter(MEAN_COST == max(MEAN_COST)) %>%
    sample_n(size = 1)

  # Use the best cost for each optimisation run to find the optimal parameters.
  opt_par_df <- opt_df %>%
    group_by(OPTIM_ID, PARTICLE_ID) %>%
    filter(STEP_ID == max(STEP_ID)) %>%
    right_join(max_cost_df, by = c("OPTIM_ID", "PARTICLE_ID")) %>%
    ungroup() %>%
    select(OPTIM_ID, PAR_NAME, INDEX, PAR_VALUE) %>%
    rename(PAR_INDEX = INDEX, SAMPLE_ID = OPTIM_ID)

  # Load in the default parameter set and record the vaccination rate that was used.
  default_par_df <- fread(paste0(in_folder, "default_pars_", sprintf("%.3d", i - 1), ".csv"),
    header = TRUE, data.table = FALSE
  ) %>%
    filter(PAR_NAME %in% default_par_list) %>%
    expand(SAMPLE_ID = opt_par_df$SAMPLE_ID, nesting(PAR_NAME, PAR_INDEX, PAR_VALUE))

  # Ammend this information to the optimal parameter data frame
  # and store the unique optimisation scenario ID.
  opt_par_df <- opt_par_df %>%
    bind_rows(default_par_df) %>%
    arrange(SAMPLE_ID, PAR_NAME, PAR_INDEX)

  # Store the the optimal parameter data frame as a csv file
  # ready for reading in and running more scenarios on top of the optimal distributions.
  write.table(opt_par_df,
    sep = ",",
    file = paste0(out_folder, "optimal_", sprintf("%.3d", i - 1), ".csv"),
    col.names = TRUE, row.names = FALSE
  )
}
