 // parameters.h: Defines the class named Parameters which
// are used with in the model.
#include <string>   // Using strings.

#ifndef PARAMETERS_H
#define PARAMETERS_H

class Parameters
{
public:
    int n_age;                      // Number of livestock age-groups (grouped by year).
    int n_age_vac;                  // Maximum age at which to vaccinate.
    int n_age_move;                 // Maximum age group which can be moved.
    int n_island;                   // Number of islands.
    int n_steps;                    // Number of time steps to simulate (epi-weeks).
    double import_start;            // Start time of import into Grande Comore.
    double import_duration;         // Duration of import into Grande Comore.
    double import_size;             // Number of infectious imports into Grande Comore per year.
    double* n_pop;                  // Population size of each island.
    double* p_immune_start;         // Proportion of animals that are immune at time zero on each island.
    double* trans_scale;            // Scalar for transmission rate on each island.
    double* ndvi_rate;              // Scalar for how NDVI changes force of infection for each island.
    double** move;                  // Number of annual movements between each island.
    double* pop_structure;          // Population structure (proportion of population in
                                    // each age group).
    double* mortality;              // Age-specific mortality rates (proportion per epi-week).
    int vac_identifiable;           // Whether (> 0) or not (== 0) only unvaccinated animals are vaccinated.
    double vac_efficacy;            // The efficacy of the vaccine wrt infection prevention.
    double vac_rate;                // Vaccination rate across the archipelago.
    double* vac_prop;               // Distribution of vaccines across the archipelago (per island).
    double vac_protect_duration;    // Duration of protection once vaccinated
    int vac_t_first;                // The first time point at which the vaccine should be administered.
    int vac_t_year_offset;          // The first time point of each epidemiological year which to give the vaccine.
    int vac_t_year_length;          // The length of the vaccine time window in epi-weeks.
    int vac_t_freq;                 // Frequency of the vaccination window (in years).
    int reg_import_freq;            // Frequiency of external imports into Grande Comore.
    int reg_import_t_start;         // The first time point where regular imports can be introduced into Grande Comore.
    double reg_import_scale;        // Scalar for the number of regular imports to make per week (of import_size).

    // Function to compute the total population size.
    double CalculateTotalPop();

    // Function which writes the default parameters to a file appended with an integer.
    void WriteDefaultPars(int);    

    // Constructor and destructor.
    Parameters();
    ~Parameters();

    // Constructor for a deep copy implementation.
    Parameters(const Parameters& default_pars);
        
};

#endif