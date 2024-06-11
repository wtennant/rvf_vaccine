// main.cpp: describes the workflow of running through either
// a series of scenarios in the RVF Comoros model, or optimises
// the best vaccination strategy (or other strategies if implemented)
// to control RVFv in the Comoros archipelago.
#include <format>           // Nice output formatting.
#include <iostream>         // Input-output information to console.
#include "config.h"         // Number of MCMC chains to run.
#include "data.h"           // Data class definition.
#include "parameters.h"     // Parameter class definition.
#include "optimise.h"       // Optimisation class definition.
#include "posterior.h"      // Posterior class definition.
#include "scenario.h"       // Scenario class definition.

// Function describing the workflow for running a range of scenarios
// in an RVF vaccination model.
// This function has the flexibility to run the simulation for a range
// of parameter values, *BUT ONLY FOR THE GIVEN PARAMETERS*.
// Note that the default parameters may be manipulated within the scenario tests.
// These manipulations will be noted in the scenario parameters file, however.
void RunScenarios(Parameters* default_pars, Data* data, Posterior* posterior)
{
    // Write the default parameters to file.
    default_pars->WriteDefaultPars(0);

    // Define the number of parameters to perform scenarios over.
    int n_theta = 6;

    // Space to store scenario parameters to override.
    Theta** override_theta = new Theta*[n_theta];

    // Create a space to store 3-value characters.
    char str_id[4];

    // Define the number of different scenario parameters to run.
    const int n_scenarios = 24;
    int scenario_ids[n_scenarios];
    for (int i = 0; i < n_scenarios; ++i)
    {
        scenario_ids[i] = i;
    }

    // Distributions for parameters in scenarios that are to be tested.
    for (int i = 0; i < n_scenarios; ++i)
    {
        // Format the IDs for the custom scenario parameters.
        sprintf(str_id, "%.3d", scenario_ids[i]);
        
        // Create the custom scenario parameters.
        Empirical vac_rate{"vac_rate", "optimal_" + static_cast<std::string>(str_id) + ".csv", 0};     // Vaccination rate (vaccines adminsitered per year).
        Empirical vac_dist{"vac_dist", "optimal_" + static_cast<std::string>(str_id) + ".csv", 2};     // Vaccine distribution between islands.
        Empirical vac_protect_duration{"vac_protect_duration", "optimal_" + static_cast<std::string>(str_id) + ".csv", 3};
        Empirical vac_identifiable{"vac_identifiable", "optimal_" + static_cast<std::string>(str_id) + ".csv", 9};

        // Store the distributions into the overriding theta distribution.
        override_theta[0] = &vac_rate;
        override_theta[1] = &vac_dist;
        override_theta[2] = &vac_protect_duration;
        override_theta[3] = &vac_identifiable;

        // Define the number and types of scenarios to run.
        const int n_n_age_vac = 2;
        int n_age_vacs[n_n_age_vac] = {2, 10};
        const int n_vac_t_year_length = 2;
        int vac_t_year_lengths[n_vac_t_year_length] = {4, 48};

        // For each desired vaccine frequency to run scenarios for...
        for (int n_age_vac_idx = 0; n_age_vac_idx < n_n_age_vac; ++n_age_vac_idx)
        {
            // Set up a distribution at the desired scenario distribution.
            DiscreteUniform n_age_vac{"n_age_vac",  n_age_vacs[n_age_vac_idx], n_age_vacs[n_age_vac_idx], 6};

            // Store the distribution into the overriding scenario distribution.
            override_theta[4] = &n_age_vac;

            // For each vaccine age group to run scenarios for...
            for (int vac_t_year_length_idx = 0; vac_t_year_length_idx < n_vac_t_year_length; ++vac_t_year_length_idx)
            {
                // Set up the dsitribution for maximum age group to vaccinate.
                DiscreteUniform vac_t_year_length{"vac_t_year_length", vac_t_year_lengths[vac_t_year_length_idx], vac_t_year_lengths[vac_t_year_length_idx], 5};

                // Store the distribution into the overriding scenario distribution.
                override_theta[5] = &vac_t_year_length;

                // Define an index for the scenario.
                int scen_idx = i*n_n_age_vac*n_vac_t_year_length + n_age_vac_idx*n_vac_t_year_length + vac_t_year_length_idx;

                // Output some information for the user.
                std::cout << std::format("\nRunning scenario {:03d} of {:03d}...\n", scen_idx + 1, n_scenarios*n_n_age_vac*n_vac_t_year_length);

                // Setup the scenario parameter file.
                Scenario scenario(default_pars, override_theta, n_theta, scen_idx);

                // Run the scenario tests.
                scenario.Run(data, posterior);
            }
        }
    }
}

// Function describing the workflow for optimising a set of vaccination strategies
// in the RVF vaccination model.
// This function can be customised in order to run through a range of scenarios
// over which to optimise.
// Currently the optimisation is fixed to determine the best way to distribute vaccines across
// the four islands.
void RunOptimisation(Parameters* default_pars, Data* data, Posterior* posterior)
{
    // Define the settings for running a series of optimisations.
    // Discretisation of the first parameter to perform optimisations over.    
    double start_par = 0.05; double end_par = 0.3; int n_series = 6;
    for (int series_id = 0; series_id < n_series; ++series_id)
    {
        // The first parameter to get the optimal vaccine distribution for is...
        default_pars->vac_rate = (n_series > 1) ? start_par + series_id*(end_par - start_par) / static_cast<double>(n_series - 1) : start_par;
        default_pars->vac_rate = default_pars->vac_rate*default_pars->CalculateTotalPop();

        // Store the set of default parameters used in the optimisation.
        default_pars->WriteDefaultPars(series_id);

        // Setup space for the optimisation for the desired parameter.
        Optimise optimise(default_pars, data, posterior);

        // Run the scenario tests.
        optimise.Run(series_id);
    }  
}

// Function which is the first to be invoked.
int main()
{
    // Output the configuration of the model to file.
    WriteConfig();

    // Setup the default parameters for the model.
    Parameters default_pars;
    
    // Set up the data for the model.
    // Adjust the number of time steps in the set of default parameters
    // to the number of time steps that appear in the NDVI data.
    Data data(&default_pars);

    // Set up the posterior distribution to run scenarios from.
    // This is where one source of uncertainty comes from in the model.
    // The posterior distribution came from fitting the original model
    // to serosurvey data (Tennant, 2021).
    Posterior posterior;

    // Set up the space for scenario testing or optimisation.
    switch(RVF_OPTIMISE)
    {
        case 0: // Scenario testing.
            RunScenarios(&default_pars, &data, &posterior);
            break;
        case 1: // Optimisation.
            RunOptimisation(&default_pars, &data, &posterior);
            break;
        default: // Run scenarios by default.
            RunScenarios(&default_pars, &data, &posterior);
    }

    // Return an error-free code.
    return 0;
}

