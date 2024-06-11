// Scenario.h: Defines the class Scenario used to run different scenario tests
// given the posterior from a previous model fit.
#include "data.h"           // Definition of Data class.
#include "posterior.h"      // Definition of the Posterior class.
#include "simulation.h"     // Definition of Simulation class.
#include "theta.h"          // Definition of MCMC parameter class.

#ifndef SCENARIO_H
#define SCENARIO_H

class Scenario
{
public:
    Scenario(Parameters* default_pars, Theta**, int, int);          // Default constructor.
    ~Scenario();                                                    // Default deconstructor.

    // Function to start running the scenario tests.
    void Run(Data* data, Posterior* post);

private:
    // Declare space to store a series of simulations and parameters.
    Simulation** sim;
    Parameters** pars;

    // Declare a pointer to a series of random number generators.
    gsl_rng** rng;

    // Define the number of simulations to run simultaneously.
    int n_sims;

    // Define a copy of the number of islands.
    int n_island;

    // Sample IDs for the posteriors used.
    int* post_idx;

    // Number of scenario-related parameters which are sampled at the same time.
    int n_scenarios;

    // These parameters are transformed versions of the scenario parameters.
    double* vac_prop;

    // Space to store all scenario and default scenario parameters.
    Theta** default_theta;
    Theta** theta;
    
    // Unique ID number associated with the scenario class.
    int scenario_id;

    // Define the list of possible scenario parameters.
    Theta* vac_rate;
    Theta* vac_protect_duration;
    Theta* vac_dist;
    Theta* vac_efficacy;
    Theta* vac_t_year_offset;
    Theta* vac_t_year_length;
    Theta* n_age_vac;
    Theta* vac_t_freq;
    Theta* reg_import_scale;
    Theta* vac_identifiable;
    
    // Functions to carry out the scenario tests.
    void SetSimPars(int sample_id, Posterior* post, int sim_id);        // Set the parameters of the simulation from the posterior distribution.
    void UpdateSimPars(int sample_id, int sim_id);                      // Update the simulation parameter sets based on sampled scenarios.
    void PrintProgress(int sample_id);                                  // Print progress of scenario tests to console.
    void WritePosteriorSamples(Posterior* post);                        // Output sampled posterior parameters to file.
    void WriteScenarioPars();                                           // Output scenario parameters to file.
};

#endif