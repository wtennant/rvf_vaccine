// scenario.cpp: Defines the function to run all the scenario tests.
#include <chrono>               // Time used as the random number generator seed.
#include <fstream>              // Writing to file.
#include <iostream>             // Output information to console.
#include <iomanip>              // Set precision of file output.
#include "omp.h"                // Parallelisation API.
#include "gsl/gsl_randist.h"    // Distributions to evaluate likelihood.
#include "config.h"             // Configuration of file output.
#include "scenario.h"           // Definition of the Scenario class.

// Constructor of the Scenario class.
Scenario::Scenario(Parameters* default_pars, Theta** override_theta, int n_theta, int id)
{
    // Set the scenario ID.
    scenario_id = id;

    // Define the number of parameters involed in the scenario testing.
    n_scenarios = 10;    // The number of parameters independent of scenario choice.

    // Setup the scenario parameters to default parameter values initially.
    vac_rate = new Uniform{"vac_rate", default_pars->vac_rate, default_pars->vac_rate, 0};
    vac_protect_duration = new Uniform{"vac_protect_duration", default_pars->vac_protect_duration,
                                       default_pars->vac_protect_duration, 1};
    vac_dist = new Dirichlet{"vac_dist", 1, 2};
    vac_efficacy = new Uniform{"vac_efficacy", default_pars->vac_efficacy, default_pars->vac_efficacy, 3};                    
    vac_t_year_offset = new DiscreteUniform{"vac_t_year_offset", default_pars->vac_t_year_offset,
                                            default_pars->vac_t_year_offset, 4};  
    vac_t_year_length = new DiscreteUniform{"vac_t_year_length", default_pars->vac_t_year_length,
                                            default_pars->vac_t_year_length, 5};
    n_age_vac = new DiscreteUniform{"n_age_vac", default_pars->n_age_vac, default_pars->n_age_vac, 6};
    vac_t_freq = new DiscreteUniform{"vac_t_freq", default_pars->vac_t_freq, default_pars->vac_t_freq, 7};
    reg_import_scale = new Uniform{"reg_import_scale", default_pars->reg_import_scale, default_pars->reg_import_scale, 8}; 
    vac_identifiable = new DiscreteUniform{"vac_identifiable", default_pars->vac_identifiable, default_pars->vac_identifiable, 9};

    // Set up a vector of default scenario parameter distributions, and one
    // that is actually used.
    default_theta = new Theta*[n_scenarios];
    theta = new Theta*[n_scenarios];

    // To the default vector of distributions, assign the addresses
    // of the default distributions in a specific order.
    default_theta[0] = vac_rate;
    default_theta[1] = vac_protect_duration;
    default_theta[2] = vac_dist;
    default_theta[3] = vac_efficacy;
    default_theta[4] = vac_t_year_offset;
    default_theta[5] = vac_t_year_length;
    default_theta[6] = n_age_vac;
    default_theta[7] = vac_t_freq;
    default_theta[8] = reg_import_scale;
    default_theta[9] = vac_identifiable;

    // Loop through each scenario parameter, and check if there is an override.
    // Note we are reliant that the overriding thetas are not deallocated prior to the 
    // usage of the scenario class.
    // This can be ammended using a deep-copy of the overriding thetas.
    for (int i = 0; i < n_scenarios; ++i)
    {
        // Check the parameter name against the current theta parameter names.
        bool is_matching_theta = false;
        int j = 0;

        // Until a match is found...
        while ((!is_matching_theta) && (j < n_theta))
        {
            // Check the parameter names against one another.
            if (override_theta[j]->par_name == default_theta[i]->par_name)
            {
                // If they match, deallocate the space of the default parameter set,
                // and assign the address of the overriding distribution.
                is_matching_theta = true;
                theta[i] = override_theta[j];
            }
            ++j;
        }

        // If no match was found in the overriding theta, then use the default.
        if (!is_matching_theta)
        {
            theta[i] = default_theta[i];
        }
    }
    
    // Hard copy back the addresses stored in theta to their corresponding parameter variables.
    vac_rate = theta[0];
    vac_protect_duration = theta[1];
    vac_dist = theta[2];
    vac_efficacy = theta[3];
    vac_t_year_offset = theta[4];
    vac_t_year_length = theta[5];
    n_age_vac = theta[6];
    vac_t_freq = theta[7];
    reg_import_scale = theta[8];
    vac_identifiable = theta[9];

    // Copy across the number of islands in the simulation.
    n_island = default_pars->n_island;

    // Define space for the proportion of vaccines per island.
    vac_prop = new double[n_island*RVF_NSAMPLES];

    // Randomly sample and store all scenario parameters.
    for (int i = 0; i < RVF_NSAMPLES; ++i)
    {
        for (int j = 0; j < n_scenarios; ++j)
        {
            theta[j]->Sample(i);
        }

         // Set the default proportion of vaccines per island.
         for (int j = 0; j < n_island; ++j)
         {
            vac_prop[i*n_island + j] = default_pars->vac_prop[j];
         }         
    }

    // Declare space for the IDs of posterior samples.
    post_idx = new int[RVF_NSAMPLES*RVF_POST_NSAMPLES];

    // Define the number of simulations to run simultaneously.
    n_sims = RVF_NTHREADS;

    // Set up the random number generator based on time.
    rng = new gsl_rng*[n_sims];
    for (int sim_id = 0; sim_id < n_sims; ++sim_id)
    {
        rng[sim_id] = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng[sim_id], std::chrono::system_clock::now().time_since_epoch().count() + sim_id);
    }

    // Allocate space for the model parameters.
    // Use the set of default parameters to deep copy parameter values across.
    pars = new Parameters*[n_sims];
    for (int sim_id = 0; sim_id < n_sims; ++sim_id)
    {
        pars[sim_id] = new Parameters(*default_pars);
    }

    // Allocate space for running the simulations.
    sim = new Simulation*[n_sims];
    for (int sim_id = 0; sim_id < n_sims; ++sim_id)
    {
        sim[sim_id] = new Simulation(pars[sim_id], scenario_id, sim_id);
    }
}

// Deconstructor of the Scenario class.
Scenario::~Scenario()
{
    // De-allocate the random number generator, simulation space,
    // and parameter space.
    for (int sim_id = 0; sim_id < n_sims; ++sim_id)
    {
        delete pars[sim_id];
        delete sim[sim_id];
        gsl_rng_free(rng[sim_id]);
    }
    delete[] rng;

    // De-allocate the default scenario parameters.
    // Note it is not our responsibility to deallocate the overriding theta parameters
    // as input of the constructor.
    for (int i = 0; i < n_scenarios; ++i)
    {
        delete default_theta[i];
    }
    delete[] default_theta;
    delete[] theta;
    
    // Free up the transformed parameters.
    delete[] vac_prop;

    // Free up the posterior sample IDs.
    delete[] post_idx;
}

// Scenario testing algorithm:
// Randomly sample from the posterior distribution, then randomly sample 
// from the distribution of each scenario parameter, then run the simulation,
// then store the result.
void Scenario::Run(Data* data, Posterior* post)
{
    // Set up a progress counter.
    int progress_idx = 0;

    // For each requested scenario sample:
    #pragma omp parallel for num_threads(n_sims)
    for (int i = 0; i < RVF_NSAMPLES; ++i)
    {
        // Define which simulation space to use.
        int sim_id = omp_get_thread_num();

        // For each scenario sample, run the scenario with different posterior samples.
        for (int j = 0; j < RVF_POST_NSAMPLES; ++j)
        {
            // Randomly sample which posterior sample to use.
            int post_sample_id = gsl_rng_uniform_int(rng[sim_id], post->n_samples);
            post_idx[i*RVF_POST_NSAMPLES + j] = post_sample_id;

            // Set the simulation parameters from this posterior sample.
            SetSimPars(post_sample_id, post, sim_id);
            
            // If a simulation without vaccination is requested, run one.
            if (RVF_NO_VAC)
            {
                pars[sim_id]->vac_rate = 0.0;
                sim[sim_id]->Simulate(data, pars[sim_id], 0, pars[sim_id]->n_steps);
                if (RVF_WRITE_FULL)
                {
                    #pragma omp critical
                    sim[sim_id]->WriteNoVac(i, j);
                }
                if (RVF_WRITE_SUMMARY)
                {
                    #pragma omp critical
                    sim[sim_id]->WriteNoVacSummary(i, j);
                }                
            }

            // Update the simulation parameters based on sampled scenario parameters.
            UpdateSimPars(i, sim_id);

            // Run the simulation for all time steps.
            sim[sim_id]->Simulate(data, pars[sim_id], 0, pars[sim_id]->n_steps);

            // Save the simulation result.
            // For some reason, this cannot be done in parallel, even if
            // the results are written to a string first before writing to file.
            if (RVF_WRITE_FULL)
            {
                #pragma omp critical
                sim[sim_id]->WriteOutput(i, j);
            }
            if (RVF_WRITE_SUMMARY)
            {
                #pragma omp critical
                sim[sim_id]->WriteSummary(i, j);
            }
        }

        // Update the progress counter.
        #pragma omp atomic update
        ++progress_idx;

        // Update the user on progress of the scenario tests.
        PrintProgress(progress_idx);
    }

    // Write the sampled scenario parameters to file.
    WriteScenarioPars();

    // Write the sampled posterior values to file.
    WritePosteriorSamples(post);
}

// Function to set the simulation parameters from the posterior.
void Scenario::SetSimPars(int sample_id, Posterior* post, int sim_id)
{
    // Set the parameters in the simulation.
    pars[sim_id]->p_immune_start[0] = post->p_immune_start_anj.chain[sample_id];
    pars[sim_id]->p_immune_start[1] = post->p_immune_start_gra.chain[sample_id];
    pars[sim_id]->p_immune_start[2] = post->p_immune_start_may.chain[sample_id];
    pars[sim_id]->p_immune_start[3] = post->p_immune_start_moh.chain[sample_id];
    pars[sim_id]->move[0][1] = post->move_anj_gra.chain[sample_id];
    pars[sim_id]->move[0][2] = post->move_anj_may.chain[sample_id];
    pars[sim_id]->move[0][3] = post->move_anj_moh.chain[sample_id];
    pars[sim_id]->move[1][0] = post->move_gra_anj.chain[sample_id];
    pars[sim_id]->move[1][3] = post->move_gra_moh.chain[sample_id];
    pars[sim_id]->move[3][0] = post->move_moh_anj.chain[sample_id];
    pars[sim_id]->move[3][1] = post->move_moh_gra.chain[sample_id];
    pars[sim_id]->import_start = post->import_start.chain[sample_id];
    pars[sim_id]->import_duration = post->import_duration.chain[sample_id];
    pars[sim_id]->import_size = post->import_size.chain[sample_id];

    // The model parameters used depend on configuration.
    switch(RVF_DIFF_NDVI)
    {
        case 0:     // Same NDVI per island.
            pars[sim_id]->ndvi_rate[0] = post->ndvi_rate.chain[sample_id];
            pars[sim_id]->ndvi_rate[1] = post->ndvi_rate.chain[sample_id];
            pars[sim_id]->ndvi_rate[2] = post->ndvi_rate.chain[sample_id];
            pars[sim_id]->ndvi_rate[3] = post->ndvi_rate.chain[sample_id];
            break;
        case 1:     // Different NDVI per island.
            pars[sim_id]->ndvi_rate[0] = post->ndvi_rate_anj.chain[sample_id];
            pars[sim_id]->ndvi_rate[1] = post->ndvi_rate_gra.chain[sample_id];
            pars[sim_id]->ndvi_rate[2] = post->ndvi_rate_may.chain[sample_id];
            pars[sim_id]->ndvi_rate[3] = post->ndvi_rate_moh.chain[sample_id];
            break;
        default:    // Default is same.
            pars[sim_id]->ndvi_rate[0] = post->ndvi_rate.chain[sample_id];
            pars[sim_id]->ndvi_rate[1] = post->ndvi_rate.chain[sample_id];
            pars[sim_id]->ndvi_rate[2] = post->ndvi_rate.chain[sample_id];
            pars[sim_id]->ndvi_rate[3] = post->ndvi_rate.chain[sample_id];
    }
    switch(RVF_DIFF_SCALE)
    {
        case 0:     // Same transmission scale per island.
            pars[sim_id]->trans_scale[0] = post->trans_scale.chain[sample_id];
            pars[sim_id]->trans_scale[1] = post->trans_scale.chain[sample_id];
            pars[sim_id]->trans_scale[2] = post->trans_scale.chain[sample_id];
            pars[sim_id]->trans_scale[3] = post->trans_scale.chain[sample_id];
            break;
        case 1:     // Different transmission scale per island.
            pars[sim_id]->trans_scale[0] = post->trans_scale_anj.chain[sample_id];
            pars[sim_id]->trans_scale[1] = post->trans_scale_gra.chain[sample_id];
            pars[sim_id]->trans_scale[2] = post->trans_scale_may.chain[sample_id];
            pars[sim_id]->trans_scale[3] = post->trans_scale_moh.chain[sample_id];
            break;
        default:    // Default is same.
            pars[sim_id]->trans_scale[0] = post->trans_scale.chain[sample_id];
            pars[sim_id]->trans_scale[1] = post->trans_scale.chain[sample_id];
            pars[sim_id]->trans_scale[2] = post->trans_scale.chain[sample_id];
            pars[sim_id]->trans_scale[3] = post->trans_scale.chain[sample_id];
    }
}

// Function to update the simulation parameters based on scenario tests.
void Scenario::UpdateSimPars(int sample_id, int sim_id)
{
    // Parameters to update depends on the scenarios requested.
    pars[sim_id]->vac_rate = vac_rate->chain[sample_id];
    pars[sim_id]->vac_protect_duration = vac_protect_duration->chain[sample_id];
    pars[sim_id]->vac_efficacy = vac_efficacy->chain[sample_id];
    pars[sim_id]->reg_import_scale = reg_import_scale->chain[sample_id];
    pars[sim_id]->vac_t_year_length = static_cast<int>(vac_t_year_length->chain[sample_id]); // Could be made more clean in theta.h
    pars[sim_id]->vac_t_year_offset = static_cast<int>(vac_t_year_offset->chain[sample_id]);
    pars[sim_id]->n_age_vac = static_cast<int>(n_age_vac->chain[sample_id]);
    pars[sim_id]->vac_t_freq = static_cast<int>(vac_t_freq->chain[sample_id]);
    pars[sim_id]->vac_identifiable = static_cast<int>(vac_identifiable->chain[sample_id]);

    // Set the distribution of vaccines. This distribution will depend on how 
    // many parameters are in the Dirchlet distribution.
    int k = vac_dist->GetNumDims(); // Number dimensions of the Dirichlet (or Empirical) distribution.
    if (k == 2) {
        // The scenario distribution determines the proportions for GC+Moh (by pop) and Anj+May (by pop).
        // Set the population weight of Anjouan between Anjouan and Mayotte.
        double anj_weight = pars[sim_id]->n_pop[0] / (pars[sim_id]->n_pop[0] + pars[sim_id]->n_pop[2]);
        double gra_weight = pars[sim_id]->n_pop[1] / (pars[sim_id]->n_pop[1] + pars[sim_id]->n_pop[3]);

        // Order of the islands is Anjouan, Grande Comore, Mayotte and Moheli.
        vac_prop[sample_id*n_island] = anj_weight*vac_dist->chain[sample_id*k];
        vac_prop[sample_id*n_island + 1] = gra_weight*vac_dist->chain[sample_id*k + 1];
        vac_prop[sample_id*n_island + 2] = (1.0 - anj_weight)*vac_dist->chain[sample_id*k];
        vac_prop[sample_id*n_island + 3] = (1.0 - gra_weight)*vac_dist->chain[sample_id*k + 1];

    } else if (k == 3) {
        // The scenario distribution determines the proportions for Grande Comore, Moheli and Anj+May (by pop).
        // Set the population weight of Anjouan between Anjouan and Mayotte.
        double anj_weight = pars[sim_id]->n_pop[0] / (pars[sim_id]->n_pop[0] + pars[sim_id]->n_pop[2]);

        // Order of the islands is Anjouan, Grande Comore, Mayotte and Moheli.
        vac_prop[sample_id*n_island] = anj_weight*vac_dist->chain[sample_id*k];
        vac_prop[sample_id*n_island + 1] = vac_dist->chain[sample_id*k + 1];
        vac_prop[sample_id*n_island + 2] = (1.0 - anj_weight)*vac_dist->chain[sample_id*k];
        vac_prop[sample_id*n_island + 3] = vac_dist->chain[sample_id*k + 2];

    } else if (k == 4) {
        // This is simple. The proportions at which to ditribute the vaccine are determined by the distribution itself.
        for (int i = 0; i < n_island; ++i){ vac_prop[sample_id*n_island + i] = vac_dist->chain[sample_id*k + i]; }
    } else if (k == 1) {
        // By default, k == 1, which means the default vaccine proportions should be used (proportional to pop size).
        for (int i = 0; i < n_island; ++i){ vac_prop[sample_id*n_island + i] = pars[sim_id]->n_pop[i] / pars[sim_id]->CalculateTotalPop(); }
    }

    // Copy across the vaccine proportions to the transformed scenario parameters.
    for (int i = 0; i < n_island; ++i){ pars[sim_id]->vac_prop[i] = vac_prop[sample_id*n_island + i]; }
}

// Function to print the progress of the scenario testing to console.
void Scenario::PrintProgress(int i)
{
    if (i % RVF_UPDATE_FREQ == 0)
    {
        // Note how many iterations have been completed.
        std::cout << "\rScenario replicates completed: " << i << " / " << RVF_NSAMPLES;
    }
}

// Function to write the sampled posterior parameters to file.
void Scenario::WritePosteriorSamples(Posterior* post)
{
    // Declare the file stream.
    std::fstream file;

    // Format the scenario_id into at least 3 digits.
    char scenario_id_buffer[4];
    sprintf(scenario_id_buffer, "%.3d", scenario_id);

    // Open a file for writing.
    file.open(static_cast<std::string>(RVF_ODIR) + "posterior_" + scenario_id_buffer + ".csv", std::fstream::out);

    // If we're at the start, open the file without appending.
    if (file.is_open())
    {
        // Output headers to the file!
        file << "SAMPLE_ID,REP_ID,PAR_NAME,PAR_VALUE,POST_ID";

        // For each scenario tested, save its ID, the scenario parameter names
        // and their values.
        for (int i = 0; i < RVF_NSAMPLES; ++i)
        {
            for (int j = 0; j < RVF_POST_NSAMPLES; ++j)
            {
                for (int k = 0; k < post->n_pars; ++k)
                {
                    file << "\n" << i << "," << j << "," << post->theta[k]->par_name;
                    file << "," << post->theta[k]->chain[post_idx[i*RVF_POST_NSAMPLES + j]];
                    file << "," << post_idx[i*RVF_POST_NSAMPLES + j];
                }
            }
        }
    }
    else
    {
        std::cout << "Error: there was a problem opening the file for writing";
        std::cout << " scenario output." << std::endl;
        exit(1);
    }

    // Close the file after writing has finished.
    file.close();    
}

// Function to write the scenario parameters to file.
void Scenario::WriteScenarioPars()
{
    // Declare the file stream.
    std::fstream file;

    // Format the scenario_id into at least 3 digits.
    char scenario_id_buffer[4];
    sprintf(scenario_id_buffer, "%.3d", scenario_id);

    // Open a file for writing.
    file.open(static_cast<std::string>(RVF_ODIR) + "scenario_pars_" + scenario_id_buffer + ".csv", std::fstream::out);

    // If we're at the start, open the file without appending.
    if (file.is_open())
    {
        // Output headers to the file!
        file << "SAMPLE_ID,PAR_NAME,INDEX,PAR_VALUE";

        // For each scenario tested, save its ID, the scenario parameter names
        // and their values.
        for (int i = 0; i < RVF_NSAMPLES; ++i)
        {
            // Scenario parameters.
            for (int j = 0; j < n_scenarios; ++j)
            {
                file << "\n" << theta[j]->OutputString(i);
            }

            // Transformed scenario parameters.
            for (int k = 0; k < n_island; ++k)
            {
                file << "\n" << i << ",vac_prop," << k << "," << vac_prop[i*n_island + k];
            }
        }
    }
    else
    {
        std::cout << "Error: there was a problem opening the file for writing";
        std::cout << " scenario output." << std::endl;
        exit(1);
    }

    // Close the file after writing has finished.
    file.close();    
}
