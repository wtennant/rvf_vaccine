// simulation.cpp: Definition of the function that simulates RVF.
#include <cstdio>
#include <chrono>               // Seed for random number generation.
#include <fstream>              // Writing data to file.
#include <iostream>             // Input-output to console.
#include <math.h>               // Mathematical operations.
#include "gsl/gsl_blas.h"
#include "omp.h"                // Parallelism.
#include "config.h"             // Location to output file.
#include "Simulation.h"         // Definition of the Simulation class.
#include <cstdio>

// Constructor to allocate memory for each infection state, age and time step.
Simulation::Simulation(Parameters* pars, int sens_id, int sim_id)
{
    // Store local copies of the number of age groups, islands and time steps.
    n_age = pars->n_age;
    n_age_move = pars->n_age_move;
    n_island = pars->n_island;
    n_steps = pars->n_steps;

    // Calculate the scaling factor for imports in each age group.
    imp_scaling = 0.0;
    for (int age = 0; age < n_age_move; ++age)
    {
        imp_scaling += pars->pop_structure[age]; 
    }

    // Allocate space for each time step...
    compartments = new double**[n_steps + 1];

    // For each time step...
    for (int t = 0; t <= n_steps; ++t)
    {
        // Allocate space for each island.
        compartments[t] = new double*[n_island];

        // For each island...
        for (int i = 0; i < n_island; ++i)
        {
            // Allocate space for each age and epidemiological group.
            compartments[t][i] = new double[n_age*n_epi];
        }
    }

    // Allocate time step to counters.
    doses = new double**[n_steps + 1];

    // For each time step, allocate space for each island.
    for (int t = 0; t <= n_steps; ++t)
    {
        doses[t] = new double*[n_island];
            
        // For each island, allocate age groups.
        for (int i = 0; i < n_island; ++i)
        {
            doses[t][i] = new double[n_age];

            for (int age = 0; age < n_age; ++age)
            {
                doses[t][i][age] = 0.0;
            }
        }
    }
    
    // Allocate space for each source island.
    compartments_move = new double**[n_island];

    // For each source island...
    for (int i = 0; i < n_island; ++i)
    {
        // Allocate space for each destination island.
        compartments_move[i] = new double*[n_island];

        // For each destination island...
        for (int j = 0; j < n_island; ++j)
        {
            // Allocate space for each age group that's allowed to move.
            compartments_move[i][j] = new double[n_age_move*n_epi];

            // For each age group...
            for (int age_epi = 0; age_epi < n_age_move*n_epi; ++age_epi)
            {
                compartments_move[i][j][age_epi] = 0.0;
            }
        }
    }

    // Keep track of the net movement numbers on each island.
    net_move = new double[n_island];

    // Create a counter to keep track of the number of deaths per island.
    deaths = new double[n_island];

    // Set up counter to store the number of individuals that are up for
    // moving on each island.
    total_move_N = new double[n_island];

    // Define space for a buffer used to write to file.
    file_buffer = new char[RVF_MAX_BUFFER];
    p_buffer = file_buffer;

    // Format the sensitivity ID into 3 digits.
    char sens_id_buffer[4];
    sprintf(sens_id_buffer, "%.3d", sens_id);

    // Write all compartments for all islands, ages and time points.
    if (RVF_WRITE_FULL)
    {
        // Open up a file stream for recording simulated data.
        
        file.open(static_cast<std::string>(RVF_ODIR) + "simulation_raw_" + sens_id_buffer + "_" + std::to_string(sim_id) + ".csv", std::fstream::out);

        // Check that the file is ready for writing.
        if (!file.is_open())
        {
            std::cout << "Error: there was a problem opening the file for writing";
            std::cout << " simulation output." << std::endl;
            perror("Error code: ");
            exit(1);
        }
        else
        {
            // Write headers to file.
            file << "SAMPLE_ID,REP_ID,TIME,ISLAND_ID,AGE,S_U,E_U,I_U,R_U,S_V1,E_V1,I_V1,R_V1,S_V2,E_V2,I_V2,R_V2,S_W,E_W,I_W,R_W,DOSES";
        }
    }

    // Write all compartments for all islands and time points.
    if (RVF_WRITE_SUMMARY)
    {
        // Open up a file stream for recording simulated data.
        file_summary.open(static_cast<std::string>(RVF_ODIR) + "simulation_summary_" + sens_id_buffer + "_" + std::to_string(sim_id) + ".csv", std::fstream::out);

        // Check that the file is ready for writing.
        if (!file_summary.is_open())
        {
            std::cout << "Error: there was a problem opening the summary file for writing";
            std::cout << " simulation output." << std::endl;
            perror("Error code: ");
            exit(1);
        }
        else
        {
            // Write headers to file.
            file_summary << "SAMPLE_ID,REP_ID,TIME_ID,MIN_TIME,MAX_TIME,ISLAND_ID,S_U,E_U,I_U,R_U,S_V1,E_V1,I_V1,R_V1,S_V2,E_V2,I_V2,R_V2,S_W,E_W,I_W,R_W,DOSES";
        }
    }

    // Write output of a simulation without any vaccintion.
    if (RVF_NO_VAC)
    {
        // Write the output as is.
        if (RVF_WRITE_FULL)
        {
            // Open up a file stream for recording simulated data.
            file_no_vac.open(static_cast<std::string>(RVF_ODIR) + "simulation_no_vac_raw_" + sens_id_buffer + "_" + std::to_string(sim_id) + ".csv", std::fstream::out);

            // Check that the file is ready for writing.
            if (!file_no_vac.is_open())
            {
                std::cout << "Error: there was a problem opening the no vaccination raw file for writing";
                std::cout << " simulation output." << std::endl;
                perror("Error code: ");
                exit(1);
            }
            else
            {
                // Write headers to file.
                file_no_vac << "SAMPLE_ID,REP_ID,TIME,ISLAND_ID,AGE,S_U,E_U,I_U,R_U";
            }
        }

        // For a summary of the output information.
        if (RVF_WRITE_SUMMARY)
        {
            // Open up a file stream for recording simulated data.
            file_no_vac_summary.open(static_cast<std::string>(RVF_ODIR) + "simulation_no_vac_summary_" + sens_id_buffer + "_" + std::to_string(sim_id) + ".csv", std::fstream::out);

            // Check that the file is ready for writing.
            if (!file_no_vac_summary.is_open())
            {
                std::cout << "Error: there was a problem opening the no vaccination summary file for writing";
                std::cout << " simulation output." << std::endl;
                perror("Error code: ");
                exit(1);
            }
            else
            {
                // Write headers to file.
                file_no_vac_summary << "SAMPLE_ID,REP_ID,TIME_ID,MIN_TIME,MAX_TIME,ISLAND_ID,S_U,E_U,I_U,R_U";
            }
        }
    }
}

// Destructor to de-allocate memory for each infection state, age and time step.
Simulation::~Simulation()
{
    // Free-up memory that was allocated for movement.
    for (int i = 0; i < n_island; ++i)
    {
        for (int j = 0; j < n_island; ++j)
        {
            delete[] compartments_move[i][j];
        }
        delete[] compartments_move[i];
    }
    delete[] compartments_move;

    // Free-up memory that was allocated for the compartments.
    for (int t = 0; t <= n_steps; ++t)
    {
        for (int i = 0; i < n_island; ++i)
        {
            delete[] compartments[t][i];
        }
        delete[] compartments[t];
    }
    delete[] compartments;

    // Free-up the memory that was allocated for each time step.
    for (int t = 0; t <= n_steps; ++t)
    {
        for (int i = 0; i < n_island; ++i)
        {
            delete[] doses[t][i];
        }
        delete[] doses[t];
    }
    delete[] doses;

    // Close the files used to record simulation output if they were open.
    if (file.is_open()){ file.close(); }
    if (file_summary.is_open()){ file_summary.close(); }
    if (file_no_vac.is_open()){ file_no_vac.close(); }
    if (file_no_vac_summary.is_open()){ file_no_vac_summary.close(); }
    delete[] file_buffer;

    // De-allocate memory for dounters.
    delete[] deaths;
    delete[] net_move;
    delete[] total_move_N;
}

// Function to simulate the RVF given model parameters and NDVI data.
void Simulation::Simulate(Data* data, Parameters* pars, int start_time, int end_time)
{
    // If at the first time point, initialise the compartments.
    if (start_time == 0)
    {
        // First zero-out all compartments since we don't know if this model is for
        // a vaccine strategy or not.
        for (int i = 0; i < n_island; ++i)
        {
            for (int age_epi = 0; age_epi < n_age*n_epi; ++age_epi)
            {
                compartments[0][i][age_epi] = 0.0;
            }
        }

        // For each island in the metapopulation...
        for (int i = 0; i < n_island; ++i)
        {
            // Sample the number of individuals that are immune at the start of the simulation.
            double init_immune = pars->p_immune_start[i]*pars->n_pop[i];

            // Sample the age distriubtion of the initial livestock populations from a multinomial distribution.
            for (int age = 0; age < n_age; ++age)
            {
                compartments[0][i][age*n_epi + SU] = pars->pop_structure[age]*(pars->n_pop[i] - 10.0 - init_immune);
                compartments[0][i][age*n_epi + EU] = pars->pop_structure[age]*5.0;
                compartments[0][i][age*n_epi + IU] = pars->pop_structure[age]*5.0;
                compartments[0][i][age*n_epi + RU] = pars->pop_structure[age]*init_immune;
            }
        }
    }

    // Proportion of vaccinated individuals at each time step for which
    // immunity has waned.
    double prob_vac_wane = 1.0 - exp(-1.0 / pars->vac_protect_duration);

    // Calculate the probability of the incubation and infectious period ending.
    // This saves the calculation of the exponential function an excessive
    // number of times.
    double prob_inc_end = 1.0;
    double prob_inf_end = 1.0;

    // Calculate the probability of the first and second phases of vaccine
    // induced immunity moving to the second phase and protected class respectively.
    double prob_vac_first_end = 1.0;
    double prob_vac_second_end = 1.0;

    // Run the model for the desired number of time steps.
    for (int t = start_time; t < end_time; ++t)
    {
        // Check if vaccination is possible at the current time step.
        bool is_post_vac_t_first = (t >= (pars->vac_t_first + pars->vac_t_year_offset));
        bool is_pre_vac_t_year_end = (((t - pars->vac_t_year_offset) % 48) < pars->vac_t_year_length);
        bool is_vac_year = (((t - pars->vac_t_first - pars->vac_t_year_offset) / static_cast<int>(48)) % pars->vac_t_freq) == 0;

        // For each island, run the infection and ageing process.
        for (int i = 0; i < n_island; ++i)
        {
            // Reset the total number of individuals counters.
            net_move[i] = 0.0;

            // Determine which individuals survive the time step using the survival matrix.
            for (int age = 0; age < n_age; ++age)
            {
                for (int epi = 0; epi < n_epi; ++epi)
                {
                    compartments[t + 1][i][age*n_epi + epi] = (1.0 - pars->mortality[age]) * compartments[t][i][age*n_epi + epi];
                }
            }

            // Calculate how many individuals are currently alive, and how many deaths there were.
            deaths[i] = cblas_dasum(n_epi*n_age, compartments[t][i], 1) - cblas_dasum(n_epi*n_age, compartments[t + 1][i], 1);

            // Calculate how many individuals are eligable for vaccination.
            // This depends on whether animals are identifiable or not (i.e. do we know if they have already been vaccinated).
            double total_vaccine_eligible = 0.0;
            if (pars->vac_identifiable > 0){
                // If we know whether or not individuals have been vaccinated before, then only
                // vaccinate those who are the appropriate age and are unvaccinated.
                for (int age = 0; age < pars->n_age_vac; ++age)
                {
                    total_vaccine_eligible += compartments[t + 1][i][age*n_epi + SU] + compartments[t + 1][i][age*n_epi + EU] + compartments[t + 1][i][age*n_epi + IU] + compartments[t + 1][i][age*n_epi + RU];
                }
            } else {
                // If we don't know, vaccine doses are assumed to be spread proportionally
                // between animals of the correct age group.
                total_vaccine_eligible = cblas_dasum(n_epi*pars->n_age_vac, compartments[t + 1][i], 1);
            }

            // Calculate the total population on the island.
            double N = cblas_dasum(n_epi*n_age, compartments[t + 1][i], 1);

            // Calculate the force of infection.
            double foi = exp(pars->ndvi_rate[i]*(data->ndvi.ndvi[i][t] - data->ndvi.min_local_ndvi[i]) + pars->trans_scale[i]);

            // Calculate the number of infectious individuals. By construction, there are 4*n_age
            // infectious categories to sum over. Each of them are separated by 4 doubles, start in position 2 of
            // the compartment array.
            double total_infectious = cblas_dasum(4*n_age, &compartments[t + 1][i][2], 4);

            // Then calculate the probability of becoming infected assuming a Poisson process.
            double prob_inf = 1.0 - exp(-foi*total_infectious / N);
            double prob_vac_inf = prob_inf*(1.0 - pars->vac_efficacy);

            // For each age group, calculate the probability of the unvaccinated age groups.
            for (int age = 0; age < pars->n_age; ++age)
            {
                // Check if vaccination is possible with the current age group.
                bool is_vac_age = (age < pars->n_age_vac);

                // Check all vaccination conditions...
                double prob_vac = 0.0;
                if (is_post_vac_t_first && is_pre_vac_t_year_end && is_vac_age && is_vac_year)
                {
                    // The actual number of doses will depend on the total number of individuals remaining in each age group.
                    // One cannot vaccinate more individuals than actually exist!
                    prob_vac = pars->vac_rate * pars->vac_prop[i] / total_vaccine_eligible / static_cast<double>(pars->vac_t_year_length);
                    prob_vac = ((prob_vac > 1.0) || isinf(prob_vac)) ? 1.0 : prob_vac; // Catch infinity too, although in C: Inf > 1.0 == true.
                }

                // Calculate the total number of individuals involved in all 40 possible transitions.
                // Infection process on unvaccinated individuals.
                double SU_to_EU = (1.0 - prob_vac)*prob_inf*compartments[t + 1][i][age*n_epi + SU];
                double EU_to_IU = (1.0 - prob_vac)*prob_inc_end*compartments[t + 1][i][age*n_epi + EU];
                double IU_to_RU = (1.0 - prob_vac)*prob_inf_end*compartments[t + 1][i][age*n_epi + IU];

                // Unvaccinated individuals that become vaccinated.
                double SU_to_SV1 = prob_vac*(1.0 - prob_inf)*compartments[t + 1][i][age*n_epi + SU];
                double EU_to_EV1 = prob_vac*(1.0 - prob_inc_end)*compartments[t + 1][i][age*n_epi + EU];
                double IU_to_IV1 = prob_vac*(1.0 - prob_inf_end)*compartments[t + 1][i][age*n_epi + IU];
                double RU_to_RV1 = prob_vac*compartments[t + 1][i][age*n_epi + RU];

                // Infection process of unvaccinated individuals that are also vaccinated.
                double SU_to_EV1 = prob_vac*prob_inf*compartments[t + 1][i][age*n_epi + SU];
                double EU_to_IV1 = prob_vac*prob_inc_end*compartments[t + 1][i][age*n_epi + EU];
                double IU_to_RV1 = prob_vac*prob_inf_end*compartments[t + 1][i][age*n_epi + IU];

                // Infection process of vaccinated individuals in phase 1 that remain in phase 1.
                double SV1_to_EV1 = (1.0 - prob_vac_first_end)*prob_inf*compartments[t + 1][i][age*n_epi + SV1];
                double EV1_to_IV1 = (1.0 - prob_vac_first_end)*prob_inc_end*compartments[t + 1][i][age*n_epi + EV1];
                double IV1_to_RV1 = (1.0 - prob_vac_first_end)*prob_inf_end*compartments[t + 1][i][age*n_epi + IV1];

                // Vaccinated individuals in phase 1 that move to phase 2.
                double SV1_to_SV2 = prob_vac_first_end*(1.0 - prob_inf)*compartments[t + 1][i][age*n_epi + SV1];
                double EV1_to_EV2 = prob_vac_first_end*(1.0 - prob_inc_end)*compartments[t + 1][i][age*n_epi + EV1];
                double IV1_to_IV2 = prob_vac_first_end*(1.0 - prob_inf_end)*compartments[t + 1][i][age*n_epi + IV1];
                double RV1_to_RV2 = prob_vac_first_end*compartments[t + 1][i][age*n_epi + RV1];

                // Infection process of vaccinated phase 1 that move to phase 2.
                double SV1_to_EV2 = prob_vac_first_end*prob_inf*compartments[t + 1][i][age*n_epi + SV1];
                double EV1_to_IV2 = prob_vac_first_end*prob_inc_end*compartments[t + 1][i][age*n_epi + EV1];
                double IV1_to_RV2 = prob_vac_first_end*prob_inf_end*compartments[t + 1][i][age*n_epi + IV1];

                // Infection process of vaccinated individuals in phase 2 that remain in phase 2.
                double SV2_to_EV2 = (1.0 - prob_vac_second_end)*prob_inf*compartments[t + 1][i][age*n_epi + SV2];
                double EV2_to_IV2 = (1.0 - prob_vac_second_end)*prob_inc_end*compartments[t + 1][i][age*n_epi + EV2];
                double IV2_to_RV2 = (1.0 - prob_vac_second_end)*prob_inf_end*compartments[t + 1][i][age*n_epi + IV2];

                // Vaccinated individuals in phase 2 that move to protected.
                double SV2_to_SW = prob_vac_second_end*(1.0 - prob_inf)*compartments[t + 1][i][age*n_epi + SV2];
                double EV2_to_EW = prob_vac_second_end*(1.0 - prob_inc_end)*compartments[t + 1][i][age*n_epi + EV2];
                double IV2_to_IW = prob_vac_second_end*(1.0 - prob_inf_end)*compartments[t + 1][i][age*n_epi + IV2];
                double RV2_to_RW = prob_vac_second_end*compartments[t + 1][i][age*n_epi + RV2];

                // Infection process of vaccinated phase 2 that move to protected.
                double SV2_to_EW = prob_vac_second_end*prob_inf*compartments[t + 1][i][age*n_epi + SV2];
                double EV2_to_IW = prob_vac_second_end*prob_inc_end*compartments[t + 1][i][age*n_epi + EV2];
                double IV2_to_RW = prob_vac_second_end*prob_inf_end*compartments[t + 1][i][age*n_epi + IV2];

                // Infection process of vaccinated individuals.
                double SW_to_EW = (1.0 - prob_vac_wane)*prob_vac_inf*compartments[t + 1][i][age*n_epi + SW];
                double EW_to_IW = (1.0 - prob_vac_wane)*prob_inc_end*compartments[t + 1][i][age*n_epi + EW];
                double IW_to_RW = (1.0 - prob_vac_wane)*prob_inf_end*compartments[t + 1][i][age*n_epi + IW];

                // Individuals whose vaccine wanes.
                double SW_to_SU = prob_vac_wane*(1.0 - prob_vac_inf)*compartments[t + 1][i][age*n_epi + SW];
                double EW_to_EU = prob_vac_wane*(1.0 - prob_inc_end)*compartments[t + 1][i][age*n_epi + EW];
                double IW_to_IU = prob_vac_wane*(1.0 - prob_inf_end)*compartments[t + 1][i][age*n_epi + IW];
                double RW_to_RU = prob_vac_wane*compartments[t + 1][i][age*n_epi + RW];

                // Individuals whose vaccine wanes and move through the infection process.
                double SW_to_EU = prob_vac_wane*prob_vac_inf*compartments[t + 1][i][age*n_epi + SW];
                double EW_to_IU = prob_vac_wane*prob_inc_end*compartments[t + 1][i][age*n_epi + EW];
                double IW_to_RU = prob_vac_wane*prob_inf_end*compartments[t + 1][i][age*n_epi + IW];

                // Update all compartments.
                compartments[t + 1][i][age*n_epi + SU] += SW_to_SU - SU_to_EU - SU_to_EV1 - SU_to_SV1;
                compartments[t + 1][i][age*n_epi + EU] += SW_to_EU + EW_to_EU + SU_to_EU - EU_to_EV1 - EU_to_IU - EU_to_IV1;
                compartments[t + 1][i][age*n_epi + IU] += EW_to_IU + IW_to_IU + EU_to_IU - IU_to_IV1 - IU_to_RU - IU_to_RV1;
                compartments[t + 1][i][age*n_epi + RU] += IW_to_RU + RW_to_RU + IU_to_RU - RU_to_RV1;
                compartments[t + 1][i][age*n_epi + SV1] += SU_to_SV1 - SV1_to_EV1 - SV1_to_EV2 - SV1_to_SV2;
                compartments[t + 1][i][age*n_epi + EV1] += SU_to_EV1 + EU_to_EV1 + SV1_to_EV1 - EV1_to_EV2 - EV1_to_IV1 - EV1_to_IV2;
                compartments[t + 1][i][age*n_epi + IV1] += EU_to_IV1 + IU_to_IV1 + EV1_to_IV1 - IV1_to_IV2 - IV1_to_RV1 - IV1_to_RV2;
                compartments[t + 1][i][age*n_epi + RV1] += IU_to_RV1 + RU_to_RV1 + IV1_to_RV1 - RV1_to_RV2;
                compartments[t + 1][i][age*n_epi + SV2] += SV1_to_SV2 - SV2_to_EV2 - SV2_to_EW - SV2_to_SW;
                compartments[t + 1][i][age*n_epi + EV2] += SV1_to_EV2 + EV1_to_EV2 + SV2_to_EV2 - EV2_to_EW - EV2_to_IV2 - EV2_to_IW;
                compartments[t + 1][i][age*n_epi + IV2] += EV1_to_IV2 + IV1_to_IV2 + EV2_to_IV2 - IV2_to_IW - IV2_to_RV2 - IV2_to_RW;
                compartments[t + 1][i][age*n_epi + RV2] += IV1_to_RV2 + RV1_to_RV2 + IV2_to_RV2 - RV2_to_RW;
                compartments[t + 1][i][age*n_epi + SW] += SV2_to_SW - SW_to_EW - SW_to_EU - SW_to_SU;
                compartments[t + 1][i][age*n_epi + EW] += SV2_to_EW + EV2_to_EW + SW_to_EW - EW_to_EU - EW_to_IW - EW_to_IU;
                compartments[t + 1][i][age*n_epi + IW] += EV2_to_IW + IV2_to_IW + EW_to_IW - IW_to_IU - IW_to_RW - IW_to_RU;
                compartments[t + 1][i][age*n_epi + RW] += IV2_to_RW + RV2_to_RW + IW_to_RW - RW_to_RU;

                // Record the number of doses adminsitered in time window [t, t + 1].
                doses[t][i][age] = SU_to_SV1 + SU_to_EV1 + EU_to_EV1 + EU_to_IV1 + IU_to_IV1 + IU_to_RV1 + RU_to_RV1;
            }
            
            // Transition individuals between age groups using CBLAS and the age transition matrix.
            for (int age = n_age - 1; age > 0; --age)
            {
                for (int epi = 0; epi < n_epi; ++epi)
                {
                    compartments[t + 1][i][age*n_epi + epi] += compartments[t + 1][i][(age - 1)*n_epi + epi] / 48.0;
                    compartments[t + 1][i][(age - 1)*n_epi + epi] *= (47.0 / 48.0);
                }
            }      
        }

        // Calculate the number of movements between islands.
        CalculateMove(pars, t+1);

        // Carry out the movements between each island for all age groups.
        for (int i = 0; i < n_island; ++i)
        {
            // We do not need to check i = j so split the destination island j into two for loops.
            for (int j = 0; j < i; ++j)
            {
                // Check if any movement is even permitted between islands i and j.
                if (pars->move[i][j] > 0.0)
                {
                    // Keep track of all movements from island i to island j.
                    double movements = 0;

                    // Move all animals between the islands.
                    // Do this for each compartment in the infection process.
                    for (int age_epi = 0; age_epi < n_age_move*n_epi; ++age_epi)
                    {
                        compartments[t + 1][j][age_epi] += compartments_move[i][j][age_epi];
                        compartments[t + 1][i][age_epi] -= compartments_move[i][j][age_epi];

                        // Add up all the movements from island i to island j.
                        movements += compartments_move[i][j][age_epi];
                    }

                    // Calculate the net movements between island i and j.
                    net_move[i] -= movements;
                    net_move[j] += movements;
                }
            }
            for (int j = i + 1; j < n_island; ++j)
            {
                // Check if any movement is even permitted between islands i and j.
                if (pars->move[i][j] > 0.0)
                {
                    // Keep track of all movements from island i to island j.
                    double movements = 0;

                    // Move all animals between the islands.
                    // Do this for each compartment in the infection process.
                    for (int age_epi = 0; age_epi < n_age_move*n_epi; ++age_epi)
                    {
                        compartments[t + 1][j][age_epi] += compartments_move[i][j][age_epi];
                        compartments[t + 1][i][age_epi] -= compartments_move[i][j][age_epi];

                        // Add up all the movements from island i to island j.
                        movements += compartments_move[i][j][age_epi];
                    }

                    // Calculate the net movements between island i and j.
                    net_move[i] -= movements;
                    net_move[j] += movements;
                }
            }
        }
        
        // Carry out movement imports if in the appropriate time window.
        if ((t >= pars->import_start) && (t < pars->import_start + pars->import_duration))
        {               
            // Import infectious animals into Grande Comore.
            for (int age = 0; age < n_age_move; ++age)
            {
                compartments[t + 1][1][age*n_epi + IU] += pars->import_size*pars->pop_structure[age] / imp_scaling;
                net_move[1] += pars->import_size*pars->pop_structure[age] / imp_scaling;
            }
        }

        // This corresponds to the original fitted 
        // Determine if (i) regular imports can begin, and (ii) regular imports are occurring.
        bool is_reg_import = t > pars->reg_import_t_start;
        bool is_reg_import_before_end = ((t - static_cast<int>(pars->import_start)) % pars->reg_import_freq) < pars->import_duration;
        if ((is_reg_import) && (is_reg_import_before_end))
        {               
            // Import infectious animals into Grande Comore.
            for (int age = 0; age < n_age_move; ++age)
            {
                compartments[t + 1][1][age*n_epi + IU] += pars->reg_import_scale*pars->import_size*pars->pop_structure[age] / imp_scaling;
                net_move[1] += pars->import_size*pars->pop_structure[age] / imp_scaling;
            }
        }

        // For every island...
        for (int i = 0; i < n_island; ++i)
        {
            // ...rebalance the total number of individuals on island i.
            compartments[t + 1][i][0] += (deaths[i] - net_move[i]);
        }
    }
}

// Function to calculate the number of movements between islands at
// a given time step.
void Simulation::CalculateMove(Parameters* pars,
                               const int t)
{
    // Calculate total number of individuals who can be moved on each island.
    for (int i = 0; i < n_island; ++i)
    {
        total_move_N[i] = cblas_dasum(n_age_move*n_epi, compartments[t][i], 1);
    }

    // Move individuals from island i...
    for (int i = 0; i < n_island; ++i)
    {
        // Calculate the number of animals to pass to other islands.
        // We do not need to check i = j so split the destination island j into two for loops.
        for (int j = 0; j < i; ++j)
        {
            // If there are movements between island i and j, then...
            if (pars->move[i][j] > 0.0)
            {
                // ...calculate the total number of movements between island i and island j.
                double movements = pars->move[i][j];

                // For each age group, transfer the infected individuals from island i to island j.
                for (int age_epi = 0; age_epi < n_age_move*n_epi; ++age_epi)
                {
                    compartments_move[i][j][age_epi]  = movements * compartments[t][i][age_epi] / total_move_N[i];
                }
            }
        }
        for (int j = i + 1; j < n_island; ++j)
        {
            // If there are movements between island i and j, then...
            if (pars->move[i][j] > 0.0)
            {
                // ...calculate the total number of movements between island i and island j.
                double movements = pars->move[i][j];

                // For each age group, transfer the infected individuals from island i to island j.
                for (int age_epi = 0; age_epi < n_age_move*n_epi; ++age_epi)
                {
                    compartments_move[i][j][age_epi]  = movements * compartments[t][i][age_epi] / total_move_N[i];
                }
            }
        }
    }
}

// Function which formats values and places them in the buffer. The buffer is flushed to the file and reset
// if the buffer would be overloaded by an update.
template<typename... Args>
void Simulation::write_buffer(std::fstream& my_file, int* buffer_left, const char* format, Args... args)
{
    // Define variables to keep track of whether the data has been written to the buffer,
    // and the space needed in the buffer to write the data.
    int next_buffer_size;
    bool is_written_in_buffer = false;
    
    // Until the contents of the file are written fully in the buffer...
    while (!is_written_in_buffer)
    {
        // Write the number of individuals in each infection compartment
        // to the file.
        next_buffer_size = snprintf(p_buffer, *buffer_left, format, args...);

        // Check that all characters (including the terminating character) were written to the buffer.
        if (next_buffer_size + 1 <= *buffer_left)
        {
            // If there was enough room to write to the buffer, then:
            // (1) update how much of the buffer is left for writing, excluding the terminating character.
            *buffer_left -= next_buffer_size;
            
            // (2) move the pointer to the buffer along the number of characters just written
            //     excluding the terminating character.
            p_buffer += next_buffer_size;

            // (3) mark that the file has been written - yay!
            is_written_in_buffer = true;
        }
        else
        {
            // If there was not enough space in the buffer, then:
            // (1) write the contents of the buffer to the file.
            my_file.write(file_buffer, p_buffer - file_buffer);
            my_file.flush();

            // (2) reset the pointer to the start of the file buffer.
            p_buffer = file_buffer;

            // Reset the size of the buffer remaining.
            *buffer_left = RVF_MAX_BUFFER;
        }
    }
}

// Function to write simulation output to file.
void Simulation::WriteOutput(int iter, int rep)
{
    // Counter for how much of the buffer has been used so far.
    int buffer_left = RVF_MAX_BUFFER;

    // For every time step...
    for (int t = 0; t <= n_steps; ++t)
    {
        // For every island...
        for (int i = 0; i < n_island; ++i)
        {
            // ...for every age group...
            for (int age = 0; age < n_age; ++age)
            {
                // Write the scenario id, posterior id, time step, island id and age of the
                // population to a buffer, and flush to a file if the buffer would be overloaded.
                write_buffer(file, &buffer_left, "\n%d,%d,%d,%d,%d", iter, rep, t, i, age);

                // ...for each compartment.
                for (int epi = 0; epi < n_epi; ++epi)
                {
                    // Get the level of precision required to record the compartments value to the desired precision.
                    // Always record to at least the requested number of decimal places or significant figures.
                    int precision = static_cast<int>(floor(log10(compartments[t][i][age*n_epi + epi]))) + 1;
                    precision = (precision <= 0) ? RVF_WRITE_MIN_PRECISION: precision + RVF_WRITE_MIN_PRECISION;

                    // Write the the number of individuals in the compartment, and flush to a file if the buffer would be overloaded.
                    write_buffer(file, &buffer_left, ",%.*g", precision, compartments[t][i][age*n_epi + epi]);
                }

                // Get the level of precision required to record the doses value to the desired precision.
                // Always record to at least the requested number of decimal places or significant figures.
                int precision = static_cast<int>(floor(log10(doses[t][i][age]))) + 1;
                precision = (precision <= 0) ? RVF_WRITE_MIN_PRECISION: precision + RVF_WRITE_MIN_PRECISION;

                // Write the the number of doses administered, and flush to a file if the buffer would be overloaded.
                write_buffer(file, &buffer_left, ",%.*g", precision, doses[t][i][age]);
            }
        }
    }

    // Send the contents to the filestream to the file.
    file.write(file_buffer, p_buffer - file_buffer);
    file.flush();

    // Reset the pointer to the buffer.
    p_buffer = file_buffer;
}

// Function to write a summary of simulation output to file.
void Simulation::WriteSummary(int iter, int rep)
{
    // Counter for how much of the buffer has been used so far.
    int buffer_left = RVF_MAX_BUFFER;

    // Define the number of grouped time-steps to summarise the results over.
    int n_groups = n_steps / RVF_SUMMARY_SIZE + 1;

    // For every grouped time step...
    for (int g = 0; g < n_groups; ++g)
    {
        // Define the minimum and maximum timesteps to group the results over.
        int min_t = g*RVF_SUMMARY_SIZE;
        int max_t = (g + 1)*RVF_SUMMARY_SIZE - 1;
        max_t = max_t > n_steps ? n_steps : max_t; 

        // ...for every island...
        for (int i = 0; i < n_island; ++i)
        {
            // Write the scenario id, posterior id, grouped time step and island id
            // a buffer, and flush to a file if the buffer would be overloaded.
            write_buffer(file_summary, &buffer_left, "\n%d,%d,%d,%d,%d,%d", iter, rep, g, min_t, max_t, i);

            // ...for each compartment.
            for (int epi = 0; epi < n_epi; ++epi)
            {   
                // Coutn the number of indiivduals in the compartment across
                // the desired group. Initialise a counter to one.
                double individuals = 0.0;

                // For each time point in the grouped time steps...
                for (int t = min_t; t <= max_t; ++t)
                {
                    // ...and for every age group, add up the number of individuals in each compartment.
                    for (int age = 0; age < n_age; ++age)
                    {
                        individuals += compartments[t][i][age*n_epi + epi];
                    }
                }

                // Get the level of precision required to record the compartments value to the desired precision.
                // Always record to at least the requested number of decimal places or significant figures.
                int precision = static_cast<int>(floor(log10(individuals))) + 1;
                precision = (precision <= 0) ? RVF_WRITE_MIN_PRECISION: precision + RVF_WRITE_MIN_PRECISION;

                // Write the the number of individuals in the group, and flush to a file if the buffer would be overloaded.
                write_buffer(file_summary, &buffer_left, ",%.*g", precision, individuals);
            }

            // Add up the number of doses for each group as well.
            double individuals = 0.0;
            for (int t = min_t; t <= max_t; ++t)
            {
                for (int age = 0; age < n_age; ++age)
                {
                    individuals += doses[t][i][age];
                }
            }

            // Get the level of precision required to record the compartments value to the desired precision.
            // Always record to at least the requested number of decimal places or significant figures.
            int precision = static_cast<int>(floor(log10(individuals))) + 1;
            precision = (precision <= 0) ? RVF_WRITE_MIN_PRECISION: precision + RVF_WRITE_MIN_PRECISION;

            // Write the the number of individuals in the group, and flush to a file if the buffer would be overloaded.
            write_buffer(file_summary, &buffer_left, ",%.*g", precision, individuals);
        }
    }

    // Send the contents to the filestream to the file.
    file_summary.write(file_buffer, p_buffer - file_buffer);
    file_summary.flush();

    // Reset the pointer to the buffer.
    p_buffer = file_buffer;
}

// Function to write simulation of the no-vaccination scenario to file.
void Simulation::WriteNoVac(int iter, int rep)
{
    // Counter for how much of the buffer has been used so far.
    int buffer_left = RVF_MAX_BUFFER;

    // For every time step...
    for (int t = 0; t <= n_steps; ++t)
    {
        // For every island...
        for (int i = 0; i < n_island; ++i)
        {
            // ...for every age group...
            for (int age = 0; age < n_age; ++age)
            {
                // Write the scenario id, posterior id, time step, island id and age of the
                // population to a buffer, and flush to a file if the buffer would be overloaded.
                write_buffer(file_no_vac, &buffer_left, "\n%d,%d,%d,%d,%d", iter, rep, t, i, age);
   
                for (int epi = 0; epi < n_no_vac_epi; ++epi)
                {
                    // Get the level of precision required to record the compartments value to the desired precision.
                    // Always record to at least the requested number of decimal places or significant figures.
                    int precision = static_cast<int>(floor(log10(compartments[t][i][age*n_epi + epi]))) + 1;
                    precision = (precision <= 0) ? RVF_WRITE_MIN_PRECISION: precision + RVF_WRITE_MIN_PRECISION;

                    // Write the the number of individuals in the compartment, and flush to a file if the buffer would be overloaded.
                    write_buffer(file_no_vac, &buffer_left, ",%.*g", precision, compartments[t][i][age*n_epi + epi]);
                }
            }
        }
    }

    // Send the remaining contents to the file.
    file_no_vac.write(file_buffer, p_buffer - file_buffer);
    file_no_vac.flush();

    // Reset the pointer to the start of the buffer.
    p_buffer = file_buffer;
}

// Function to write a summary of simulation output to file.
void Simulation::WriteNoVacSummary(int iter, int rep)
{
    // Counter for how much of the buffer has been used so far.
    int buffer_left = RVF_MAX_BUFFER;

    // Define the number of grouped time-steps to summarise the results over.
    int n_groups = n_steps / RVF_SUMMARY_SIZE + 1;

    // For every grouped time step...
    for (int g = 0; g < n_groups; ++g)
    {
        // Define the minimum and maximum timesteps to group the results over.
        int min_t = g*RVF_SUMMARY_SIZE;
        int max_t = (g + 1)*RVF_SUMMARY_SIZE - 1;
        max_t = max_t > n_steps ? n_steps : max_t; 

        // ...for every island...
        for (int i = 0; i < n_island; ++i)
        {
            // Write the scenario id, posterior id, grouped time step and island id
            // a buffer, and flush to a file if the buffer would be overloaded.
            write_buffer(file_no_vac_summary, &buffer_left, "\n%d,%d,%d,%d,%d,%d", iter, rep, g, min_t, max_t, i);

            // ...for each compartment.
            for (int epi = 0; epi < n_no_vac_epi; ++epi)
            {   
                // Coutn the number of indiivduals in the compartment across
                // the desired group. Initialise a counter to one.
                double individuals = 0.0;

                // For each time point in the grouped time steps...
                for (int t = min_t; t <= max_t; ++t)
                {
                    // ...and for every age group, add up the number of individuals in each compartment.
                    for (int age = 0; age < n_age; ++age)
                    {
                        individuals += compartments[t][i][age*n_epi + epi];
                    }
                }

                // Get the level of precision required to record the compartments value to the desired precision.
                // Always record to at least the requested number of decimal places or significant figures.
                int precision = static_cast<int>(floor(log10(individuals))) + 1;
                precision = (precision <= 0) ? RVF_WRITE_MIN_PRECISION: precision + RVF_WRITE_MIN_PRECISION;

                // Write the the number of individuals in the group, and flush to a file if the buffer would be overloaded.
                write_buffer(file_no_vac_summary, &buffer_left, ",%.*g", precision, individuals);
            }
        }
    }

    // Send the contents to the filestream to the file.
    file_no_vac_summary.write(file_buffer, p_buffer - file_buffer);
    file_no_vac_summary.flush();

    // Reset the pointer to the buffer.
    p_buffer = file_buffer;
}
