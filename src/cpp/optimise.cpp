// optimise.cpp: Defines the functions used to optimise a given summary function
// across a defined set of parameters (optimise.h) using optimal stochastic annealing (OSA).
#include <chrono>               // Time used as the random number generator seed.
#include <cmath>                // Square root function.
#include <format>               // Formatting output strings.
#include <fstream>              // Writing to file.
#include <iostream>             // Output information to console.
#include <iomanip>              // Set precision of file output.
#include "omp.h"                // Parallelisation API.
#include "gsl/gsl_cblas.h"      // CBLAS operations.
#include "gsl/gsl_linalg.h"     // Cholesky decompositions.
#include "gsl/gsl_randist.h"    // Distributions to evaluate likelihood.
#include "gsl/gsl_statistics.h" // Calculating the covariance.
#include "config.h"             // Configuration of file output.
#include "optimise.h"           // Definition of the Scenario class.

// Constructor of the Optimise class.
// Takes the input of the default parameter set.
Optimise::Optimise(Parameters *default_pars, Data *default_data, Posterior *in_posterior)
{
    // Define the number of independent optimisation runs to take.
    n_optim = RVF_OPTIM_RUNS;

    // Define the number of threads available for parallelisation.
    // Take one less than the maximum number of threads in order to not lock out the machine.
    n_threads = RVF_NTHREADS;

    // Define the number of particles per optimisation run.
    n_particles = RVF_OPTIM_NPARTICLES;

    // Define the number of samples of the cost function to take per particle.
    n_cost_samples = RVF_OPTIM_NCOST;

    // Define the number of steps per optimisation run.
    n_steps = RVF_OPTIM_MAX_STEPS;

    // Define the number of parameter distributions involed in optimisation.
    n_dist = 1; // The number of parameters independent of scenario choice.

    // Allocate space for the optimisation distributions.
    vac_dist = new Dirichlet **[n_optim];
    for (int optim_id = 0; optim_id < n_optim; ++optim_id)
    {
        vac_dist[optim_id] = new Dirichlet *[n_steps + 2]; // An extra step for proposals.
        for (int step_id = 0; step_id < n_steps + 2; ++step_id)
        {
            vac_dist[optim_id][step_id] = new Dirichlet("vac_dist", 4, 3);
        }
    }

    // Space to store all the scenario parameters.
    theta = new Theta ***[n_optim];
    for (int optim_id = 0; optim_id < n_optim; ++optim_id)
    {
        theta[optim_id] = new Theta **[n_steps + 2]; // An extra step for proposals.
        for (int step_id = 0; step_id < n_steps + 2; ++step_id)
        {
            theta[optim_id][step_id] = new Theta *[n_dist]; // List all distributions.
            theta[optim_id][step_id][0] = vac_dist[optim_id][step_id];
        }
    }

    // Initialise the number of free optimisation parameters to zero.
    n_pars = 0;

    // Calculate the number of free optimisation parameters.
    for (int dist_id = 0; dist_id < n_dist; ++dist_id)
    {
        n_pars += theta[0][0][dist_id]->GetNumFreeDims();
    }

    // Define the list of possible tempering factors.
    n_delta = RVF_OPTIM_NDELTA;
    delta_opt = new double[n_delta];
    proposed_weights = new double *[n_delta]; // A list of proposed weights for each new possible delta value.
    for (int delta_idx = 0; delta_idx < n_delta; ++delta_idx)
    {
        delta_opt[delta_idx] = exp(-20.0 + delta_idx * 20.0 / static_cast<double>(RVF_OPTIM_NDELTA - 1));
        proposed_weights[delta_idx] = new double[n_particles]; // There needs to be a new weight for each particle.
    }

    // Define space to store:
    cost = new double ***[n_optim];             // Samples of the cost function for a particular run, particle and step in the SMC.
    summary_cost = new double **[n_optim];      // Summary of the cost samples for a particular run, particle and step in the SMC.
    tempered_cost = new double **[n_optim];     // Tempered cost function which permits a smooth transition from the initial distribution to the cost distribution.
    weight = new double **[n_optim];            // Normalised weights of each particle at each step in the SMC run.
    post_idx = new int ***[n_optim];            // The indices associated with the posterior distribution used to generate stochasticity in the cost function.
    tempering_delta = new double *[n_optim];    // Tempering factor used for each optimisation run and step in SMC.
    tempering_delta_idx = new int *[n_optim];   // The index of the list of possible tempering factors that can be used in the SMC.
    max_step_reached = new int[n_optim];        // Maximum step per optimisation run reached in the tempering SMC algorithm.
    log_proposal_ratio = new double *[n_optim]; // Logged ratio of proposals during particle rejuvenation step.
    for (int optim_id = 0; optim_id < n_optim; ++optim_id)
    {
        // The length of each chain should be the total number of steps in the optimisation,
        // plus one for initial states, and another to use for proposals.
        cost[optim_id] = new double **[n_steps + 2];
        summary_cost[optim_id] = new double *[n_steps + 2];
        tempered_cost[optim_id] = new double *[n_steps + 2];
        weight[optim_id] = new double *[n_steps + 2];
        post_idx[optim_id] = new int **[n_steps + 2];

        // The tempering parameter does not depend on the number of particles.
        tempering_delta[optim_id] = new double[n_steps + 2];
        tempering_delta_idx[optim_id] = new int[n_steps + 2];
        tempering_delta[optim_id][0] = 0.0;
        tempering_delta_idx[optim_id][0] = -1;

        // Each proposal ratio per optimisation only needs space for each particle.
        log_proposal_ratio[optim_id] = new double[n_particles];

        // For each step...
        for (int step_id = 0; step_id < n_steps + 2; ++step_id)
        {
            // The number of particles in each step of the chain is n_particles.
            cost[optim_id][step_id] = new double *[n_particles];
            summary_cost[optim_id][step_id] = new double[n_particles];
            tempered_cost[optim_id][step_id] = new double[n_particles];
            weight[optim_id][step_id] = new double[n_particles];
            post_idx[optim_id][step_id] = new int *[n_particles];

            // For each particle in an optimisation run...
            for (int particle_id = 0; particle_id < n_particles; ++particle_id)
            {
                // The number of samples to take from the cost function.
                cost[optim_id][step_id][particle_id] = new double[n_cost_samples];
                post_idx[optim_id][step_id][particle_id] = new int[n_cost_samples];

                // The initial weight of each particle is 1 over the number of particles.
                weight[optim_id][step_id][particle_id] = 1.0 / static_cast<double>(n_particles);
            }
        }
    }

    // Define whether or not a no-vaccination scenario has already been costed.
    is_no_vac_costed = new bool[in_posterior->n_samples];
    no_vac_cost = new double[in_posterior->n_samples * default_pars->n_island];
    for (int post_id = 0; post_id < in_posterior->n_samples; ++post_id)
    {
        is_no_vac_costed[post_id] = false;
        for (int i = 0; i < default_pars->n_island; ++i)
        {
            no_vac_cost[post_id * default_pars->n_island + i] = 0.0;
        }
    }

    // Define space for storing vaccination costs per CPU thread.
    vac_cost = new double[n_threads * default_pars->n_island];

    // Declare space for the covariance matrix, zero vector and result vector
    // for each active CPU thread.
    mean_vec = new gsl_vector *[n_threads];
    result_vec = new gsl_vector *[n_threads];
    work_vec = new gsl_vector *[n_threads];
    prev_vec = new gsl_vector *[n_threads];
    zero_vec = new gsl_vector *[n_threads];
    covariance = new gsl_matrix *[n_threads];
    rw_covariance = new gsl_matrix *[n_threads];
    for (int thread_id = 0; thread_id < n_threads; ++thread_id)
    {
        mean_vec[thread_id] = gsl_vector_alloc(n_pars);
        result_vec[thread_id] = gsl_vector_alloc(n_pars);
        work_vec[thread_id] = gsl_vector_alloc(n_pars);
        prev_vec[thread_id] = gsl_vector_alloc(n_pars);
        zero_vec[thread_id] = gsl_vector_alloc(n_pars);
        gsl_vector_set_zero(zero_vec[thread_id]);
        covariance[thread_id] = gsl_matrix_alloc(n_pars, n_pars);
        rw_covariance[thread_id] = gsl_matrix_alloc(n_pars, n_pars);
    }

    // Set up the random number generator based on time.
    // Store one for each active CPU thread.
    rng = new gsl_rng *[n_threads];
    for (int thread_id = 0; thread_id < n_threads; ++thread_id)
    {
        rng[thread_id] = gsl_rng_alloc(gsl_rng_mt19937);
        gsl_rng_set(rng[thread_id], std::chrono::system_clock::now().time_since_epoch().count() + thread_id);
    }

    // Allocate space for the model parameters.
    // Use the set of default parameters to deep copy parameter values across.
    // Store one for each active CPU thread.
    pars = new Parameters *[n_threads];
    for (int thread_id = 0; thread_id < n_threads; ++thread_id)
    {
        pars[thread_id] = new Parameters(*default_pars);
    }

    // Allocate space for the NDVI data.
    // Use the set of default data to deep copy data across.
    // Store one for each active CPU thread.
    data = new Data *[n_threads];
    for (int thread_id = 0; thread_id < n_threads; ++thread_id)
    {
        data[thread_id] = new Data(*default_data);
    }

    // Allocate space for running the simulations.
    // Store one for each active CPU thread.
    sim = new Simulation *[n_threads];
    for (int thread_id = 0; thread_id < n_threads; ++thread_id)
    {
        sim[thread_id] = new Simulation(pars[thread_id], 0, thread_id);
    }

    // Link up the input posterior to an internal address.
    posterior = in_posterior;
}

// Deconstructor of the Optimise class.
// De-allocates all dynamically allocated arrays.
Optimise::~Optimise()
{
    // De-allocate variables associated with optimisation id,
    // step in the SMC and the particle in the SMC.
    for (int optim_id = 0; optim_id < n_optim; ++optim_id)
    {
        for (int step_id = 0; step_id < n_steps; ++step_id)
        {
            for (int particle_id = 0; particle_id < n_particles; ++particle_id)
            {
                delete[] cost[optim_id][step_id][particle_id];
                delete[] post_idx[optim_id][step_id][particle_id];
            }
            delete[] summary_cost[optim_id][step_id];
            delete[] tempered_cost[optim_id][step_id];
            delete[] weight[optim_id][step_id];
            delete[] post_idx[optim_id][step_id];
            delete[] cost[optim_id][step_id];
        }
        delete[] summary_cost[optim_id];
        delete[] tempered_cost[optim_id];
        delete[] weight[optim_id];
        delete[] post_idx[optim_id];
        delete[] cost[optim_id];
        delete[] tempering_delta[optim_id];
        delete[] tempering_delta_idx[optim_id];
        delete[] log_proposal_ratio[optim_id];
    }
    delete[] summary_cost;
    delete[] tempered_cost;
    delete[] weight;
    delete[] post_idx;
    delete[] cost;
    delete[] tempering_delta;
    delete[] tempering_delta_idx;
    delete[] max_step_reached;
    delete[] log_proposal_ratio;

    // De-allocate space associated with the tempering parameter.
    for (int delta_id = 0; delta_id < n_delta; ++delta_id)
    {
        delete[] proposed_weights[delta_id];
    }
    delete[] delta_opt;

    // De-allocate space for no-vaccination costs.
    delete[] vac_cost;
    delete[] no_vac_cost;
    delete[] is_no_vac_costed;

    // De-allocate space for variables associated per thread.
    for (int thread_id = 0; thread_id < n_threads; ++thread_id)
    {
        gsl_matrix_free(covariance[thread_id]);
        gsl_matrix_free(rw_covariance[thread_id]);
        gsl_vector_free(zero_vec[thread_id]);
        gsl_vector_free(result_vec[thread_id]);
        gsl_vector_free(prev_vec[thread_id]);
        gsl_vector_free(work_vec[thread_id]);
        gsl_vector_free(mean_vec[thread_id]);
        gsl_rng_free(rng[thread_id]);
    }
    delete[] covariance;
    delete[] rw_covariance;
    delete[] zero_vec;
    delete[] result_vec;
    delete[] prev_vec;
    delete[] work_vec;
    delete[] mean_vec;
    delete[] rng;

    // De-allocate the optimising parameter space.
    for (int i = 0; i < n_dist; ++i)
    {
        delete[] theta[i];
    }
    delete[] theta;

    // Some memory hasn't been specifically deallocated.
}

// Optimisation algorithm (simulated annealing):
// Begin with a randomly sampled point in parameter space.
// Propose a candidate parameter set and evaluate a cost metric for the candidate
// parameter set and previously accepted parameter set.
// Accept or reject that parameter set based on how much the cost metric has improved between
// the candidate parameter set and the previously accepted parameter set.
// Due to noise in the cost metric, multiple cost metrics will be required to make a decision.
// This will dependent on the standard deviation of the cost metric.
// The standard deviation of the cost metric is approximated periodically in the algorithm.
// Randomly sample from the posterior distribution, then randomly sample
// from the distribution of each scenario parameter, then run the simulation,
// then store the result.
void Optimise::Run(int series_id)
{
    // Output the information to console. This is thread-safe.
    std::cout << "Starting optimisation...";

// For each sample of the optimal parameter which maximises the cost function...
#pragma omp parallel for num_threads(n_threads) schedule(dynamic, 1)
    for (int optim_id = 0; optim_id < n_optim; ++optim_id)
    {
        // Get the thread that the optimisation is running on.
        int thread_id = omp_get_thread_num();

        // Create a local stringstream to store output information into.
        std::stringstream local_output_string;

        // Output the information to console. This is thread-safe.
        std::cout << local_output_string.str();

        // For each particle in the optimisation...
        for (int particle_id = 0; particle_id < n_particles; ++particle_id)
        {
            // Begin by randomly sampling a set of parameters for the distribution to optimise over.
            for (int dist_id = 0; dist_id < n_dist; ++dist_id)
            {
                // Sample from the distribution, and store the result in the first position of the chain.
                // Note that the sample is assumed to be valid (i.e. satisfy the support for the distribution).
                theta[optim_id][0][dist_id]->Sample(particle_id);
            }

            // Calculate the cost function for the initial proposed parameter set.
            // Take multiple samples from the cost function.
            for (int cost_id = 0; cost_id < n_cost_samples; ++cost_id)
            {
                cost[optim_id][0][particle_id][cost_id] = CalculateCost(optim_id, 0, particle_id, cost_id);
            }

            // Summarise the cost function across the multiple samples.
            // This may be the minimum, mean, or maximum of the cost samples.
            summary_cost[optim_id][0][particle_id] = SummariseCost(optim_id, 0, particle_id);

            // Calculate the tempered-cost value that is used in the SMC.
            // This will take an initial tempering parameter of zero.
            // (i.e. all cost is associated with the initial sampling distribution).
            tempered_cost[optim_id][0][particle_id] = TemperCost(optim_id, 0, particle_id, 0.0);
        }

        // Create information string that initialisation has completed.
        local_output_string.str(std::string());
        local_output_string << std::format("\nInitialisation of optimisation {:{}d} of {:{}d} completed (series {:02d})",
                                           optim_id + 1, (int)floor(log10(n_optim)) + 1, n_optim, (int)floor(log10(n_optim)) + 1, series_id);

        // Output the information to console. This is thread-safe.
        std::cout << local_output_string.str();

        // Keep track of the number of rounds of SMC that have been completed.
        int step_id = 0;

        // For each step in the optimisation algorithm, propose a new optimised parameter set, and
        // decide whether to accept or reject the proposal.
        while ((tempering_delta_idx[optim_id][step_id] < (n_delta - 1)) && (step_id < RVF_OPTIM_MAX_STEPS))
        {
            // Calculate the next best tempering parmeter to use. We seek an effective
            // sample size of at least 50% of the total number of particles.
            // Initialise the next best tempering parameter to use, and the associated effective sample size (ESS).
            tempering_delta_idx[optim_id][step_id + 1] = tempering_delta_idx[optim_id][step_id];
            double best_ess = n_particles;

            // Now try every tempering parameter greater than the current one.
            for (int delta_id = tempering_delta_idx[optim_id][step_id] + 1; delta_id < n_delta; ++delta_id)
            {
                // We want to calculate the ESS using unnormalised weights.
                // We will keep track of a sum of the weights (n_sum) and the sum of square weights.
                double weight_sum = 0.0;
                double weight_pow_2_sum = 0.0;

                // Sum across all particles at the current step in the SMC.
                for (int particle_id = 0; particle_id < n_particles; ++particle_id)
                {
                    // Calculate the tempered cost of the particle with the new proposed tempering factor.
                    double updated_tempered_cost = TemperCost(optim_id, step_id, particle_id, delta_opt[delta_id]);

                    // Calculate what the weights would be with the proposed tempering factor.
                    // This is the previous weight multiplied by the new tempered cost over the tempered cost with the old tempering value.
                    proposed_weights[delta_id][particle_id] = weight[optim_id][step_id][particle_id] * updated_tempered_cost / tempered_cost[optim_id][step_id][particle_id];

                    // Add these weights (and squared weights) to the appropriate sums.
                    weight_sum += proposed_weights[delta_id][particle_id];
                    weight_pow_2_sum += gsl_pow_2(proposed_weights[delta_id][particle_id]);
                }

                // Square the weighted sum, and then calculate the ESS.
                weight_sum = gsl_pow_2(weight_sum);
                double ess = weight_sum / weight_pow_2_sum;

                // If that ESS is larger than half the number of particles, and
                // is the smallest ESS found so far, then record it as the best found so far.
                if ((ess >= 0.5 * n_particles) && (ess < best_ess))
                {
                    best_ess = ess;
                    tempering_delta_idx[optim_id][step_id + 1] = delta_id;
                }
            }

            // Make a note of the tempering value.
            tempering_delta[optim_id][step_id + 1] = delta_opt[tempering_delta_idx[optim_id][step_id + 1]];

            // Using the tempering value which gives the desired ESS, copy over the weights.
            cblas_dcopy(n_particles, proposed_weights[tempering_delta_idx[optim_id][step_id + 1]], 1, weight[optim_id][step_id + 1], 1);

            // Take the absolute sum (or just the same as the weights should be positive) of all weights.
            // And then normalise the weights by this value.
            double normalising_const = cblas_dasum(n_particles, weight[optim_id][step_id + 1], 1);
            cblas_dscal(n_particles, normalising_const, weight[optim_id][step_id + 1], 1);

            // Shuffle all the particles based on their weight (duplicates of the same particle are allowed).
            ResampleParticles(optim_id, step_id);

            // Pass the set of particles through a Metropolis-Hastings kernel until
            // a desired cumulative acceptance rate is reached. Also make sure that
            // the chain doesn't continue indefinitely.
            int mh_chain_id = 0;
            double cumulative_acceptance_rate = 0.0;
            while ((cumulative_acceptance_rate < 2.0) && (mh_chain_id < RVF_OPTIM_MAX_MH))
            {
                // Keep track of whether jittered particles are valid (within the domain of the optimal)
                // parameters or not.
                double particle_validity = 0.0;

                // Begin by proposing a new parameter set.
                // This has a small chance to randomly explore the entire parameter space.
                // Otherwise, it uses previously accepted values to suggest a potentially better
                // set of particles. Note that the the candidate parameter may be outside the
                // validity of the underlying parameter distribution.
                // Proposed particles are stored in step_id + 2.
                JitterParticles(optim_id, step_id + 1);

                // For each particle...
                for (int particle_id = 0; particle_id < n_particles; ++particle_id)
                {
                    // Check if the particle is valid.
                    // (i.e. the proposed particle is within the support of the distribution).
                    bool is_valid = IsParticleValid(optim_id, step_id + 2, particle_id);

                    // Define a boolean on whether the parameter is accepted or not.
                    bool is_accepted = false;

                    // If the proposal is valid, then...
                    if (is_valid)
                    {
                        // Increment the counter for whether a jittered particle was valid or not.
                        particle_validity += 1.0 / static_cast<double>(n_particles);

                        // Calculate the cost function for the initial proposed parameter set.
                        // Take multiple samples from the cost function.
                        for (int cost_id = 0; cost_id < n_cost_samples; ++cost_id)
                        {
                            cost[optim_id][step_id + 2][particle_id][cost_id] = CalculateCost(optim_id, step_id + 2, particle_id, cost_id);
                        }

                        // Summarise the cost function across the multiple samples.
                        // This may be the minimum, mean, or maximum of the cost samples.
                        summary_cost[optim_id][step_id + 2][particle_id] = SummariseCost(optim_id, step_id + 2, particle_id);

                        // Calculate the tempered-cost value that is used in the SMC.
                        // This will take an initial tempering parameter of zero.
                        // (i.e. all cost is associated with the initial sampling distribution).
                        tempered_cost[optim_id][step_id + 2][particle_id] = TemperCost(optim_id, step_id + 2, particle_id, tempering_delta[optim_id][step_id + 1]);

                        // Calculate the acceptance probability.
                        // As the proposal is symmetric, we do not need to factor this in.
                        double acceptance_prob = tempered_cost[optim_id][step_id + 2][particle_id] / tempered_cost[optim_id][step_id + 1][particle_id];

                        // Factor in the ratio of proposals to maintain detailed balance:
                        // that is, include q(theta | theta') / q(theta' | theta).
                        acceptance_prob *= exp(log_proposal_ratio[optim_id][particle_id]);

                        // Determine whether the sample is accepted or not.
                        is_accepted = gsl_ran_flat(rng[thread_id], 0.0, 1.0) < acceptance_prob;
                    }

                    // If the sample is rejected, only update the acceptance rate.
                    // If the sample is accepted, copy across all the proposed values
                    // (costs and parameters) to the current values.
                    if (is_accepted)
                    {
                        AcceptCandidate(optim_id, step_id + 1, particle_id);
                        cumulative_acceptance_rate += 1.0 / static_cast<double>(n_particles);
                    }
                }

                // Output information on the particle validity and cumulative acceptance rate at the current step.
                local_output_string.str(std::string());
                local_output_string << std::format("\nValidity and cumulative acceptance rate for MH-kernel in optimisation {:{}d}",
                                                   optim_id + 1, (int)floor(log10(n_optim)) + 1);
                local_output_string << std::format(", SMC-step {:{}d}, MH-step {:{}d}: {:3.0f}% {:3.0f}%",
                                                   step_id + 1, (int)floor(log10(RVF_OPTIM_MAX_STEPS)) + 1,
                                                   mh_chain_id + 1, (int)floor(log10(RVF_OPTIM_MAX_MH)) + 1,
                                                   particle_validity*100, cumulative_acceptance_rate*100);

                // Output the information to console. This is thread-safe.
                std::cout << local_output_string.str();

                // Move onto the next step of the Metropolis-Hastings kernel.
                ++mh_chain_id;
            }

            // If the limit of the Metropolis Hastings kernel is reached, then it should be flagged up
            // as this means the acceptance rate was poor.
            if (mh_chain_id == RVF_OPTIM_MAX_MH)
            {
                std::cout << "Metropolis-Hastings jittering kernel reached it's maximum number of steps.";
                std::cout << "\nConsider increasing the number of allowed steps in the kernel.";
                std::cout << "\nCumulative acceptance rate: " << cumulative_acceptance_rate << std::endl;
            }

            // Move onto the next step of the SMC algorithm.
            ++step_id;
        }

        // If the limit of steps has been reached...
        if (step_id == RVF_OPTIM_MAX_STEPS)
        {
            std::cout << "Tempering density Sequential Monte Carlo reached it's maximum number of steps.";
            std::cout << "\nConsider increasing the number of allowed steps in the SMC.";
        }

        // Record the maximum step that was reached in the current optimisation run.
        max_step_reached[optim_id] = step_id;
    }

    // Write the sampled optimal parameter chains to file.
    WriteOptimPars(series_id);

    // Write the sampled posterior values to file.
    // WritePostPars(series_id);

    // Write the cost values to file.
    WriteCost(series_id);
}

// Function to calculate the cost of a given parameter set.
// Inputs: ID of the optimisation run, step in the optimisation chain, ID of the particle, sample of the cost function.
// Output: Write total cost to the cost vector.
double Optimise::CalculateCost(int optim_id, int step_id, int particle_id, int sample_id)
{
    // Get the thread that the function is operating on.
    int thread_id = omp_get_thread_num();

    // Randomly sample which posterior sample to use.
    int post_sample_id = gsl_rng_uniform_int(rng[thread_id], posterior->n_samples);
    post_idx[optim_id][step_id][particle_id][sample_id] = post_sample_id;

    // Set the simulation parameters from the posterior sample.
    SetSimPostPars(thread_id, post_sample_id);

    if (!is_no_vac_costed[post_sample_id])
    {
        // Run a simulation without vaccination first.
        double vac_rate = pars[thread_id]->vac_rate;
        pars[thread_id]->vac_rate = 0.0;
        sim[thread_id]->Simulate(data[thread_id], pars[thread_id], 0, pars[thread_id]->n_steps);

        // Summarise the simulation with no vaccination.
        for (int i = 0; i < pars[thread_id]->n_island; ++i)
        {
            double no_vac_cost_sum = 0.0;
            for (int t = pars[thread_id]->vac_t_first + 1; t <= pars[thread_id]->n_steps; ++t)
            {

                no_vac_cost_sum += cblas_dasum(pars[thread_id]->n_age, sim[thread_id]->compartments[t][i] + 2, sim[thread_id]->n_epi);
            }

// Store the cost for that specific island.
#pragma omp critical
            no_vac_cost[post_sample_id * pars[thread_id]->n_island + i] = no_vac_cost_sum;
        }

        // Reset the vaccination rate.
        pars[thread_id]->vac_rate = vac_rate;

// Note that this scenario has now been costed.
#pragma omp critical
        is_no_vac_costed[post_sample_id] = true;
    }

    // Update the optimisation parameters in the simulation.
    SetSimOptimPars(thread_id, optim_id, step_id, particle_id);

    // Run the simulation for all time steps.
    sim[thread_id]->Simulate(data[thread_id], pars[thread_id], 0, pars[thread_id]->n_steps);

    // Summarise the simulation with no vaccination.
    for (int i = 0; i < pars[thread_id]->n_island; ++i)
    {
        double vac_cost_sum = 0.0;
        for (int t = pars[thread_id]->vac_t_first + 1; t <= pars[thread_id]->n_steps; ++t)
        {
            vac_cost_sum += cblas_dasum(4 * pars[thread_id]->n_age, sim[thread_id]->compartments[t][i] + 2, 4);
        }
        vac_cost[thread_id * pars[thread_id]->n_island + i] = vac_cost_sum;
    }

    // There are two choices for optimisation:
    // One is maximising the proportion of cases averted across all islands (RVF_OPTIM_FUNC == 0),
    // and the other is maximising the proportion of cases averted in the worst affected island (RVF_OPTIM_FUNC == 1).
    // Default cases are the first optimisation.
    double total_cost = 0.0;
    switch (RVF_OPTIM_FUNC)
    {
    case 0:
    {
        // Initialise the variables for the total number of infections across all islands under both no-vaccination and vaccination scenarios.
        double no_vac_total_cost = 0.0;
        double vac_total_cost = 0.0;

        // Add up all the infections for the no vaccination and vaccination scenarios.
        for (int i = 0; i < pars[thread_id]->n_island; ++i)
        {
            no_vac_total_cost += no_vac_cost[post_sample_id * pars[thread_id]->n_island + i];
            vac_total_cost += vac_cost[thread_id * pars[thread_id]->n_island + i];
        }

        // Calculate the proportion of cases averted with the vaccination scenario.
        total_cost = (no_vac_total_cost - vac_total_cost) / no_vac_total_cost;
        break;
    }
    case 1:
    {
        // Define the initial cost total cost for this scenario as 100%.
        total_cost = 1.0; // This means it will get caught when trying to find the island with the minimum number of cases averted.

        // Initialise a counter for the total cost for a specific island.
        double total_cost_island = 0.0;
        for (int i = 0; i < pars[thread_id]->n_island; ++i)
        {
            // Calculate the proportion of cases averted for island i.
            total_cost_island = (no_vac_cost[post_sample_id * pars[thread_id]->n_island + i] - vac_cost[thread_id * pars[thread_id]->n_island + i]);
            total_cost_island = total_cost_island / no_vac_cost[post_sample_id * pars[thread_id]->n_island + i];

            // If this new cost is smalled than the previously found smallest cost, update the final total cost for the entire scenario.
            if (total_cost > total_cost_island)
            {
                total_cost = total_cost_island;
            }
        }
        break;
    }
    default:
    {
        // Initialise the variables for the total number of infections across all islands under both no-vaccination and vaccination scenarios.
        double no_vac_total_cost = 0.0;
        double vac_total_cost = 0.0;

        // Add up all the infections for the no vaccination and vaccination scenarios.
        for (int i = 0; i < pars[thread_id]->n_island; ++i)
        {
            no_vac_total_cost += no_vac_cost[post_sample_id * pars[thread_id]->n_island + i];
            vac_total_cost += vac_cost[thread_id * pars[thread_id]->n_island + i];
        }

        // Calculate the proportion of cases averted with the vaccination scenario.
        total_cost = (no_vac_total_cost - vac_total_cost) / no_vac_total_cost;
        break;
    }
    }

    // Return the cost of the vaccination campaign.
    return total_cost;
}

// Function to summarise multiple samples taken from the cost function.
double Optimise::SummariseCost(int optim_id, int step_id, int particle_id)
{
    // The summary function to use will depend on the pre-compiled option in config.h.
    switch (RVF_OPTIM_SUMMARY)
    {
    case 0: // Take the mean of all cost samples.
        return gsl_stats_mean(cost[optim_id][step_id][particle_id], 1, n_cost_samples);
    case 1: // Take the minimum of all cost samples.
        return gsl_stats_min(cost[optim_id][step_id][particle_id], 1, n_cost_samples);
    case 2: // Take the maximum of all cost samples.
        return gsl_stats_max(cost[optim_id][step_id][particle_id], 1, n_cost_samples);
    default: // The default option should be the mean of all samples.
        return gsl_stats_mean(cost[optim_id][step_id][particle_id], 1, n_cost_samples);
    }
}

// Function which tempers the cost of a given strategy using a pre-set
// tempering parameter delta. The tempered cost is:
// f(theta)^delta * i(theta)^(1 - delta)
// where f is the pre-tempered cost and i is the initial distribution of theta.
double Optimise::TemperCost(int optim_id, int step_id, int particle_id, double delta)
{
    // Initialise a variable containing the final tempered cost of the particle.
    double tempered_final_cost = 1.0;

    // First add on the tempered intial cost.
    for (int dist_id = 0; dist_id < n_dist; ++dist_id)
    {
        tempered_final_cost *= theta[optim_id][step_id][dist_id]->Pdf(particle_id);
    }

    // Temper the initial cost.
    tempered_final_cost = pow(tempered_final_cost, 1.0 - delta);

    // Now temper the strategy cost.
    tempered_final_cost *= pow(exp(RVF_OPTIM_POWER * summary_cost[optim_id][step_id][particle_id]), delta);

    // Return this tempered overall cost.
    return tempered_final_cost;
}

// Function to resample a set of particles based on their current weight
// in the optimisation run. The corresponding weights of the particles
// must already be computed.
void Optimise::ResampleParticles(int optim_id, int step_id)
{
    // Re-allocate the weights in the discrete random number generator.
    // This random number generator, will generate a value between 0 and K-1.
    // The first argument is the number of unique values (particles), and
    // the second argument is the point to the vector of weights.
    gsl_ran_discrete_t *discrete_dist = gsl_ran_discrete_preproc(n_particles, weight[optim_id][step_id + 1]);

    // Get the ID of the active working thread.
    int thread_id = omp_get_thread_num();

    // For every single particle, sample a corresponding particle ID, and place it
    // into the next step in the algorithm.
    for (int particle_id = 0; particle_id < n_particles; ++particle_id)
    {
        // Use the discrete sampler to get a particle ID to place in the next set of particle IDs.
        int sampled_particle_id = gsl_ran_discrete(rng[thread_id], discrete_dist);

        // For the sampled parameters, go through each distribution we are optimising over.
        for (int dist_id = 0; dist_id < n_dist; ++dist_id)
        {
            // Get the dimensions of the distribution with ID dist_id.
            int n_dim = theta[optim_id][step_id][dist_id]->GetNumDims();

            // Copy across the information for each dimension...
            for (int dim_id = 0; dim_id < n_dim; ++dim_id)
            {
                // The chains contain all the particles, where a single particles (multiple, if applicable) parameters
                // from the distribution are next to each other in memory.
                theta[optim_id][step_id + 1][dist_id]->chain[particle_id * n_dim + dim_id] = theta[optim_id][step_id][dist_id]->chain[sampled_particle_id * n_dim + dim_id];
            }
        }

        // Use the sampled particle ID to create the next set of tempered costs, sampled parameters,
        // summary costs, and overall costs.
        summary_cost[optim_id][step_id + 1][particle_id] = summary_cost[optim_id][step_id][sampled_particle_id];
        tempered_cost[optim_id][step_id + 1][particle_id] = TemperCost(optim_id, step_id + 1, particle_id, tempering_delta[optim_id][step_id + 1]);
        cblas_dcopy(n_cost_samples, cost[optim_id][step_id][sampled_particle_id], 1, cost[optim_id][step_id + 1][particle_id], 1);
    }

    // De-allocate the discrete distribution.
    gsl_ran_discrete_free(discrete_dist);
}

// Function used when accepting a candidate parameter set.
void Optimise::AcceptCandidate(int optim_id, int step_id, int particle_id)
{
    // If a candidate parameter set is accepted, then replace it with
    // the proposed parameter set (one ahead in step).
    for (int dist_id = 0; dist_id < n_dist; ++dist_id)
    {
        // Get the dimension of the optimisation parameter distribution.
        int n_dim = theta[optim_id][step_id][dist_id]->GetNumDims();

        // Set each sample to the previously accepted value.
        for (int dim_id = 0; dim_id < n_dim; ++dim_id)
        {
            theta[optim_id][step_id][dist_id]->chain[particle_id * n_dim + dim_id] = theta[optim_id][step_id + 1][dist_id]->chain[particle_id * n_dim + dim_id];
        }
    }

    // Also copy across the associated cost values.
    tempered_cost[optim_id][step_id][particle_id] = tempered_cost[optim_id][step_id + 1][particle_id];
    summary_cost[optim_id][step_id][particle_id] = summary_cost[optim_id][step_id + 1][particle_id];
    for (int cost_id = 0; cost_id < n_cost_samples; ++cost_id)
    {
        cost[optim_id][step_id][particle_id][cost_id] = cost[optim_id][step_id + 1][particle_id][cost_id];
    }
}

// Function which proposes a new candidate parameter set.
// The only inputs are the ID of the optimisation chain,
// and the step the algorithm is in the optimisation chain
// (i.e. the parameter idx to be updated).
void Optimise::JitterParticles(int optim_id, int step_id)
{
    // Get the thread that the function is operating on.
    int thread_id = omp_get_thread_num();

    // First construct the covariance matrix of the current set of particles.
    // Keep track of the first parameter ID in the covariance matrix.
    int par1_id = 0;

    // For every distribution we are proposing over...
    for (int dist1_id = 0; dist1_id < n_dist; ++dist1_id)
    {
        // Get the number of dimensions in the first distribution.
        // We will only propose values for the free dimensions.
        int n_free_dim1 = theta[optim_id][step_id][dist1_id]->GetNumFreeDims();
        int n_dim1 = theta[optim_id][step_id][dist1_id]->GetNumDims();

        // For free parameter in the first distribution...
        for (int dim1_id = 0; dim1_id < n_free_dim1; ++dim1_id)
        {
            // Keep track of the second parameter ID in the covariance matrix.
            int par2_id = 0;

            // Calculate the mean of the current set of accepted parameters.
            double mean_val = gsl_stats_mean(&(theta[optim_id][step_id][dist1_id]->chain[dim1_id]), n_dim1, n_particles);
            gsl_vector_set(mean_vec[thread_id], par1_id, mean_val);

            // For every remaining distribution that we are optimising over...
            for (int dist2_id = 0; dist2_id < n_dist; ++dist2_id)
            {
                // Get the number of dimensions in the second distribution.
                int n_free_dim2 = theta[optim_id][step_id][dist2_id]->GetNumFreeDims();
                int n_dim2 = theta[optim_id][step_id][dist2_id]->GetNumDims();

                // For every free parameter in the second distribution...
                for (int dim2_id = 0; dim2_id < n_free_dim2; ++dim2_id)
                {
                    // Calculate the covariance between the first and second parameter.
                    // The parameters in the chain are separated by n_dim1 and n_dim2 values.
                    // E.g. if there are 4 dimensions to the parameter, then each value of the same parameter is
                    // repeated every 4 spaces in the chain.
                    double cov_val = gsl_stats_covariance(&(theta[optim_id][step_id][dist1_id]->chain[dim1_id]), n_dim1,
                                                          &(theta[optim_id][step_id][dist2_id]->chain[dim2_id]), n_dim2,
                                                          n_particles);

                    // Set the covariance matrix ready for sampling.
                    gsl_matrix_set(covariance[thread_id], par1_id, par2_id, cov_val);
                    gsl_matrix_set(rw_covariance[thread_id], par1_id, par2_id, cov_val);

                    // Move onto the next parameter ID.
                    ++par2_id;
                }
            }

            // Move onto the next first parameter ID.
            ++par1_id;
        }
    }

    // Scale the covariance matrix for the random-walk by our favourite  number of parameters we are trying to estimate.
    gsl_matrix_scale(rw_covariance[thread_id], 1.0 / static_cast<double>(n_pars));

    // Decompose the covariance matrix using Cholesky decomposition.
    gsl_linalg_cholesky_decomp1(rw_covariance[thread_id]);
    gsl_linalg_cholesky_decomp1(covariance[thread_id]);

    // Take the current particle, and either jitter it according to the covariance matrix
    // of the entire particle set, or propose a brand new particle sampled from the original underlying distribution.
    for (int particle_id = 0; particle_id < n_particles; ++particle_id)
    {
        // Set the log-ratio of proposed parameter to current parameter to zero.
        log_proposal_ratio[optim_id][particle_id] = 0.0;

        // Each particle has a small chance of being proposed from the original underlying distribution.
        // Also have a small chance of sampling from the multivariate gaussian as an independent sampler.
        double ran_unif = gsl_ran_flat(rng[thread_id], 0.0, 1.0);
        if (ran_unif < 1.0 / 3.0) // Sampling from the initialisation distribution.
        {
            for (int dist_id = 0; dist_id < n_dist; ++dist_id)
            {
                // Sample from the initialisation distribution.
                theta[optim_id][step_id + 1][dist_id]->Sample(particle_id);

                // Calculate the ratio of proposals---old particle - new particle.
                log_proposal_ratio[optim_id][particle_id] -= log(theta[optim_id][step_id + 1][dist_id]->Pdf(particle_id));
                log_proposal_ratio[optim_id][particle_id] += log(theta[optim_id][step_id][dist_id]->Pdf(particle_id));
            }
        }
        else if (ran_unif < 2.0 / 3.0) // Independence sampler (i.e. does not depend on current particle position).
        {
            // Otherwise create an independent sample using a multivariate gaussian.
            gsl_ran_multivariate_gaussian(rng[thread_id], mean_vec[thread_id], covariance[thread_id], result_vec[thread_id]);

            // Retrieve the sample and ammend to the next proposed sample.
            int par_id = 0;
            for (int dist_id = 0; dist_id < n_dist; ++dist_id)
            {
                // Get the number of dimensions in the first distribution.
                int n_free_dim1 = theta[optim_id][step_id][dist_id]->GetNumFreeDims();
                int n_dim1 = theta[optim_id][step_id][dist_id]->GetNumDims();
                for (int dim1_id = 0; dim1_id < n_free_dim1; ++dim1_id)
                {
                    theta[optim_id][step_id + 1][dist_id]->chain[particle_id * n_dim1 + dim1_id] = gsl_vector_get(result_vec[thread_id], par_id);

                    // At the same time retrieve the previously accepted parameter.
                    gsl_vector_set(prev_vec[thread_id], par_id, theta[optim_id][step_id][dist_id]->chain[particle_id * n_dim1 + dim1_id]);
                    ++par_id;
                }

                // Update the ratio of proposals.
                double lpdf;
                gsl_ran_multivariate_gaussian_log_pdf(result_vec[thread_id], mean_vec[thread_id], covariance[thread_id], &lpdf, work_vec[thread_id]);
                log_proposal_ratio[optim_id][particle_id] -= lpdf;
                gsl_ran_multivariate_gaussian_log_pdf(prev_vec[thread_id], mean_vec[thread_id], covariance[thread_id], &lpdf, work_vec[thread_id]);
                log_proposal_ratio[optim_id][particle_id] += lpdf;

                // Compute any of the dependent parameters for the distribution (if applicable).
                theta[optim_id][step_id + 1][dist_id]->CompleteSample(particle_id);
            }
        }
        else // Random walk sampler (i.e. does depend on current particle position).
        {
            // Otherwise create a random walk from the current particle using a multivariate gaussian.
            gsl_ran_multivariate_gaussian(rng[thread_id], zero_vec[thread_id], rw_covariance[thread_id], result_vec[thread_id]);

            // Retrieve the sample and ammend to the next proposed sample.
            int par_id = 0;
            for (int dist_id = 0; dist_id < n_dist; ++dist_id)
            {
                // Get the number of dimensions in the first distribution.
                int n_free_dim1 = theta[optim_id][step_id][dist_id]->GetNumFreeDims();
                int n_dim1 = theta[optim_id][step_id][dist_id]->GetNumDims();
                for (int dim1_id = 0; dim1_id < n_free_dim1; ++dim1_id)
                {
                    theta[optim_id][step_id + 1][dist_id]->chain[particle_id * n_dim1 + dim1_id] = gsl_vector_get(result_vec[thread_id], par_id);
                    theta[optim_id][step_id + 1][dist_id]->chain[particle_id * n_dim1 + dim1_id] += theta[optim_id][step_id][dist_id]->chain[particle_id * n_dim1 + dim1_id];
                    ++par_id;
                }

                // Compute any of the dependent parameters for the distribution (if applicable).
                theta[optim_id][step_id + 1][dist_id]->CompleteSample(particle_id);

                // Note since we are using a random walk for this particle we do not need to update the
                // ratio of proposals.
            }
        }
    }
}

// Function to determine is a particle is valid.
// i.e. within the support of its underlying distribution.
bool Optimise::IsParticleValid(int optim_id, int step_id, int particle_id)
{
    // Begin by assuming that the particle is valid.
    bool is_valid = true;

    // For each underlying distribution in the particle, check if
    // it is valid using the distributions member function.
    for (int dist_id = 0; dist_id < n_dist; ++dist_id)
    {
        is_valid &= theta[optim_id][step_id][dist_id]->IsValidSample(particle_id);
    }

    // Return whether or not the particle is valid.
    return is_valid;
}

// Function to set the simulation parameters from the posterior.
void Optimise::SetSimPostPars(int thread_id, int posterior_id)
{
    // Set the parameters in the simulation.
    pars[thread_id]->p_immune_start[0] = posterior->p_immune_start_anj.chain[posterior_id];
    pars[thread_id]->p_immune_start[1] = posterior->p_immune_start_gra.chain[posterior_id];
    pars[thread_id]->p_immune_start[2] = posterior->p_immune_start_may.chain[posterior_id];
    pars[thread_id]->p_immune_start[3] = posterior->p_immune_start_moh.chain[posterior_id];
    pars[thread_id]->move[0][1] = posterior->move_anj_gra.chain[posterior_id];
    pars[thread_id]->move[0][2] = posterior->move_anj_may.chain[posterior_id];
    pars[thread_id]->move[0][3] = posterior->move_anj_moh.chain[posterior_id];
    pars[thread_id]->move[1][0] = posterior->move_gra_anj.chain[posterior_id];
    pars[thread_id]->move[1][3] = posterior->move_gra_moh.chain[posterior_id];
    pars[thread_id]->move[3][0] = posterior->move_moh_anj.chain[posterior_id];
    pars[thread_id]->move[3][1] = posterior->move_moh_gra.chain[posterior_id];
    pars[thread_id]->import_start = posterior->import_start.chain[posterior_id];
    pars[thread_id]->import_duration = posterior->import_duration.chain[posterior_id];
    pars[thread_id]->import_size = posterior->import_size.chain[posterior_id];

    // The model parameters used depend on configuration.
    switch (RVF_DIFF_NDVI)
    {
    case 0: // Same NDVI per island.
        pars[thread_id]->ndvi_rate[0] = posterior->ndvi_rate.chain[posterior_id];
        pars[thread_id]->ndvi_rate[1] = posterior->ndvi_rate.chain[posterior_id];
        pars[thread_id]->ndvi_rate[2] = posterior->ndvi_rate.chain[posterior_id];
        pars[thread_id]->ndvi_rate[3] = posterior->ndvi_rate.chain[posterior_id];
        break;
    case 1: // Different NDVI per island.
        pars[thread_id]->ndvi_rate[0] = posterior->ndvi_rate_anj.chain[posterior_id];
        pars[thread_id]->ndvi_rate[1] = posterior->ndvi_rate_gra.chain[posterior_id];
        pars[thread_id]->ndvi_rate[2] = posterior->ndvi_rate_may.chain[posterior_id];
        pars[thread_id]->ndvi_rate[3] = posterior->ndvi_rate_moh.chain[posterior_id];
        break;
    default: // Default is same.
        pars[thread_id]->ndvi_rate[0] = posterior->ndvi_rate.chain[posterior_id];
        pars[thread_id]->ndvi_rate[1] = posterior->ndvi_rate.chain[posterior_id];
        pars[thread_id]->ndvi_rate[2] = posterior->ndvi_rate.chain[posterior_id];
        pars[thread_id]->ndvi_rate[3] = posterior->ndvi_rate.chain[posterior_id];
    }
    switch (RVF_DIFF_SCALE)
    {
    case 0: // Same transmission scale per island.
        pars[thread_id]->trans_scale[0] = posterior->trans_scale.chain[posterior_id];
        pars[thread_id]->trans_scale[1] = posterior->trans_scale.chain[posterior_id];
        pars[thread_id]->trans_scale[2] = posterior->trans_scale.chain[posterior_id];
        pars[thread_id]->trans_scale[3] = posterior->trans_scale.chain[posterior_id];
        break;
    case 1: // Different transmission scale per island.
        pars[thread_id]->trans_scale[0] = posterior->trans_scale_anj.chain[posterior_id];
        pars[thread_id]->trans_scale[1] = posterior->trans_scale_gra.chain[posterior_id];
        pars[thread_id]->trans_scale[2] = posterior->trans_scale_may.chain[posterior_id];
        pars[thread_id]->trans_scale[3] = posterior->trans_scale_moh.chain[posterior_id];
        break;
    default: // Default is same.
        pars[thread_id]->trans_scale[0] = posterior->trans_scale.chain[posterior_id];
        pars[thread_id]->trans_scale[1] = posterior->trans_scale.chain[posterior_id];
        pars[thread_id]->trans_scale[2] = posterior->trans_scale.chain[posterior_id];
        pars[thread_id]->trans_scale[3] = posterior->trans_scale.chain[posterior_id];
    }
}

// Function to update the simulation parameters with the desired
// optimisation parameters on a given optimisation chain.
void Optimise::SetSimOptimPars(int thread_id, int optim_id, int step_id, int particle_id)
{
    // Declare the proportion of vaccines to assign to each island.
    double vac_prop_anj, vac_prop_gra, vac_prop_may, vac_prop_moh;

    // Set the distribution of vaccines. This distribution will depend on how
    // many parameters are in the Dirchlet distribution.
    int k = vac_dist[optim_id][step_id]->GetPar(0); // Number of parameters in the Dirichlet distribution.
    if (k == 2)
    {
        // The scenario distribution determines the proportions for GC+Moh (by pop) and Anj+May (by pop).
        // Set the population weight of Anjouan between Anjouan and Mayotte.
        double anj_weight = pars[thread_id]->n_pop[0] / (pars[thread_id]->n_pop[0] + pars[thread_id]->n_pop[2]);
        double gra_weight = pars[thread_id]->n_pop[1] / (pars[thread_id]->n_pop[1] + pars[thread_id]->n_pop[3]);

        // Order of the islands is Anjouan, Grande Comore, Mayotte and Moheli.
        vac_prop_anj = anj_weight * vac_dist[optim_id][step_id]->chain[particle_id * k];
        vac_prop_gra = gra_weight * vac_dist[optim_id][step_id]->chain[particle_id * k + 1];
        vac_prop_may = (1.0 - anj_weight) * vac_dist[optim_id][step_id]->chain[particle_id * k];
        vac_prop_moh = (1.0 - gra_weight) * vac_dist[optim_id][step_id]->chain[particle_id * k + 1];
    }
    else if (k == 3)
    {
        // The scenario distribution determines the proportions for Grande Comore, Moheli and Anj+May (by pop).
        // Set the population weight of Anjouan between Anjouan and Mayotte.
        double anj_weight = pars[thread_id]->n_pop[0] / (pars[thread_id]->n_pop[0] + pars[thread_id]->n_pop[2]);

        // Order of the islands is Anjouan, Grande Comore, Mayotte and Moheli.
        vac_prop_anj = anj_weight * vac_dist[optim_id][step_id]->chain[particle_id * k];
        vac_prop_gra = vac_dist[optim_id][step_id]->chain[particle_id * k + 1];
        vac_prop_may = (1.0 - anj_weight) * vac_dist[optim_id][step_id]->chain[particle_id * k];
        vac_prop_moh = vac_dist[optim_id][step_id]->chain[particle_id * k + 2];
    }
    else if (k == 4)
    {
        // This is simple. The proportions at which to ditribute the vaccine are determined by the distribution itself.
        vac_prop_anj = vac_dist[optim_id][step_id]->chain[particle_id * k];
        vac_prop_gra = vac_dist[optim_id][step_id]->chain[particle_id * k + 1];
        vac_prop_may = vac_dist[optim_id][step_id]->chain[particle_id * k + 2];
        vac_prop_moh = vac_dist[optim_id][step_id]->chain[particle_id * k + 3];
    }
    else
    {
        std::cout << "Requested optimisation of a Dirichlet distribution with dimension that is incompatible with simulation." << std::endl;
        std::exit(-10);
    }

    // Copy across the vaccine proportions to the transformed scenario parameters.
    pars[thread_id]->vac_prop[0] = vac_prop_anj;
    pars[thread_id]->vac_prop[1] = vac_prop_gra;
    pars[thread_id]->vac_prop[2] = vac_prop_may;
    pars[thread_id]->vac_prop[3] = vac_prop_moh;
}

/*
// Function to write the sampled posterior parameters to file.
void Optimise::WritePostPars(int series_id)
{
    // Declare the file stream.
    std::fstream file;

    // Open a file for writing.
    char buffer[4]; // Three digits and the terminating character.
    sprintf(buffer, "%.3d", series_id);  // Format each number with some leading zero.
    file.open(static_cast<std::string>(RVF_ODIR) + "posterior_" + buffer + ".csv", std::fstream::out);

    // If we're at the start, open the file without appending.
    if (file.is_open())
    {
        // Output headers to the file!
        file << "OPTIM_ID,STEP_ID,SAMPLE_ID,PAR_NAME,PAR_VALUE";

        // For each scenario tested, save its ID, the scenario parameter names
        // and their values.
        for (int optim_id = 0; optim_id < n_optim; ++optim_id)
        {
            for (int i = 0; i < RVF_NSAMPLES; ++i)
            {
                for (int j = 0; j < n_samples[optim_id][i]; ++j)
                {
                    for (int k = 0; k < posterior->n_pars; ++k)
                    {
                        file << "\n" << std::fixed << std::setprecision(3);
                        file << optim_id << "," << i << "," << j << "," << posterior->theta[k]->par_name;
                        file << "," << posterior->theta[k]->chain[post_idx[optim_id][i][j]];
                    }
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
*/

// Function to write the cost values to file.
void Optimise::WriteCost(int series_id)
{
    // Declare the file stream.
    std::fstream file;

    // Open a file for writing.
    char buffer[4];                     // Three digits and the terminating character.
    sprintf(buffer, "%.3d", series_id); // Format each number with some leading zero.
    file.open(static_cast<std::string>(RVF_ODIR) + "cost_" + buffer + ".csv", std::fstream::out);

    // If we're at the start, open the file without appending.
    if (file.is_open())
    {
        // Output headers to the file!
        file << "OPTIM_ID,STEP_ID,PARTICLE_ID,SAMPLE_ID,COST_VALUE";

        // For each scenario tested, save its ID, the scenario parameter names
        // and their values.
        for (int optim_id = 0; optim_id < n_optim; ++optim_id)
        {
            for (int step_id = 0; step_id <= max_step_reached[optim_id]; ++step_id)
            {
                for (int particle_id = 0; particle_id < n_particles; ++particle_id)
                {
                    for (int cost_id = 0; cost_id < n_cost_samples; ++cost_id)
                    {
                        file << "\n"
                             << optim_id << "," << step_id << "," << particle_id << "," << cost_id;
                        file << "," << cost[optim_id][step_id][particle_id][cost_id];
                    }
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
void Optimise::WriteOptimPars(int series_id)
{
    // Declare the file stream.
    std::fstream file;

    // Open a file for writing.
    char buffer[7];                     // Five digits and the terminating character.
    sprintf(buffer, "%.3d", series_id); // Format each number with some leading zero.
    file.open(static_cast<std::string>(RVF_ODIR) + "optim_pars_" + buffer + ".csv", std::fstream::out);

    // If we're at the start, open the file without appending.
    if (file.is_open())
    {
        // Output headers to the file!
        file << "OPTIM_ID,STEP_ID,PARTICLE_ID,PAR_NAME,INDEX,PAR_VALUE";

        // For each optimisation, save all particle information.
        for (int optim_id = 0; optim_id < n_optim; ++optim_id)
        {
            // For each step in the SMC...
            for (int step_id = 0; step_id <= max_step_reached[optim_id]; ++step_id)
            {
                // For every particle in that distribution...
                for (int particle_id = 0; particle_id < n_particles; ++particle_id)
                {
                    // For each distribution in theta...
                    for (int dist_id = 0; dist_id < n_dist; ++dist_id)
                    {
                        // Get the number of dimensions in the distribution.
                        int n_dim = theta[optim_id][step_id][dist_id]->GetNumDims();
                        for (int dim_id = 0; dim_id < n_dim; ++dim_id)
                        {
                            file << "\n"
                                 << optim_id << "," << step_id << "," << particle_id;
                            file << "," << theta[optim_id][step_id][dist_id]->par_name << "," << dim_id;
                            file << "," << theta[optim_id][step_id][dist_id]->chain[particle_id * n_dim + dim_id];
                        }
                    }
                }
            }
        }
    }
    else
    {
        std::cout << "Error: there was a problem opening the file for writing";
        std::cout << " optimisation parameter output." << std::endl;
        exit(1);
    }

    // Close the file after writing has finished.
    file.close();
}