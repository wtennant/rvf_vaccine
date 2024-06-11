
#include "gsl/gsl_rng.h"
#include "gsl/gsl_vector.h"
#include "gsl/gsl_matrix.h"
#include "parameters.h"
#include "posterior.h"
#include "simulation.h"
#include "theta.h"

#ifndef OPTIMISE_H
#define OPTIMISE_H

class Optimise
{
public:
    // Run the optimisation algorithm on the desired parameter distributions.
    void Run(int);

    // Constructor and deconstructor of the optimise class.
    Optimise(Parameters*, Data*, Posterior*);
    ~Optimise();

private:

    // Distributions of parameters to be optimised.
    // Parameters for these distributions are defined in optimise.cpp.
    Dirichlet*** vac_dist;

    // Space to store all optimisation parameter distributions.
    Theta**** theta;

    // Declare space to store a series of simulations, data and parameters.
    Simulation** sim;
    Parameters** pars;
    Data** data;

    // Declare space to store the addresses of the posterior distributions.
    Posterior* posterior;

    // Declare a pointer to a series of random number generators.
    gsl_rng** rng;

    // A covariance matrix for the previously accepted values of the
    // optimisation algorithm. We also need vectors for a zero vector,
    // sample vector, work vector and previous sample vector from a
    // multivariate-gaussian.
    gsl_matrix** covariance;
    gsl_matrix** rw_covariance;
    gsl_vector** mean_vec;
    gsl_vector** zero_vec;
    gsl_vector** result_vec;
    gsl_vector** prev_vec;
    gsl_vector** work_vec;

    // An array to keep track of the evaluation of the cost function of accepted
    // Indexing: optimisation run, particle in optimisation run,
    //           step in optimisation run, sample id of cost function.
    double**** cost;
    double*** summary_cost;
    double*** tempered_cost;

    // An array to keep track of the logged ratio of proposals during
    // the particle rejuvination step.
    // Indexing: optimisation run, particle in optimisation run.
    double** log_proposal_ratio;

    // An array of the weights of each particle per step in the optimisation run.
    double *** weight;

    // ID of the posterior parameters used in the optimisation.
    // Indexing: optimisation run, particle in optimisation run,
    //           step in optimisation run, sample id of cost function.
    int**** post_idx;

    // Vector to keep track of the tempering parameter and the associated index
    // in the list of possible tempering parameters for each optimisation run and step.
    double** tempering_delta;
    int** tempering_delta_idx;

    // An array for the adjusted weights of each particle in the optimisation run
    // for all possible tempering parameter options.
    double** proposed_weights;

    // Vector of possible tempering parameters.
    double* delta_opt;

    // The maximum step reached by each optimisation run.
    int* max_step_reached;

    // Fixed integers used to define:
    int n_dist;         // The number of parameter distributions to optimise.
    int n_particles;    // The number of particles to use in a single optimisation.
    int n_pars;         // The total numbr of free optimisation paraemters.
    int n_optim;        // The number of optimisation chains to run in parallel.
    int n_threads;      // The number of CPU-based threads available for parallelisation.
    int n_cost_samples; // The number of repeated sampels from the cost function for a fixed particle.
    int n_steps;        // The number of steps in a single optimisation run.
    int n_delta;        // The number of possible tempering factors (excluding zero).

    // Whether or not the no-vaccination scenario has been already calculated
    // for a given posterior-index.
    // Also create space to temporarily store vaccination costs.
    bool* is_no_vac_costed;
    double* no_vac_cost;
    double* vac_cost;

    // Function to calculate a metric given simulation data.
    // This is our cost function.
    double CalculateCost(int, int, int, int);

    // Function to set the posterior parameters of the simulation
    // on a given optimisation chain.
    void SetSimPostPars(int, int);

    // Function to update the optimisation parameters of 
    // the simulation on a given optimisation chain.
    void SetSimOptimPars(int, int, int, int);

    // Function to summarise multiple samples from the cost function.
    double SummariseCost(int, int, int);

    // Function which tempers the cost of a given strategy using a pre-set
    // tempering parameter
    double TemperCost(int, int, int, double);

    // Function to shuffle the particles at the current step based on their
    // updated weight.
    void ResampleParticles(int, int);

    // Function to jitter a set of particles using a Metropolis-Hastings kernel.
    void JitterParticles(int, int);

    // Function to determine is a proposed set of parameters is valid.
    bool IsParticleValid(int, int, int);

    // Function to accept a candidate parameter set, and overwrite the candidate
    // parameter set in the chain by the proposed parameter set.
    void AcceptCandidate(int, int, int);

    // Function to write the used posterior samples to file.
    void WritePostPars(int);

    // Function to write the accepted parameters in the simulated
    // annealing algorithm to file.
    void WriteOptimPars(int);

    // Function to write the cost values to file
    // for a particular optimisation chain.
    void WriteCost(int);
};

#endif