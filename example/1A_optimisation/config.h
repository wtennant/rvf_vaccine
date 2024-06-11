// config.h: Configures macros for directories to read files in from.

#ifndef CONFIG_H
#define CONFIG_H

#define RVF_IDIR "in/"                      // Relative path (from executable) to the data.
#define RVF_POSTDIR "in/"                   // Relative path (from executable) to the folder containing the posterior.
#define RVF_ODIR "out/"                     // Relative path (from executable) to the output folder.
#define RVF_UPDATE_FREQ 10                  // Frequency to update user on scenario testing progress in console.
#define RVF_NSAMPLES 500                    // Define the number of scenarios/independent optimisations to run.
#define RVF_POST_NSAMPLES 50                // How many posterior samples to use per scenario.
#define RVF_MIN_NDVI 0                      // Minimum NDVI is local (0) or global (1) or nothing (2).
#define RVF_TRANSMISSION 2                  // Use a constant (0), linear (1), or exponential (2) transmission model.
#define RVF_DIFF_SCALE 1                    // Use the same (0) or different (1) transmission scalar for islands.
#define RVF_DIFF_NDVI 0                     // Use the same (0) or different (1) NDVI scalars for each island.
#define RVF_MAX_BUFFER 5*1024*1024          // Maximum size of a file buffer in bytes.
#define RVF_WRITE_MIN_PRECISION 3           // The minimum precision to record data values as in a file (at least this many decimals places + significant figures). 
#define RVF_WRITE_FULL 0                    // Should the values of all compartments be written?
#define RVF_WRITE_SUMMARY 0                 // Should summary values (sum across ages) of all compartments be written?
#define RVF_SUMMARY_SIZE 48                 // Ideally divides 48 (the number of epi-weeks in a year).
#define RVF_NO_VAC 1                        // Should a scenario without vaccination be run?
#define RVF_LONG_TERM 1                     // Should the long-term effects of vaccination be simulated?
#define RVF_OPTIMISE 1                      // Define whether we will optimise a defined summary metric over a given set of parameters.
#define RVF_OPTIM_RUNS 48                   // The number of optimisation runs to make.
#define RVF_OPTIM_NCOST 10                  // The number of repeats of the cost function to take.
#define RVF_OPTIM_NPARTICLES 50             // The number of particles to use per optimisiation run.
#define RVF_OPTIM_MAX_STEPS 20              // Define the maximum number of steps permitted in the optimisation algorithm.
#define RVF_OPTIM_SUMMARY 0                 // Summary function to use: mean (0), minimum (1) or maximum (2).
#define RVF_OPTIM_NDELTA 100                // Number of tempering parameters that can be used.
#define RVF_OPTIM_MAX_MH 25                 // Maximum number of steps in the Metropolis-Hastings kernels for jittering particles.
#define RVF_OPTIM_FUNC 0                    // Optimisation function to use: maximise overall infections averted (0), maximise infections averted on worst affected island (1).
#define RVF_OPTIM_POWER 4                   // Power to raise costs to in order to widen/narrow regions of optimality.
#define RVF_NTHREADS 1                     // Number of threads to use in parallelism.

// Define a function which outputs the configuration to file.
void WriteConfig();

#endif