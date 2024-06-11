// simulation.h: Function delcaration for simulating the RVF model.
#include <fstream>          // File stream for recording data.
#include "gsl/gsl_rng.h"    // Random number generation.
#include "data.h"           // Definition of Data class.
#include "parameters.h"     // Definition of Parameters class.

#ifndef SIMULATION_H
#define SIMULATION_H

class Simulation
{
public:
    // Declare variables for the livestock population, split by infectious status:
    // Indexing: time, island, age group, and epidemiological group.
    double*** compartments;

    // Define the number of epidemiological groups.
    static const int n_epi = 16;
    static const int n_no_vac_epi = 4;

    // Counters.
    double*** doses;    // The number of administered vaccine doses.

    // Declare a file stream for recording simulation output.
    std::fstream file;
    std::fstream file_summary;
    std::fstream file_no_vac;
    std::fstream file_no_vac_summary;

    // Function to simulate the dynamics, including intiialisation,
    // between the desired time points.
    void Simulate(Data* data, Parameters* pars, int start_time, int end_time);

    // Function to write the output of the simulation.
    void WriteOutput(int iter, int rep);
    void WriteSummary(int iter, int rep);
    void WriteNoVac(int iter, int rep);
    void WriteNoVacSummary(int iter, int rep);

    // Constructor and destructor to allocate and de-allocate space for model
    // simulation.
    Simulation(Parameters* pars, int sens_id, int sim_id);
    ~Simulation();

private:
    int n_age;          // Store copy of the number of age groups.
    int n_age_move;     // Store copy of the maximum age group that can be moved.
    int n_island;       // Store copy of the number of islands.
    int n_steps;        // Store copy of the number of time-steps.
    double imp_scaling; // Scaling factor for working out the proportion of imports each age group.

    // Variable for storing the animals that are moved between islands.
    // Index ordering: source island, destination island, age+epi group.
    double*** compartments_move;

    // Unscoped enumeration delcaration for order of epidemiological compartments.
    enum {SU = 0, EU = 1, IU = 2, RU = 3, SV1 = 4, EV1 = 5, IV1 = 6, RV1 = 7, SV2 = 8, EV2 = 9, IV2 = 10, RV2 = 11, SW = 12, EW = 13, IW = 14, RW = 15};

    // Keep track of...
    double* net_move;           // Net movement per island.
    double* deaths;             // Number of deaths per island.
    double* total_move_N;       // Number of individuals who are permitted to move per island.

    // Buffer used to write to file.
    char* file_buffer;
    char* p_buffer;

    // Function to calculate the number of movements between islands at
    // a given time step.
    void CalculateMove(Parameters* pars,
                       const int t);

    // Function which formats values and places them in the buffer. The buffer is flushed to the file and reset
    // if the buffer would be overloaded by an update.
    template<typename... Args>
    void write_buffer(std::fstream& my_file, int* buffer_left, const char* format, Args... args);
};

#endif