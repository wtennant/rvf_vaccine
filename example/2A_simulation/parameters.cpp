// parameters.cpp: Contains the constructor and destructor of the
// class, Parameters.
#include <fstream>          // File stream.
#include <iostream>         // Input-output stream.
#include "config.h"         // Output directory.
#include "parameters.h"     // Definition of the Parameters class.

// Default constructor for the Parameters class.
Parameters::Parameters()
{
    // Number of livestock age-groups (grouped by year).
    // Note: Currently undefined behaviour if n_age is not equal to 10.
    n_age = 10;

    // Maximum age group to vaccinate.
    n_age_vac = 10;

    // Maximum age group which can be moved.
    n_age_move = 2;

    // Define the number of islands in the metapopulation model.
    n_island = 4;

    // Number of livestock in the population.
    n_pop = new double[n_island];
    n_pop[0] = 93616.0;     // Anjouan.
    n_pop[1] = 224353.0;    // Grande Comore.
    n_pop[2] = 20052.0;     // Mayotte.
    n_pop[3] = 31872.0;     // Moheli.

    // Number of time steps to simulate (epi-weeks).
    n_steps = 576;

    // Proportion of animals that are immune at time zero on each island.
    p_immune_start = new double[n_island]{0.1298};

    // Scalars for how NDVI influences the force of infection (FoI).
    // FoI = 1 - exp(-exp(-rate*NDVI + constant)*I/N).
    trans_scale = new double[n_island]{-2.19};
    ndvi_rate = new double[n_island]{3.15};

    // Declare space for the number of movements between islands.
    move = new double*[n_island];
    for (int i = 0; i < n_island; ++i)
    {
        move[i] = new double[n_island];
        for (int j = 0; j < n_island; ++j)
        {
            move[i][j] = 0.0;
        }
    }

    // Declare space for each of the age-dependent parameters.
    // Proportion of the population in each age group.
    pop_structure = new double[n_age]{0.294, 0.206, 0.144, 0.101, 0.071,
                                      0.050, 0.035, 0.024, 0.017, 0.058};

    // Mortality rate for each age group.
    mortality = new double[n_age]{0.0088, 0.0088, 0.0088, 0.0088, 0.0088,
                                  0.0088, 0.0088, 0.0088, 0.0088, 0.0062};

    // Whether (> 0) or not (== 0) only unvaccinated animals are vaccinated.
    vac_identifiable = 0;           

    // Efficacy of the vaccine (as a proportion).
    vac_efficacy = 0.9;

    // Define the default distribution of vaccines across the archipelago.
    vac_prop = new double[n_island];
    double total_pop = 0;
    for (int i = 0; i < n_island; ++i){ total_pop += n_pop[i]; }
    for (int i = 0; i < n_island; ++i){ vac_prop[i] = n_pop[i] / total_pop; }

    // Define the default number of vaccines to administer across the archipelago.
    vac_rate = 0.15*total_pop;

    // Duration of vaccine protection.
    vac_protect_duration = 2*48;

    // Define the vaccine timing parameters.
    vac_t_first = 11*48;    // The earliest time point to implement the vaccine.
    vac_t_year_offset = 0;  // The first time point each year to implement the vaccine (0--47).
    vac_t_year_length = 48; // The length of the annual vaccine implementation (1--48).
    vac_t_freq = 1;         // Frequency of the vaccination window.

    // Define the frequency of imports into Grande Comore.
    reg_import_t_start = 11*48;
    reg_import_freq = 10*48;
    reg_import_scale = 1.0;
}

// Default deconstructor for the Parameter class.
Parameters::~Parameters()
{
    // Free dynamically allocated memory.
    delete[] vac_prop;
    delete[] mortality;
    delete[] pop_structure;
    for (int i = 0; i < n_island; ++i)
    {
        delete[] move[i];
    }
    delete[] move;
    delete[] ndvi_rate;
    delete[] trans_scale;
    delete[] p_immune_start;
    delete[] n_pop;
}

// Deep copy constructor for a given set of parameters.
Parameters::Parameters(const Parameters& default_pars)
{
    // Number of livestock age-groups (grouped by year).
    // Note: Currently undefined behaviour if n_age is not equal to 10.
    n_age = default_pars.n_age;

    // Maximum age to vaccinate.
    n_age_vac = default_pars.n_age_vac;

    // Maximum age group which can be moved.
    n_age_move = default_pars.n_age_move;

    // Define the number of islands in the metapopulation model.
    n_island = default_pars.n_island;

    // Number of livestock in the population.
    n_pop = new double[n_island];
    for (int i = 0; i < n_island; ++i)
    {
        n_pop[i] = default_pars.n_pop[i];
    }

    // Number of time steps to simulate (epi-weeks).
    n_steps = default_pars.n_steps;

    // Proportion of animals that are immune at time zero on each island.
    p_immune_start = new double[n_island];

    // Scalars for how NDVI influences the force of infection (FoI).
    // FoI = 1 - exp(-exp(-rate*NDVI + constant)*I/N).
    trans_scale = new double[n_island];
    ndvi_rate = new double[n_island];

    // The number of vaccines and duration of protection.
    vac_rate = default_pars.vac_rate;
    vac_protect_duration = default_pars.vac_protect_duration;

    // Define space for the distribution of vaccines to administer on each
    // island per epidemiological week.
    vac_prop = new double[n_island];

    // Declare space for the number of movements between islands.
    move = new double*[n_island];
    for (int i = 0; i < n_island; ++i)
    {
        move[i] = new double[n_island];

        // Perform deep copies of other n_island-dependent objects here.
        p_immune_start[i] = default_pars.p_immune_start[i];
        trans_scale[i] = default_pars.trans_scale[i];
        ndvi_rate[i] = default_pars.ndvi_rate[i];
        vac_prop[i] = default_pars.vac_prop[i];
        for (int j = 0; j < n_island; ++j)
        {
            move[i][j] = default_pars.move[i][j];
        }
    }

    // Declare space for each of the age-dependent parameters.
    // Proportion of the population in each age group.
    pop_structure = new double[n_age];

    // Mortality rate for each age group.
    mortality = new double[n_age];

    // Deep copy across n_age-dependent variables.
    for (int age = 0; age < n_age; ++age)
    {
        pop_structure[age] = default_pars.pop_structure[age];
        mortality[age] = default_pars.mortality[age];
    }

    // Whether (> 0) or not (== 0) only unvaccinated animals are vaccinated.
    vac_identifiable = default_pars.vac_identifiable;  

    // Efficacy of the vaccine (as a proportion).
    vac_efficacy = default_pars.vac_efficacy;

    // Vaccine timings.
    vac_t_first = default_pars.vac_t_first;
    vac_t_year_length = default_pars.vac_t_year_length;
    vac_t_year_offset = default_pars.vac_t_year_offset;
    vac_t_freq = default_pars.vac_t_freq;

    // Regular imports into Grande Comore.
    reg_import_freq = default_pars.reg_import_freq;
    reg_import_scale = default_pars.reg_import_scale;
    reg_import_t_start = default_pars.reg_import_t_start;
}

// Function to calculate the total population size across all islands.
double Parameters::CalculateTotalPop()
{
    double total_pop = 0.0;
    for (int i = 0; i < n_island; ++i){ total_pop += n_pop[i]; }
    return total_pop;
}

// Write default parameters to file.
void Parameters::WriteDefaultPars(int series_id)
{
    // Declare the file stream.
    std::fstream file;

    // Open the file to record parameters.
    char buffer[4]; // Three digits and a terminating character.
    sprintf(buffer, "%.3d", series_id); // Leading zeros on the series id.
    file.open(static_cast<std::string>(RVF_ODIR) + "default_pars_" + buffer + ".csv", std::fstream::out);

    // Check that the file is open.
    if (file.is_open())
    {
        // Record the header for the default parameter file.
        file << "PAR_NAME,PAR_INDEX,PAR_VALUE" << std::endl;

        // For each parameter, output the name and value.
        file << "n_age,0," << n_age << std::endl;
        file << "n_age_vac,0," << n_age_vac << std::endl;
        file << "n_age_move,0," << n_age_move << std::endl;
        file << "n_island,0," << n_island << std::endl;
        file << "n_steps,0," << n_steps << std::endl;
        for (int i = 0; i < n_island; ++i)
        {
            file << "n_pop," << i << "," << n_pop[i] << std::endl;
        }
        for (int age = 0; age < n_age; ++age)
        {
            file << "pop_structure," << age << "," << pop_structure[age] << std::endl;
        }
        for (int age = 0; age < n_age; ++age)
        {
            file << "mortality," << age << "," << mortality[age] << std::endl;
        }
        file << "vac_identifiable,0," << vac_identifiable << std::endl;
        file << "vac_efficacy,0," << vac_efficacy << std::endl;
        file << "vac_protect_duration,0," << vac_protect_duration << std::endl;
        file << "vac_rate,0," << vac_rate << std::endl;
        for (int i = 0; i < n_island; ++i)
        {
            file << "vac_prop," << i << "," << vac_prop[i] << std::endl;
        }
        file << "vac_t_first,0," << vac_t_first << std::endl;
        file << "vac_t_year_length,0," << vac_t_year_length << std::endl;
        file << "vac_t_year_offset,0," << vac_t_year_offset << std::endl;
        file << "vac_t_freq,0," << vac_t_freq << std::endl;
        file << "reg_import_freq,0," << reg_import_freq << std::endl;
        file << "reg_import_scale,0," << reg_import_scale << std::endl;
        file << "reg_import_t_start,0," << reg_import_t_start << std::endl;
    }
    else
    {
        // If the file could not be opened, throw a wobbly.
        std::cout << "File could not be opened to record default parameter set." << std::endl;
        exit(1);
    }

    // Close the file after writing has finished.
    file.close();    
}