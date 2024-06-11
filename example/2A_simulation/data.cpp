// data.cpp: Contains the constructor and destructor of the data classes.
#include <fstream>      // For reading in the data from file.
#include <iostream>     // Input-output to console.
#include "config.h"     // Configuration for the RVF model.
#include "data.h"       // Definition of the data classes.

// Function to read in the NDVI data.
void Ndvi::ReadNDVI(Parameters* default_pars)
{
    // Create a file stream.
    std::fstream file;
    
    // Open the file to read in.
    if (RVF_LONG_TERM == 0)
    {
        file.open(static_cast<std::string>(RVF_IDIR) + "ndvi.csv", std::fstream::in);
    }
    else
    {
        file.open(static_cast<std::string>(RVF_IDIR) + "ndvi_extended.csv", std::fstream::in);
    }

    // Check that the file opened.
    if (file.is_open())
    {
        // Define the contents of a line.
        std::string line;

        // Throw away the first line of the data (the header).
        getline(file, line);

        // Get the number of lines in the file.
        int n_lines = 0;
        while (getline(file, line))
        {
            ++n_lines;
        }

        // Set the length of the simulation.
        // There are n_island times as many lines as simulation steps.
        default_pars->n_steps = n_lines / n_island;
        n_time_steps = default_pars->n_steps;

        // For each island, allocate space to store the data based on the number
        // of lines in the file. The data is in months, so multiply by four to
        // get the data in weeks.
        for (int i = 0; i < n_island; ++i)
        {
            ndvi[i] = new double[n_time_steps];
        }

        // Reset the position of the file stream to the first line.
        file.clear();
        file.seekg(0, std::ios::beg);

        // Throw away the first line of the data (the header) again.
        getline(file, line);

        // Keep track of the number of entries scanned, and the
        // number of entries read in to the data.
        int entries = 0;
        int* read_entries = new int[n_island]{0};

        // Define the number of columns in the data.
        int n_cols = 6;

        // Define the separator of entries.
        char sep = ',';

        // Variable to keep track of the island ID number.
        int island_id = 0;

        // For each line, extract the final entry in the data.
        // The final entry is the NDVI data.
        while(getline(file, line, sep))
        {       
            // Increase the counter for the number of entries that have been read.
            ++entries;

            // If this value is divisible by four, then the current line contains
            // the ndvi data.
            // On the first pass, record the month and year that the NDVI data
            // starts at.
            if (entries % n_cols == 0)
            {
                // Write the entry to the NDVI time-series for each
                // epidemiological week.
                ndvi[island_id][read_entries[island_id]] = std::stod(line);

                // Store the minimum NDVI value if necessary.
                if (ndvi[island_id][read_entries[island_id]] < min_local_ndvi[island_id])
                {
                    min_local_ndvi[island_id] = ndvi[island_id][read_entries[island_id]];

                    // Also check if this is a new global minimum.
                    if (min_local_ndvi[island_id] < min_global_ndvi)
                    {
                        min_global_ndvi = min_local_ndvi[island_id];
                    }
                }

                // Increase the number of entries that have been read.
                ++read_entries[island_id];
            }
            else if (entries == 4)
            {
                // Get the starting month of the data.
                start_month = std::stoi(line);
            }
            else if (entries == 5)
            {
                // Get the starting year of the data.
                start_year = std::stoi(line);
            }
            else if ((entries - 1) % n_cols == 0)
            {
                // Get the island ID number for the row.
                island_id = std::stoi(line);
            }

            // Define the next separator character.
            // If about to read the final entry of a line, change the separator
            // to a new line.
            sep = (entries + 1) % n_cols == 0 ? '\n' : ',';
        }

        // Close the file.
        file.close();
    }
    else
    {
        // Write an error messsage and terminate the execution of the program.
        std::cout << "Error: NDVI file could not be read." << std::endl;
        exit(1);
    }
}

// Constructor for the seroprevalence class.
Ndvi::Ndvi(Parameters* default_pars)
{
    // Store the number of islands.
    n_island = default_pars->n_island;

    // Set up space to store the NDVI data for each island.
    ndvi = new double*[n_island];

    // Set up space to store the minimum NDVI across all islands and for each island.
    min_global_ndvi = 1.5;
    min_local_ndvi = new double[n_island];

    // ... and initialise.
    for (int i = 0; i < n_island; ++i)
    {
        min_local_ndvi[i] = 1.5;
    }

    // Set up space to store the NDVI data and read in the NDVI data.
    // Also set the length of the simulation.
    ReadNDVI(default_pars);
}

// Deconstructor for the seroprevalence class.
Ndvi::~Ndvi()
{
    // De-allocate space for the seroprevalence data.
    delete[] min_local_ndvi;
    for (int i = 0; i < n_island; ++i)
    {
        delete[] ndvi[i];
    }
    delete[] ndvi;
}

Ndvi::Ndvi(const Ndvi& default_ndvi)
{
    // Store the number of islands.
    n_island = default_ndvi.n_island;

    // Set up space to store the NDVI data for each island.
    ndvi = new double*[default_ndvi.n_island];
    for (int i = 0; i < default_ndvi.n_island; ++i)
    {
        ndvi[i] = new double[default_ndvi.n_time_steps];
        for (int t = 0; t < default_ndvi.n_time_steps; ++t)
        {
            ndvi[i][t] = default_ndvi.ndvi[i][t];
        }
    }

    // Set up space to store the minimum NDVI across all islands and for each island.
    min_global_ndvi = default_ndvi.min_global_ndvi;
    min_local_ndvi = new double[default_ndvi.n_island];

    // ... and initialise.
    for (int i = 0; i < default_ndvi.n_island; ++i)
    {
        min_local_ndvi[i] = default_ndvi.min_local_ndvi[i];
    }
}

// Constructor for the data class.
Data::Data(Parameters* default_pars) : ndvi(default_pars)
{
    // Set the starting year and month equal to the first
    // entry in the NDVI data.
    start_month = ndvi.start_month;
    start_year = ndvi.start_year;
}

// Constructor for a deep copy implementation.
Data::Data(const Data& default_data) : ndvi(default_data.ndvi) 
{
    // Set the starting year and month equal to the first
    // entry in the NDVI data.
    start_month = default_data.ndvi.start_month;
    start_year = default_data.ndvi.start_year;
}
