// posterior.cpp: Defines the class to read in the posterior distribution
// of parameters from a previous model fit.
#include <algorithm>            // Iterator for finding characters in a string.
#include <fstream>              // Writing to file.
#include <iostream>             // Output errors to console.
#include <string>               // Removing portions of a string.
#include "config.h"             // Configuration of the model.
#include "posterior.h"          // Posterior class definition.

// Constructor for the Marginal class.
Marginal::Marginal(std::string in_par_name) : par_name(in_par_name)
{
    // The constructor sets the chains to the null pointer.
    chain = nullptr;
}

// Deconstructor for the Marginal class.
Marginal::~Marginal()
{
    // De-allocates the chain. This is safe on a null pointer.
    delete[] chain;
}

// Constructor of the Posterior class.
Posterior::Posterior()
{
    // Define the number of parameters.
    n_pars = 14;    // The number of parameters independent of model choice.

    // Add in the number of parameters depending on:
    // different ndvi influence on transmission per island (or not),
    // and different transmission scales per island (or not).
    switch(RVF_DIFF_NDVI)
    {
        case 0: // NDVI influence is the same on each island.
            n_pars += 1;
            break;
        case 1: // NDVI influence is different on each island.
            n_pars += 4;
            break;
        default:
            n_pars += 1;
    }
    switch(RVF_DIFF_SCALE)
    {
        case 0: // Transmission scale is the same on each island.
            n_pars += 1;
            break;
        case 1: // Transmission scale is different on each island.
            n_pars += 4;
            break;
        default:
            n_pars += 1;
    }

    // Space to store all the parameters and their chains.
    theta = new Marginal*[n_pars]{&p_immune_start_anj,
                                  &p_immune_start_gra,
                                  &p_immune_start_may,
                                  &p_immune_start_moh,
                                  &move_anj_gra,
                                  &move_anj_may,
                                  &move_anj_moh,
                                  &move_gra_anj,
                                  &move_gra_moh,
                                  &move_moh_anj,
                                  &move_moh_gra,
                                  &import_size,
                                  &import_start,
                                  &import_duration};

    // This will depend on the model: whether NDVI influence, and/or transmission
    // scale, are the same or different on each island.
    switch(RVF_DIFF_NDVI)
    {
        case 1:
            theta[14] = &ndvi_rate_anj;
            theta[15] = &ndvi_rate_gra;
            theta[16] = &ndvi_rate_may;
            theta[17] = &ndvi_rate_moh;
            switch(RVF_DIFF_SCALE)
            {
                case 1:
                    theta[18] = &trans_scale_anj;
                    theta[19] = &trans_scale_gra;
                    theta[20] = &trans_scale_may;
                    theta[21] = &trans_scale_moh;
                    break;
                default:
                    theta[18] = &trans_scale;
            }
            break;
        default:
            theta[14] = &ndvi_rate;
            switch(RVF_DIFF_SCALE)
            {
                case 1:
                    theta[15] = &trans_scale_anj;
                    theta[16] = &trans_scale_gra;
                    theta[17] = &trans_scale_may;
                    theta[18] = &trans_scale_moh;
                    break;
                default:
                    theta[15] = &trans_scale;
            }
    }

    // Output to information to console.
    std::cout << "Reading in posterior samples...";

    // Read in the posterior parameters.
    ReadPosterior();

    // Inform the user that the reading in of posterior samples is complete.
    std::cout << " complete." << std::endl;
}

// Deconstructor of the Posterior class.
Posterior::~Posterior()
{
    // Free up the space for the sampled posterior distribution.
    delete[] theta;
}

// Function to read in the posterior distribution.
void Posterior::ReadPosterior()
{
    // Create a file stream.
    std::fstream file;
    
    // Open the file to read in.
    file.open(static_cast<std::string>(RVF_POSTDIR) + "posterior.csv",
              std::fstream::in);

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

        // Define the number of samples from the posterior distribution.
        n_samples = n_lines / n_pars;

        // For each parameter in the posterior distribution, allocate space
        // to store the samples.
        for (int i = 0; i < n_pars; ++i)
        {
            theta[i]->chain = new double[n_samples];
        }

        // Reset the position of the file stream to the first line.
        file.clear();
        file.seekg(0, std::ios::beg);

        // Throw away the first line of the data (the header) again.
        getline(file, line);

        // Keep track of the number of individual entries scanned, and the
        // number of posteriors samples read in.
        int entries = 0;
        int sample = 0;

        // Define an index to keep track of which parameter from the posterior
        // is being read in.
        int theta_idx = 0;

        // Define the number of columns in the data.
        int n_cols = 3;

        // Define the separator of entries.
        char sep = ',';

        // For each line, extract the parameter name and value.
        while(getline(file, line, sep))
        {       
            // Increase the counter for the number of entries that have been read.
            ++entries;

            // If this value + 1 is divisible by 3, then the entry is the
            // parameter name. If the value itself is divisible by three,
            // then the entry is the parameter value.
            if ((entries + 1) % n_cols == 0) 
            {
                // Reset the index which tracks which parameter is being read in.
                theta_idx = -1;

                // Strip out the \" either side of the parameter name.
                line.erase(std::remove(line.begin(), line.end(), '\"'), line.end());
                
                // Check the parameter name in the file against stored
                // parameter names. Get the index of theta in which they agree.
                for (int i = 0; i < n_pars; ++i)
                {
                    if (theta[i]->par_name == line)
                    {
                        theta_idx = i;
                    }
                }

                // If the parameter name in the file does not exist
                // in theta, then something has gone wrong!
                if (theta_idx == -1)
                {
                    std::cout << "Parameter " << line << " cannot be found in class Posterior!" << std::endl;
                    exit(3);
                }
            }
            else if (entries % n_cols == 0) 
            {
                // Store the parameter value in the corresponding chain in theta.
                theta[theta_idx]->chain[sample] = std::stod(line);
            }

            // All parameters from a given sample have been read in when the
            // number of entries read is divisible by n_pars*n_cols.
            if (entries % (n_cols*n_pars) == 0)
            {
                ++sample;
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
        std::cout << "Error: posterior file could not be read." << std::endl;
        exit(1);
    }
}

// Function to re-write the posterior to file.
// This is useful for check it was read in correctly and
// for storage with the scenario test output.
void Posterior::WritePosterior()
{
    // Create a file stream.
    std::fstream file;

    // Open up a file stream for recording simulated data.
    file.open(static_cast<std::string>(RVF_ODIR) + "posterior.csv", std::fstream::out);

    // Check that the file is ready for writing.
    if (!file.is_open())
    {
        std::cout << "Error: there was a problem opening the file for writing";
        std::cout << " posterior output." << std::endl;
        exit(1);
    }
    else
    {
        // Write headers to file.
        file << "SAMPLE_ID,PAR_NAME,PAR_VALUE";

        // Go through all samples and write them to file.
        for (int i = 0; i < n_samples; ++i)
        {
            for (int j = 0; j < n_pars; ++j)
            {
                file << "\n" << i << "," << theta[j]->par_name << "," << theta[j]->chain[i];
            }
        }

        // Close the file.
        file.close();
    }
}