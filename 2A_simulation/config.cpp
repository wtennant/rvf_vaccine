// config.cpp: Contains functions related to the configuration of the model
// and fitting procedure.
#include <fstream>      // Writing to file.
#include <iostream>     // Input-output stream.
#include "config.h"     // Configuration definitions.

// Function to output the configuration to file.
void WriteConfig()
{
// Delcare the file stream.
    std::fstream file;

    // Open the file for writing.
    file.open(static_cast<std::string>(RVF_ODIR) + "config.csv", std::fstream::out);

    // Check that the file opened.
    if (file.is_open())
    {
        // If the file could be opened, write the header.
        file << "CONFIG_NAME,CONFIG_VALUE";
        
        // For each configuration, record its setting.
        file << "\nRVF_IDIR," << RVF_IDIR; 
        file << "\nRVF_POSTDIR," << RVF_POSTDIR; 
        file << "\nRVF_ODIR," << RVF_ODIR; 
        file << "\nRVF_UPDATE_FREQ," << RVF_UPDATE_FREQ; 
        file << "\nRVF_NSAMPLES," << RVF_NSAMPLES; 
        file << "\nRVF_POST_NSAMPLES," << RVF_POST_NSAMPLES; 
        file << "\nRVF_MIN_NDVI," << RVF_MIN_NDVI; 
        file << "\nRVF_TRANSMISSION," << RVF_TRANSMISSION; 
        file << "\nRVF_DIFF_SCALE," << RVF_DIFF_SCALE; 
        file << "\nRVF_DIFF_NDVI," << RVF_DIFF_NDVI; 
        file << "\nRVF_MAX_BUFFER," << RVF_MAX_BUFFER; 
        file << "\nRVF_WRITE_MIN_PRECISION," << RVF_WRITE_MIN_PRECISION; 
        file << "\nRVF_WRITE_FULL," << RVF_WRITE_FULL; 
        file << "\nRVF_WRITE_SUMMARY," << RVF_WRITE_SUMMARY; 
        file << "\nRVF_SUMMARY_SIZE," << RVF_SUMMARY_SIZE; 
        file << "\nRVF_NO_VAC," << RVF_NO_VAC; 
        file << "\nRVF_LONG_TERM," << RVF_LONG_TERM; 
        file << "\nRVF_OPTIMISE," << RVF_OPTIMISE; 
        file << "\nRVF_OPTIM_RUNS," << RVF_OPTIM_RUNS; 
        file << "\nRVF_OPTIM_NCOST," << RVF_OPTIM_NCOST; 
        file << "\nRVF_OPTIM_NPARTICLES," << RVF_OPTIM_NPARTICLES; 
        file << "\nRVF_OPTIM_MAX_STEPS," << RVF_OPTIM_MAX_STEPS; 
        file << "\nRVF_OPTIM_SUMMARY," << RVF_OPTIM_SUMMARY; 
        file << "\nRVF_OPTIM_NDELTA," << RVF_OPTIM_NDELTA; 
        file << "\nRVF_OPTIM_MAX_MH," << RVF_OPTIM_MAX_MH; 
        file << "\nRVF_OPTIM_FUNC," << RVF_OPTIM_FUNC; 
        file << "\nRVF_OPTIM_POWER," << RVF_OPTIM_POWER;

        // Close the file once completed.
        file.close();
    }
    else
    {
        std::cout << "Configuration file could not be opened for writing." << std::endl;
        exit(1);
    }    
}