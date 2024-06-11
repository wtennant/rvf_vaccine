// theta.cpp: Defines class constructors and member functions of
// the posterior distribution for each parameter.
#include <algorithm>            // Iterate through strings to erase characters.
#include <chrono>               // Generating a seed for a random number generator.
#include <cmath>                // Rounding numbers.
#include <iostream>             // Output to console.
#include <fstream>              // Reading in empirical distributions.
#include <string>               // Strings!
#include "gsl/gsl_randist.h"    // Random number distributions.
#include "gsl/gsl_sort.h"       // Sorting vectors.
#include "config.h"             // Number of scenarios to run.
#include "theta.h"              // Class defining parameters to be inferred.

// Constructor of the base class Theta.
Theta::Theta(std::string in_par_name)
{
    // Define the parameter name.
    par_name = in_par_name;
}

// Deconstructor for the base class Theta.
Theta::~Theta()
{
    // There is nothing to do.
}

// Sampling from a parameter with no underlying distribution.
void Theta::Sample(int scenario_id)
{
    std::cout << "If you see this message, then something has gone wrong with scenario " << scenario_id << "." << std::endl;
    std::cout << "You should not be sampling from a parameter with no underlying distribution!" << std::endl;
}

// Projecting a sample from a parameter with no underlying distribution.
void Theta::ProjectSample(int scenario_id)
{
    std::cout << "If you see this message, then something has gone wrong with scenario " << scenario_id << "." << std::endl;
    std::cout << "You should not be sampling from a parameter with no underlying distribution!" << std::endl;
}

// Completing a sample with no underlying distribution.
void Theta::CompleteSample(int particle_id)
{
    std::cout << "If you see this message, then something has gone wrong with particle " << particle_id << "." << std::endl;
    std::cout << "You should not be sampling from a parameter with no underlying distribution!" << std::endl;
}

// Checking if a sample is valid with no underlying distribution.
bool Theta::IsValidSample(int particle_id)
{
    std::cout << "If you see this message, then something has gone wrong with particle " << particle_id << "." << std::endl;
    std::cout << "You should not be sampling from a parameter with no underlying distribution!" << std::endl;
    return false;
}

// Getting a PDF with no underlying distribution.
double Theta::Pdf(int particle_id)
{
    std::cout << "If you see this message, then something has gone wrong with particle " << particle_id << "." << std::endl;
    std::cout << "You should not be sampling from a parameter with no underlying distribution!" << std::endl;
    return 0.0;
}

// Sampling from a parameter with no underlying distribution.
std::string Theta::OutputString(int scenario_id, char sep)
{
    std::cout << "If you see this message, then something has gone wrong with scenario " << scenario_id << "." << std::endl;
    std::cout << "You should not be outputting parameters to file that have no underlying distribution!" << std::endl;
    std::string out_str = "" + sep + sep + sep;
    return out_str;
}

// Retrieve a single parameter that specifies the underyling distribution.
double Theta::GetPar(int idx)
{
    // If the specified parameter to retrieve is out of bounds. Default to retrieving the first parameter.
    if (idx >= n_pars)
    {
        std::cout << "Index for parameter " << par_name << " is out of bounds. Retrieving first parameter." << std::endl;
        idx = n_pars;
    }
    return pars[idx];
}

// Function for number of free parameters
// that are in the support of the distribution.
int Theta::GetNumFreeDims()
{
    return n_free_dims;
}

// Function to return the dimension of the distribution.
int Theta::GetNumDims()
{
    return n_dims;
}

// Constructor of a parameter with a uniform distribution.
Uniform::Uniform(std::string in_par_name, double in_a, double in_b, int inc_seed) : Theta(in_par_name)
{
    // Define space for the sampled values.
    switch (RVF_OPTIMISE)
    {
        case 0:
            chain = new double[RVF_NSAMPLES];
            break;
        case 1:
            chain = new double[RVF_OPTIM_NPARTICLES];
            break;
        default:
            chain = new double[RVF_NSAMPLES];
            break;
    }

    // Define the number of parameters that define the distribution.
    n_pars = 2;

    // The number of parameters returned when sampling from the distribution is one.
    n_free_dims = 1;
    n_dims = 1;

    // Set up the paramters of a uniform distribution.
    pars = new double[n_pars]{in_a, in_b};

    // Set up a random number generator used to sample from the distribution.
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    long long seed = std::chrono::duration_cast<std::chrono::nanoseconds>
        (std::chrono::system_clock::now().time_since_epoch()).count();
    gsl_rng_set(rng, seed + inc_seed);

    // Set up distribution name.
    dist_name = "uniform";
}

// Deconstructor for the uniform distribution.
Uniform::~Uniform()
{
    // De-allocate parameter and chain space.
    delete[] chain;
    delete[] pars;
    
    // Also de-allocate the random number generator.
    gsl_rng_free(rng);
}

// Function to sample from a uniform ditribution and update the list of sampled
// values.
void Uniform::Sample(int scenario_id)
{
    // Randomly sample and save in the chain of parameters.
    chain[scenario_id] = gsl_ran_flat(rng, pars[0], pars[1]);
}

// Function to reproject a sample back onto a uniform distribution.
void Uniform::ProjectSample(int scenario_id)
{
    // If less than the lower bound, project to the lower bound.
    // If greater than the upper bound, project to the upper bound.
    if (chain[scenario_id] < pars[0])
    {
        chain[scenario_id] = pars[0];
    }
    else if (chain[scenario_id] > pars[1])
    {
        chain[scenario_id] = pars[1];
    }
}

// Function to check if a given sample form the distribution is valid,
// i.e. within the support of the distribution.
bool Uniform::IsValidSample(int particle_id)
{
    // There's only one thing to check..
    // If the sample is within the limits of the distribution.
    bool is_valid = (chain[particle_id] >= pars[0]) && (chain[particle_id] <= pars[1]);
    return is_valid;
}

// Function to convert a given sample to a comma-delimited string.
std::string Uniform::OutputString(int scenario_id, char sep)
{
    // Define a blank string.
    std::string out_str = "";

    // There is only one sampled value per scenario_id in a uniform distribution.
    out_str += std::to_string(scenario_id) + sep + par_name + sep + "0" + sep + std::to_string(chain[scenario_id]);
    
    // Return the output string ready for writing to file.
    return out_str;
}

// Function to evaluate the probability density function of the
// distribution at a specified index in the dsitributions chain.
double Uniform::Pdf(int scenario_id)
{
    // Use GSLs already written PDFs.
    return gsl_ran_flat_pdf(chain[scenario_id], pars[0], pars[1]);
}

// Constructor of a parameter with a Dirichlet distribution.
Dirichlet::Dirichlet(std::string in_par_name, int in_k, int inc_seed) : Theta(in_par_name)
{
    // Check that the specification of the Dirichlet distribution is valid.
    if (in_k < 1)
    {
        std::cout << "Specification of Dirichlet distribution for parameter " << in_par_name;
        std::cout << " is invalid." << std::endl;
        std::exit(2); // THIS IS BAD PRACTICE, CONSIDER EXCEPTIONS.
    }

    // Define the number of parameters in the distribution.
    // There are k alpha parameters (weight of each variable).
    n_pars = in_k + 1;

    // Define the number of free parameters in a sample of the distribution.
    // Although there are in_k probabilities, one probability depends upon the others.
    // Thus there is in_k - 1 free parameters returned from the distribution.
    n_free_dims = in_k - 1;

    // The number of dimensions to the distribution.
    n_dims = in_k;

    // Define space for the sampled values. There are as many as the number of categories.
    switch (RVF_OPTIMISE)
    {
        case 0:
            chain = new double[RVF_NSAMPLES*in_k];
            break;
        case 1:
            chain = new double[RVF_OPTIM_NPARTICLES*in_k]; 
            break;
        default:
            chain = new double[RVF_NSAMPLES*in_k];
            break;
    }

    // Define the space for a sorted sample.
    sorted_desc = new double[in_k];

    // Set up the paramters of a Dirichlet distribution.
    pars = new double[n_pars];
    pars[0] = in_k;
    for (int i = 1; i < n_pars; ++i)
    {
        pars[i] = 1.0;
    }

    // Set up a random number generator used to sample from the distribution.
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    long long seed = std::chrono::duration_cast<std::chrono::nanoseconds>
        (std::chrono::system_clock::now().time_since_epoch()).count();
    gsl_rng_set(rng, seed + inc_seed);

    // Set up distribution name.
    dist_name = "Dirichlet";
}

// Deconstructor for the uniform distribution.
Dirichlet::~Dirichlet()
{
    // De-allocate parameter and chain space.
    delete[] chain;
    delete[] pars;
    
    // Also de-allocate the random number generator.
    gsl_rng_free(rng);

    // De-allocate sorted sample space.
    delete[] sorted_desc;
}

// Function to sample from a Dirichlet ditribution and update the list of sampled values.
void Dirichlet::Sample(int scenario_id)
{
    // Ensure there is more than one category, or we cannot sample.
    if (pars[0] > 1)
    {
        // Randomly sample and save in the chain of parameters.
        gsl_ran_dirichlet(rng, static_cast<int>(pars[0]), &pars[1], chain + scenario_id*static_cast<int>(pars[0]));
    }
    else
    {
        chain[scenario_id] = 1.0;
    }
}

// Function to reproject a sample outside the support of
// the dirichlet distribution back into the support of the
// distribution.
void Dirichlet::ProjectSample(int scenario_id)
{
    // Get the dimension of the distribution.
    int k = static_cast<int>(pars[0]);

    // First, ensure that the final probability
    // in the (k-1)-simplex has been calculated.
    // Initialise the k-th probability.
    double k_prob = 1.0;

    // Subtract off the values of the (k-1) parameters in the distribution.
    for (int i = 0; i < (k - 1); ++i)
    {
        k_prob -= chain[scenario_id*k + i];
    }

    // Place in the final probability into the chain.
    chain[(scenario_id + 1)*k - 1] = k_prob;

    // First sort the samples into descending order.
    gsl_sort_largest(sorted_desc, k, &chain[scenario_id*k], 1, k);

    // Find the largest index (1:k) which satisfies
    // (sum_j=1^idx sorted_desc[j-1] - 1 / idx) < sorted_desc[idx - 1].
    int max_idx = 1;
    for (int i = 1; i <= k; ++i)
    {
        // (sum_j=1^idx sorted_desc[j-1] - 1 / idx)
        double total = 0.0;
        for (int j = 1; j <= i; ++j){ total += sorted_desc[j - 1]; }
        total = (total - 1) / i;

        // If the comparison is satisfied, then the current i will be
        // the largest index so far.
        if (sorted_desc[i - 1] > total){ max_idx = i; }
    }

    // Calculate the value of the above sum at the largest index value.
    double tau = 0.0;
    for (int j = 1; j <= max_idx; ++j){ tau += sorted_desc[j - 1]; }
    tau = (tau - 1) / max_idx;

    // Use this constant to re-scale the original set of samples.
    // Make sure that the values are non-negative.
    for (int i = 0; i < k; ++i)
    {
        chain[scenario_id*k + i] = (chain[scenario_id*k + i] - tau < 0.0) ? 0.0 : chain[scenario_id*k + i] - tau; 
    }
}

// Function to complete the sample if all other free parameters are computed.
// The final probability can be computed from all other probabilities.
void Dirichlet::CompleteSample(int particle_id)
{
    // Get the number of parameters in the distribution.
    int k = pars[0];

    // Add up the first (k-1) values for the particle sample.
    double sample_sum = 0.0;
    for (int i = 0; i < k - 1; ++i)
    {
        sample_sum += chain[particle_id*k + i];
    }
    
    // Complete the sample by subtracting the sum of the other values
    // from one. If the sum of all values is less than 1, and all
    // values are positive, then the sample is valid.
    chain[particle_id*k + k - 1] = 1.0 - sample_sum;
}

// Function to check if a given sample form the distribution is valid,
// i.e. within the support of the distribution.
bool Dirichlet::IsValidSample(int particle_id)
{
    // Get the dimension of the distribution.
    int k = pars[0];

    // If the sum of all values (minus the last) is less than 1, and
    // are themselves positive, then the sample is valid.
    // Keep track of the sum and if the values are greater than one.
    bool is_valid = true;
    double sample_sum = 0.0;
    for (int i = 0; i < k - 1; ++i)
    {
        sample_sum += chain[particle_id*k + i];
        if (chain[particle_id*k + i] < 0.0) // If satisfied, then the sample is not valid.
        {
            is_valid = false;
        }
    }

    // Check if the sample sum is less than one, and 
    // that all individual values were true.
    is_valid = ((sample_sum < 1.0) && is_valid);
    return is_valid;
}

// Function to convert a given sample to a comma-delimited string.
std::string Dirichlet::OutputString(int scenario_id, char sep)
{
    // Define a blank string.
    std::string out_str = "";

    // For each category, add in the sampled integer.
    int k = static_cast<int>(pars[0]);
    for (int i = 0; i < k; ++i) 
    {
        out_str += ((i > 0) ? "\n" : "") + std::to_string(scenario_id) + sep + par_name + sep;
        out_str += std::to_string(i) + sep + std::to_string(chain[scenario_id*k + i]);
    }
    
    // Return the output string ready for writing to file.
    return out_str;
}

// Function to evaluate the probability density function of the
// distribution at a specified index in the dsitributions chain.
double Dirichlet::Pdf(int scenario_id)
{
    // Retrieve the dimension of the dirichlet distribution.
    int k = static_cast<int>(pars[0]);

    // Use GSLs already written PDFs.
    return gsl_ran_dirichlet_pdf(k, &(pars[1]), &(chain[scenario_id*k]));
}

// Constructor of a parameter with a uniform distribution.
DiscreteUniform::DiscreteUniform(std::string in_par_name, int in_a, int in_b, int in_seed) : Theta(in_par_name)
{
    // Define space for the sampled values.
    switch (RVF_OPTIMISE)
    {
        case 0:
            chain = new double[RVF_NSAMPLES];
            break;
        case 1:
            chain = new double[RVF_OPTIM_NPARTICLES];  // Two extra spaces: one for initial values, one for proposals.
            break;
        default:
            chain = new double[RVF_NSAMPLES];
            break;
    }

    // Define the number of parameters in the distribution.
    n_pars = 2;

    // When sampling from a discrete uniform distribution there is only a single parameter returned.
    n_free_dims = 1;
    n_dims = 1;

    // Set up the paramters of a uniform distribution.
    pars = new double[n_pars]{static_cast<double>(in_a), static_cast<double>(in_b)};

    // Set up a random number generator used to sample from the distribution.
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    long long seed = std::chrono::duration_cast<std::chrono::nanoseconds>
        (std::chrono::system_clock::now().time_since_epoch()).count();
    gsl_rng_set(rng, seed + in_seed);

    // Set up distribution name.
    dist_name = "discrete_uniform";
}

// Deconstructor for the uniform distribution.
DiscreteUniform::~DiscreteUniform()
{
    // De-allocate parameter and chain space.
    delete[] chain;
    delete[] pars;
    
    // Also de-allocate the random number generator.
    gsl_rng_free(rng);
}

// Function to sample from a uniform ditribution and update the list of sampled
// values.
void DiscreteUniform::Sample(int scenario_id)
{
    // Randomly sample and save in the chain of parameters.
    chain[scenario_id] = pars[0] + static_cast<double>(gsl_rng_uniform_int(rng, static_cast<unsigned long>(pars[1]) - static_cast<unsigned long>(pars[0]) + 1));
}

// Function to reproject a sample back onto a discrete uniform distribution.
void DiscreteUniform::ProjectSample(int scenario_id)
{
    // If less than the lower bound, project to the lower bound.
    // If greater than the upper bound, project to the upper bound.
    // Round decimals to the nearest integer.
    if (chain[scenario_id] < pars[0])
    {
        chain[scenario_id] = pars[0];
    }
    else if (chain[scenario_id] > pars[1])
    {
        chain[scenario_id] = pars[1];
    }
    else
    {
        chain[scenario_id] = round(chain[scenario_id]);
    }
}

// Function to complete the sample if all other free parameters are computed.
void DiscreteUniform::CompleteSample(int particle_id)
{
    // Randomly round the sample if it is not an integer.
    // Get the fractional part of the sample.
    double frac = chain[particle_id] - floor(chain[particle_id]);

    // Round up with probability frac.
    if (gsl_ran_flat(rng, 0.0, 1.0) < frac)
    {
        chain[particle_id] = ceil(chain[particle_id]);
    }
    else
    {
        chain[particle_id] = floor(chain[particle_id]);
    }
}

// Function to check if a given sample form the distribution is valid,
// i.e. within the support of the distribution.
bool DiscreteUniform::IsValidSample(int particle_id)
{
    // There's only one thing to check..
    // If the sample is within the limits of the distribution.
    bool is_valid = (chain[particle_id] >= pars[0]) && (chain[particle_id] <= pars[1]);
    return is_valid;
}

// Function to convert a given sample to a comma-delimited string.
std::string DiscreteUniform::OutputString(int scenario_id, char sep)
{
    // Define a blank string.
    std::string out_str = "";

    // There is only one sampled value per scenario_id in a uniform distribution.
    out_str += std::to_string(scenario_id) + sep + par_name + sep + "0" + sep + std::to_string(chain[scenario_id]);
    
    // Return the output string ready for writing to file.
    return out_str;
}

// Function to evaluate the probability density function of the
// distribution at a specified index in the dsitributions chain.
double DiscreteUniform::Pdf(int scenario_id)
{
    // Use GSLs already written PDFs.
    return gsl_ran_flat_pdf(chain[scenario_id], pars[0], pars[1] + 1.0);
}

// Constructor for the empriical distribution of parameters.
Empirical::Empirical(std::string in_par_name, std::string in_file_name, int in_seed) : Theta(in_par_name)
{
    // Store the input file name containing the empirical distribution info.
    file_name = in_file_name;

    // Read in the empirical distribution from a pre-processed file.
    ReadEmpiricalFile();

    // Setup spaec to store samples from the empirical distirbution.
    switch(RVF_OPTIMISE)
    {
        case 0:
            chain = new double[RVF_NSAMPLES*n_dims]; break;
        case 1:
            chain = new double[RVF_OPTIM_NPARTICLES*n_dims]; break;
        default:
            chain = new double[RVF_NSAMPLES*n_dims]; break;
    }    
    
    // When sampling from an unknown distribution, we do not have
    // any parameters specifying the distribution except the
    // number of parameters per sample from the distribution.
    n_free_dims = n_dims;
    n_pars = 0;
    pars = nullptr;

    // Set up a random number generator used to sample from the distribution.
    rng = gsl_rng_alloc(gsl_rng_mt19937);
    long long seed = std::chrono::duration_cast<std::chrono::nanoseconds>
        (std::chrono::system_clock::now().time_since_epoch()).count();
    gsl_rng_set(rng, seed + in_seed);

    // Set up distribution name.
    dist_name = "empirical";
}

// Deconstructor for an empirical distribution.
Empirical::~Empirical()
{
    // De-allocate chain and sample space.
    delete[] chain;
    delete[] samples;
    
    // Also de-allocate the random number generator.
    gsl_rng_free(rng);
}

// Define how to sample from an empirical distribution.
void Empirical::Sample(int scenario_id)
{
    // Use the input samples as is: cycle between 0 to n_samples - 1
    // if there are repeats.
    int sample_idx = scenario_id % n_samples;

    // Set the sample in the chain of samples.
    for (int i = 0; i < n_dims; ++i)
    {
        chain[scenario_id*n_dims + i] = samples[sample_idx*n_dims + i];
    }
}

// Read in a file containing an empirical distribution of
// model parameters.
void Empirical::ReadEmpiricalFile()
{
    // Initialise the dimensions in the distribution.
    n_dims = 0;

    // Create a file stream.
    std::fstream file;
    
    // Open the file to read in.
    file.open(static_cast<std::string>(RVF_POSTDIR) + file_name,
              std::fstream::in);

    // Check that the file opened.
    if (file.is_open())
    {
        // Define the contents of a line.
        std::string line;

        // Throw away the first line of the data (the header) again.
        getline(file, line);
 
        // Define the number of colums in the data frame.
        // The column format will be ID, parameter name, index, and value.
        int n_cols = 4;

        // Define the separator of entries.
        char sep = ',';

        // Get the number of lines in the file that are
        // associated with the requested input parameter.
        int entries = 0;
        n_samples = 0;
        bool check_entry_idx = false;
        while (getline(file, line, sep))
        {
            // Note how many entries have been read.
            ++entries;

            // Strip out the \" either side of the line.
            line.erase(std::remove(line.begin(), line.end(), '\"'), line.end());

            // If there has been a request to check the index of the parameter,
            // then note down if this is the largest index seen so far.
            if (check_entry_idx)
            {
                // Retrieve the index number associated with the parameter.
                int idx = std::stoi(line);
                if (n_dims < (idx + 1))
                {
                    n_dims = idx + 1;
                }
                check_entry_idx = false;
            }

            // If the parameter name in the file matches the requested parameter, then
            // count how many lines match that parameter value.
            if (par_name == line)
            {
                ++n_samples;

                // Note that the index of the next entry should be checked.
                check_entry_idx = true;
            }

            // Define the next separator character.
            // If about to read the final entry of a line, change the separator
            // to a new line.
            sep = (entries + 1) % n_cols == 0 ? '\n' : ',';
        }

        // Divide the number of samples by the number of dimensions.
        n_samples = n_samples / n_dims;

        // For each parameter in the empirical distribution, allocate space
        // to store the samples. There are n_dims values in each sample.
        samples = new double[n_samples*n_dims];

        // Reset the position of the file stream to the first line.
        file.clear();
        file.seekg(0, std::ios::beg);

        // Throw away the first line of the data (the header) again.
        getline(file, line);

        // Keep track of the number of individual entries scanned, and the
        // number of posteriors samples read in.
        entries = 0;
        int sample_id = 0;
        int par_id = 0;

        // Reset the separator index.
        sep = ',';

        // For each line, extract the parameter name and value.
        while(getline(file, line, sep))
        {       
            // Increase the counter for the number of entries that have been read.
            ++entries;

            // If in the first column, note what the sample ID value is.
            if ((entries - 1) % n_cols == 0)
            {
                sample_id = std::stoi(line);
            }

            // Strip out any \" as part of the line.
            line.erase(std::remove(line.begin(), line.end(), '\"'), line.end());

            // Check that we are recording the correct parameter value.
            if ((entries - 2) % n_cols == 0) 
            {
                if (par_name == line)
                {
                    check_entry_idx = true;
                }
            }

            // If we are in the third column and we have the right parameter,
            // note the associated parameter index.
            if (((entries - 3) % n_cols == 0) && (check_entry_idx))
            {
                par_id = std::stoi(line);
            }

            // In the final column, we can then appropriately record the parameter value.
            if ((entries % n_cols == 0) && (check_entry_idx))
            {
                samples[sample_id*n_dims + par_id] = std::stod(line);
                check_entry_idx = false;
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
        std::cout << "Error: empirical data file could not be read." << std::endl;
        exit(1);
    }
}

// Function to convert a given sample to a comma-delimited string.
std::string Empirical::OutputString(int scenario_id, char sep)
{
    // Define a blank string.
    std::string out_str = "";

    // For each category, add in the sampled integer.
    int k = GetNumDims();
    for (int i = 0; i < GetNumDims(); ++i) 
    {
        out_str += ((i > 0) ? "\n" : "") + std::to_string(scenario_id) + sep + par_name + sep;
        out_str += std::to_string(i) + sep + std::to_string(chain[scenario_id*k + i]);
    }
    
    // Return the output string ready for writing to file.
    return out_str;
}