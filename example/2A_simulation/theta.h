// theta.h: Defines classes of a parameter to be inferred with different priors.
#include <string>           // String standard library.
#include "gsl/gsl_rng.h"    // Defines random number generator.

#ifndef THETA_H
#define THETA_H

// Base class of an inferred parameter.
class Theta
{
public:
    std::string par_name;                                   // Name of the parameter being inferred.
    double* chain;                                          // Chain of parameter values.
    virtual void Sample(int);                               // Sample from the underlying parameter distribution.
    virtual void ProjectSample(int);                        // If a sample in the chain is outside the support of the distribution, project
                                                            // it back into the support of the distribution.
    virtual void CompleteSample(int);                       // Function used to complete a sample if only the free parameters are computed. 
    virtual bool IsValidSample(int);                        // Function to check whether a given sample in the chain is within the support of the dist.                            
    virtual std::string OutputString(int, char sep = ',');  // Construct a string out of a given sample ready for output.
    virtual double Pdf(int);                                // The PDF of the distribution at a position in the chain.
    double GetPar(int idx);                                 // Retrieve a single parameter of a given distribution.
    int GetNumFreeDims();                                   // Retrieve the number of free parameters that are sampled from the distribution.
    int GetNumDims();                                       // Retrieve the number of dimensions of the distribution.
    Theta(std::string in_par_name);                         // Default constructor.
    virtual ~Theta();                                       // Default deconstructor.

protected:
    std::string dist_name;  // Name of the parameter distribution.
    int n_pars;             // Number of parameters that define the distribution.
    int n_free_dims;        // Number of free parameters that are sampled from the distribution.
    int n_dims;             // Number of dimensions in the distribution.
    double* pars;           // Parameter values of the distribution.
};

// Class of an scenario parameter with an uniform distribution.
class Uniform : public Theta
{
public:
    // Function to sample from a uniform distribution.
    void Sample(int scenario_id);

    // Function used to project a sample from outside the 
    // support of the distribution back into the support.
    void ProjectSample(int);

    // Function to check if sample is within distribution's support.
    bool IsValidSample(int);

    // Function to arrange a given scenario as an output string.
    std::string OutputString(int scenario_id, char sep = ',');

    // The PDF of a uniform distribution.
    double Pdf(int);

    // Get a parameter of the uniform distribution.
    double GetPar(int idx);

    // Constructor of an inferred parameter with a uniform distribution.
    Uniform(std::string in_par_name, double in_a, double in_b, int in_seed);

    // Deconstructor of uniform distribution.
    ~Uniform();

private:
    // Random number generator for sampling from the distribution.
    gsl_rng* rng;
};

class Dirichlet : public Theta
{
public:
    // Function to sample from a uniform distribution.
    void Sample(int scenario_id);

    // Function used to project a sample from outside the 
    // support of the distribution back into the support.
    void ProjectSample(int);

    // Function used to complete a sample if only the free parameters
    // are computed.
    void CompleteSample(int);

    // Function to check if sample is within distribution's support.
    bool IsValidSample(int);

    // The PDF of a Dirichlet distribution.
    double Pdf(int);

    // Function to arrange a given scenario as an output string.
    std::string OutputString(int scenario_id, char sep = ',');

    // Constructor of an inferred parameter with a uniform distribution.
    Dirichlet(std::string in_par_name, int in_k, int in_seed);

    // Deconstructor of the Dirichlet distribution.
    ~Dirichlet();

private:

    // Random number generator for sampling from the distribution.
    gsl_rng* rng;

    // Vector to store sorted samples. Used in reprojection.
    double* sorted_desc;
};

class DiscreteUniform : public Theta
{
public:
    // Function to sample from a uniform distribution.
    void Sample(int scenario_id);

    // Function used to project a sample from outside the 
    // support of the distribution back into the support.
    void ProjectSample(int);

    // Function used to complete a sample if only the free parameters
    // are computed.
    void CompleteSample(int);

    // Function to check if sample is within distribution's support.
    bool IsValidSample(int);

    // The PDF of a Discrete Uniform distribution.
    double Pdf(int);

    // Function to arrange a given scenario as an output string.
    std::string OutputString(int scenario_id, char sep = ',');

    // Constructor of an inferred parameter with a uniform distribution.
    DiscreteUniform(std::string in_par_name, int in_a, int in_b, int in_seed);

    // Deconstructor for the discrete uniform distribution.
    ~DiscreteUniform();

private:

    // Random number generator for sampling from the distribution.
    gsl_rng* rng;
};

// Create a distriubution that gets its values from an existing
// distribution of input values.
// This input file must contain all relevant parameters for the 
// name of the distribution.
class Empirical : public Theta
{  
public:
    // Function to sample from the empirical distribution.
    void Sample(int);

    // Function to arrange a given scenario as an output string.
    // Separate dimensions of a single sample appear on separate lines.
    std::string OutputString(int scenario_id, char sep = ',');

    // Constructor of an empirical distribution of parameters.
    Empirical(std::string in_par_name, std::string in_file_name, int in_seed);

    // Deconstructor for the empirical distribution.
    ~Empirical();

private:
    // Function to read in the empirical distribution of parameters.
    // Note that this function discards any parameters which do not
    // have a match with the requested parameter name.
    void ReadEmpiricalFile();

    // Define space to store the empirical samples.
    double* samples;

    // The file-name containing the empirical distribution samples.
    std::string file_name;

    // Define how many empirical samples there are.
    int n_samples;

    // Random number generator for sampling from the distribution.
    gsl_rng* rng;
};

#endif  