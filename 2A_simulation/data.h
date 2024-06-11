// data.h: Defines the Seroprev and Data classes. These class contains the data
// needed to run and fit the RVF simulation model.
#include "parameters.h" // Definition of the Parameter class.

#ifndef DATA_H
#define DATA_H

class Ndvi
{
public:
    int n_island;                       // The number of islands.
    int n_time_steps;                   // Number of time steps/entries in the ndvi data.
    double start_month;                 // Starting month of the NDVI data.
    double start_year;                  // Starting year of the NDVI data.
    double** ndvi;                      // Time-series of the NDVI data.
    double* min_local_ndvi;             // Get the minimum NDVI on each island.
    double min_global_ndvi;             // Get the minimum NDVI on all islands.
    Ndvi(Parameters*);                  // Constructor for the Ndvi class.
    Ndvi(const Ndvi&);                  // Deep-copy constructor.
    ~Ndvi();                            // Deconstructor.
private:
    void ReadNDVI(Parameters*);         // Function to read and store the NDVI data.
                                        // Also sets the simulation length.
};

class Data
{
public:
    double start_month;                 // Starting month of the data.
    double start_year;                  // Starting year of the data.
    Ndvi ndvi;                          // NDVI data.
    Data(Parameters*);                  // Constructor.
    Data(const Data&);                  // Constructor for a deep copy implementation.
private:

};

#endif