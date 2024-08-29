# Supporting data and code for the paper entitled *Effectiveness and equity of vaccination strategies against Rift Valley fever in a heterogeneous landscape*
This repository contains the data and code which were used to produce the manuscript entitled '[Effectiveness and equity of vaccination strategies against Rift Valley fever in a heterogeneous landscape](doi.org/10.1101/2024.07.18.604096)'.

In our manuscript, we developed a metapopulation model for RVF transmission in livestock across the Comoros archipelago which incorporates livestock vaccination in addition to heterogeneity in viral transmission rates and animal movements. We used the model to evaluate three vaccine allocation strategies—proportional allocation, optimal allocation for maximising total infections averted across the archipelago, and optimal allocation for more equitable outcomes across islands—under different vaccination coverage levels and animal identification scenarios.

Below we outline how to generate the outputs of one vaccine allocation strategy using the supplied example code, a descipription of the input data used to simulate simulate the mathematical model forward in time, and provide a description of the full source code used to produce all results in our manuscript: available as a pre-print on bioRxiv at [doi.org/10.1101/2024.07.18.604096](doi.org/10.1101/2024.07.18.604096).

### Table of Contents
- [Example code for reproducing results from the study](#example-code-for-reproducing-results-from-the-study)
    - [Step 1A: optimisation](#step-1a-optimisation)
    - [Step 1B: format optimisation output](#step-1b-format-optimisation-output)
    - [Step 2A: model simulation](#step-2a-model-simulation)
    - [Step 2B: format model simulation output](#step-2b-format-model-simulation-output)
    - [Step 3: visualisation](#step-3-visualisation)
 - [Description of the input data](#description-of-the-input-data)
    - [Normalised Difference Vegetation Index](#normalised-difference-vegetation-index)
    - [Posterior of model parameters](#posterior-of-model-parameters)
 - [Description of the source code](#description-of-the-source-code)
    - [Optimisation and model simulation](#optimisation-and-model-simulation)
    - [Data manipulation and visualisation](#data-manipulation-and-visualisation)

## Example code for reproducing results from the study
In the `example` folder is a set of example code that when executed sequentially will reproduce a (small) subset of the main results in the paper.

This example code will reproduce results for the vaccine strategy where 5%, 10%, 15%, 20%, 25% and 30% of animals are vaccinated annually across the archipelago, vaccines are distributed optimally to maximise the infections averted across thie archipelago, and it is assumed that animals are not tagged post-vaccination. This is referred to as strategy `U1` in the manuscript.

The input files necessary to run each step are already in each subfolder, so the steps may be skipped if desired. Note that the execution of all example code in this section took 2.5 hours on one CPU thread, and 12.5 minutes on twelve CPU threads on an AMD Ryzen 9 5900X processor. This runtime does not include time to read the documentation or install any pre-requisites.

### <ins>Step 1A: optimisation</ins>
The purpose of this code is to use the Sequential Monte Carlo optimisation algorithm to determine the optimal way to distribute vaccines between the islands of the Comoros archipelago. This code is contained in the `example/1A_optimisation` folder.

#### Pre-requisites
In order to run this example, your system requires a method for compiling `C++` code. Here, we used the gcc version 13.1.0 executed on Windows 11. The environment also requires a compiled library of GNU Scientific Library (GSL). The version used in this code is GSL 2.7.

#### Execution
- Prior to compilation, the number of CPU threads used to execute the code may be changed by modifying `RVF_NTHREADS` definition in the `config.h` file. The default number of threads is 1.
- Set the working directory to `example/1A_optimisation`.
- Compile all .cpp files in the `example/1A_optimisation` folder using
```cpp
    g++.exe *.cpp -Wall -Wextra -fopenmp -O3 -std=c++20 -lgsl -lgslcblas -o rvf_comoros_vaccine.exe
```
- Execute `rvf_comoros_vaccine.exe`.
Executing the code on an AMD Ryzen 9 5900X processor took 2 hours on 1 thread, or 10 minutes on 12 threads.

#### Expected output
Compiling the code and running the executable results in a series of .csv files that contain the cost and vaccine distribution of all particles in multiple runs of the optimisation algorithm.

#### Technical notes
For computational reasons, the default number of particles used in the Sequential Monte Carlo algorithm was reduced for this example code.

### <ins>Step 1B: format optimisation output</ins>
Executing the above code results in a series of a particles representing different vaccine distributions and their resultant effectiveness. The output of executing the above code needs to be formatted in order to extract the optimal particle per run of the optimisation run. This code is contained in the `example/1B_format` folder.

#### Pre-requisites
The .csv files from the above code may be copy and pasted into the `in` folder, or the supplied example .csv files may be used. In order to run this example, your system requires a method for executing `R` code. Here, we used `R` version 4.3.2 on Windows 11. The following `R` packages need to be installed prior to executing the code: data.table (version 1.14.10), tidyr (version 1.3.0) and dplyr (version 1.1.4).

#### Execution
- Set the working directory to `example/1B_optimisation`.
- Run `data_format_optimal.R` in an `R` environment.
Executing the code on an AMD Ryzen 9 5900X processor took less than 1 second.

#### Expected output
Executing this code will produce a series of .csv files containing vaccine distributions to achieve optimal strategy effectiveness in the `out` folder.

### <ins>Step 2A: model simulation</ins>
The purpose of this code is to simulate the mathematical model using the optimal vaccine distributions obtained from the above steps. This code is contained in the `example/2A_simulation` folder.

#### Pre-requisites
An example output of the previous steps have been placed in the `in` folder for convenience. The data produced by the previous steps may be used instead, but must be placed in the `in` folder. As with Step 1A, your system requires a method for compiling `C++` code. Here, we used the gcc version 13.1.0 executed on Windows 11. The environment also requires a compiled library of GNU Scientific Library (GSL). The version used in this code is GSL 2.7.

#### Execution
- Prior to compilation, the number of CPU threads used to execute the code may be changed by modifying `RVF_NTHREADS` definition in the `config.h` file. The default number of threads is 1.
- Set the working directory to `example/2A_simulation`.
- Compile all .cpp files in the `2A_simulation` folder using
```cpp
    g++.exe *.cpp -Wall -Wextra -fopenmp -O3 -std=c++20 -lgsl -lgslcblas -o rvf_comoros_vaccine.exe
```
- Run the executable `rvf_comoros_vaccine.exe`.
Executing the code on an AMD Ryzen 9 5900X processor took 24 minutes on 1 thread, or 2 minutes on 12 threads.

#### Expected output
The output is a series of .csv files containing the number of individuals in each model compartment (see paper for details) at each time step in the `out` folder. 

### <ins>Step 2B: format model simulation output</ins>
The above code will output the number of individuals in each model compartment over time. This makes it more convenient to produce visualisations of the data without having to reload the entire data set each time. This code is contained in the `example/2B_summarise` folder.

#### Pre-requisites
The .csv files from the above code may be copy and pasted into the `in` folder, or the supplied example .csv files may be used. In order to run this example, your system requires a method for executing `R` code. Here, we used `R` version 4.3.2 on Windows 11. The following `R` packages need to be installed prior to executing the code: data.table (version 1.14.10), tidyr (version 1.3.0) and dplyr (version 1.1.4).

#### Execution
- Set the working directory to `example/2B_summarise`.
- Run `data_summarise_scenario.R` in an `R` environment.
Executing the code on an AMD Ryzen 9 5900X processor took 5 seconds.

#### Expected output
The output is a series of .csv files, which summarise the model simulations produced in step 2A, located in the `out` folder.

### <ins>Step 3: visualisation</ins>
The optimisation and simulation outputs above are used to produce the figures and summary metrics presented in the paper. This code is contained in the `example/3_visualisation` folder.

#### Pre-requisites
The .csv files from the above code may be copy and pasted into the `in` folder, or the supplied example .csv files may be used. In order to run this example, your system requires a method for executing `R` code. Here, we used `R` version 4.3.2 on Windows 11. The following `R` packages need to be installed prior to executing the code: data.table (version 1.14.10), tidyr (version 1.3.0),  dplyr (version 1.1.4), egg (version 0.4.5), gridExtra (version 2.3), ggokabeito (version 0.1.0), colorspace (version 2.1), ggplot2 (version 3.4.4) and stringr (version 1.5.1).

#### Execution
- Set the working directory to `example/3_visualisation`.
- Run `plot_compare_opt_dist.R` in an `R` environment.
- Run `plot_compare_opt_effect.R` in an `R` environment.
Executing the code on an AMD Ryzen 9 5900X processor took 5 seconds.

#### Expected output
The output is a .pdf and .svg file of (subsets of) Figures 2 and Figures 3 from the paper: the percentage share of vaccines between each island, and the optimal percentage of infections averted across the archipelago, when vaccinating 5%, 10%, 15%, 20%, 25% and 30% of the livestock annually, where livestock are not tagged.

## Description of the input data
Two data sets were used as input to the mathematical model and optimisation algorithm used in the manuscript: Normalised Difference Vegetation Index over time on each island in the Comoros archipelago, and the posterior of model parameters obtained from fitting the model to serological data—see our previous work ([Tennant, WSD, et al. [2021]](https://doi.org/10.1038/s41467-021-25833-8)) and the associated repository ([wtennant/rvf_comoros](https://github.com/wtennant/rvf_comoros)) for details of the model fitting procedure.

### Normalised Difference Vegetation Index
Normalised Difference Vegetation Index (NDVI) was used to inform time-dependent transmission rates of RVF between livestock on each island in our transmission model. The data is provided as a `.csv` file in `data/ndvi_extended.csv`. This file gives the summarised NDVI data for each island in the Comoros archipelago and simulation week (note: there are 4 simulation weeks per month for ease of computation). The NDVI data from 1st July 2004 until 30th June 2020 is recycled to produce NDVI data up to 1st July 2050. The columns of this data are as follows:

| Column name | Description |
| ----------- | ----------- |
| ISLAND_ID | Unique identifier for each island in the archipelago: Anjouan (0), Grande Comore (1), Mayotte (2), Mohéli (3) |
| DATE | Approximate corresponding date for NDVI data |
| WEEK | Simulation week of a month |
| MONTH | Month of the calendar year |
| YEAR | Calendar year |
| NDVI | Value of NDVI at the given location and time point |

### Posterior of model parameters
In our previous work ([Tennant, WSD, et al. [2021]](https://doi.org/10.1038/s41467-021-25833-8)), we fitted an epidemiological model to Rift Valley fever surveillance data conducted in livestock across the Comoros archipelago between 2004 and 2015 in a Bayesian framework. Fitting the model to these data estimated the demographic and infection process parameters of the model: the probability of moving between islands per week, the island-specific disease transmission rate parameters, external importation of livestock, and initial proportion of the population immune to RVFV on each island.

We incorporated these estimated parameters into our assessment of vaccine strategy effectiveness by simulating from the model with different samples from the posterior distribution for each vaccine strategy tested. The posterior samples used in the analysis of our manuscript are provided as a `.csv` file in `data/posterior.csv`. The columns of the data are as follows:

| Column name | Description |
| ----------- | ----------- |
| ID | Unique identifier for each sample from the joint posterior distribution |
| PAR_NAME | Name of the model parameter in the posterior sample as it appears in the simulation code |
| PAR_VALUE | Value of the model parameter in the posterior sample |

## Description of the source code
To assess the effectiveness of different vaccine strategies against Rift Valley Fever, we produced custom `C++` and `R` code to optimise the allocation of vaccines across the Comoros archipelago, simulate from the mathematical model and visualise outputs of the optimisation algorithm and model simulation.

The `src` folder contains the `C++` and `R` code that was used to produce all results in the manuscript. Note that the code would need to be modified to produce each result presented in the main manuscript. The full set of results generated in the manuscript took 3 to 4 weeks on an AMD Ryzen 9 5900X processor. For this reason, a [waltkhrough example](#example-code-for-reproducing-results-from-the-study) is provided above to run a subset of the results from the main manuscript.

Below we describe the main files in the `C++` and `R` code which were used to produce the main results in the manuscript.

### Optimisation and model simulation
The `src/cpp` folder contains the source code to run the Sequential Monte Carlo optimisation algorithm and simulate from the mathematical model.

#### Pre-requisites
To produce executables from this code, your system requires a method for compiling `C++` code. Here, we used the gcc version 13.1.0 executed on Windows 11. The environment also requires a compiled library of GNU Scientific Library (GSL). The version used in this code is GSL 2.7.

#### Key files
All provided `.cpp` and `.h` files are fully commented and describe what they do. In order to reproduce the results from the study, or test different vaccine strategy scenarios using the model, the following key files would need to be modified:
- `config.h`: this file contains all configuration macros specifying technical details of the optimisation or model simulation. These include:
    - `RVF_OPTIMISE`: whether or not to perform an optimisation or model simulation,
    - `RVF_OPTIM_RUNS`: the number of runs of the optimisation algorithm to perform,
    - `RVF_NTHREADS`: the number of CPU threads to use in optimatisation or model simulation,
    - `RVF_OPTIM_FUNC`: which effectiveness function to optimise the vaccine allocation over,
    - `RVF_WRITE_FULL`: whether or not to write the full simulation outputs to file (these can consume an **extremely** large amount of space depending on how many model simulations are executed).
    - The definition of all other macros described as comments in the file.
- `main.cpp`: this file contains the calls to the optimisation or model simulation functions and under which vaccine strategy they should be performed.
    - For example, in this repository the `RunOptimisation()` function on `lines 106-126` is set up to run the optimisation algorithm for vaccination rates of 5% to 30%.
    - Similarly, in this repository the `RunScenarios()` function is set to read optimal vaccine allocation sets (as prepared by `src/r/data_format_optimal.R` below) and simulate from the mathematical model under scenarios where only young animals are vaccinated and/or vaccination only occurs at the start of each epidemiological year.
- `parameters.cpp`: this file contains the default parameter set used in the main manuscript. The parameters in this file can be modified to generate outputs similar to our study for a different set of intitial parameters. Similarly, the default parameters may be overriden elsewhere in the code (e.g. the vaccination rate) to produce the results from the manuscript.

#### Compilation and execution
After setting the working directory to `src/cpp`, we compiled the files to an executable using the following command:
```
    g++.exe *.cpp -Wall -Wextra -fopenmp -O3 -std=c++20 -lgsl -lgslcblas -o rvf_comoros_vaccine.exe
```
The code was then executed by calling `rvf_comoros_vaccine.exe`. Note that for every vaccine strategy tested and modification to the key files listed above, the code would need to be recompiled.

### Data manipulation and visualisation
The `src/r` folder contains the source code to manipulate outputs from the optimisation algorithm or model simulation, as well as code to produce all visualisations presented in the manuscript.

#### Pre-requisites
To run the scripts, your system requires a method for executing `R` code. Here, we used `R` version 4.3.2 on Windows 11. The following `R` packages need to be installed prior to executing the code: data.table (version 1.14.10), tidyr (version 1.3.0),  dplyr (version 1.1.4), egg (version 0.4.5), gridExtra (version 2.3), ggokabeito (version 0.1.0), colorspace (version 2.1), ggplot2 (version 3.4.4) and stringr (version 1.5.1). One script also requires shape files (`.shp`) for all four islands in the Comoros archipelago. These are not supplied in this repository.

#### Key files
All provided `.R` files are fully commented and describe what they do. Note that any directories set within these scripts are specific to the location of outputs produced by the the optimisation algorithm and model simulation code.

Here is a brief description of the key files used to reproduce our analysis:
- `data_format_optimal.R`: after executing the optimisation algorithm, this script produces a `.csv` file containing the optimal vaccine allocation across islands in the Comoros archipelago for the given scenario.
- `data_summarise_scenario_time.R`: after executing the model simulation, this script produces a `.csv` file containing a summary of vaccine efficiency and number of infecteds over time on each island in the Comoros archipelago.
- `data_summarise_scenario.R`: after executing the model simulation, this script produces a `.csv` file containing a summary of the number of infections with and without vaccination, and the percentage of infections averted over the simulation time period.
- `plot_compare_opt_dist.R`: after executing the optimisation algorithm and `data_format_optimal.R` for all vaccine strategies tested, this script produces a figure of the vaccine share between islands in the Comoros archipelago. This is equivalent to Figure 2 from the main text.
- `plot_compare_opt_effect.R`: after executing the optimisation algorithm, `data_format_optimal.R`, the model simulation and `data_summarise_scenario.R` for all vaccine strategies tested, this script produces a figure of the percentage of infections averted across the Comoros archipelago under each scenario. This is equivalent to Figure 3 in the main text.
- `plot_compare_opt_island_balance.R`: after executing the optimisation algorithm, `data_format_optimal.R`, the model simulation and `data_summarise_scenario.R` for all vaccine strategies tested, this script produces a figure of the percentage of infections averted on each island in the Comoros archipelago under each scenario. This is equivalent to Figure 4 in the main text. Note that this script will only execute successfully if shape files are also provided (not included in this repository).
- `plot_compare_opt_vac_efficiency.R`: after executing the optimisation algorithm, `data_format_optimal.R`, the model simulation and `data_summarise_scenario_time.R` for all vaccine strategies tested, this script produces a figure of the percentage of vaccines that are administered to susceptible livestock over time on each island in the Comoros archipelago under each scenario. This is equivalent to Figure 5 in the main text.
