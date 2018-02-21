# Predicting the fate of the Tasmanian devil with a spatial epidemiological model

Spatially explicit, individual-based framework to model the spread of DFTD in Tasmania, implemented in Matlab.

## Prerequisites

You will need to have `Matlab` installed on your computer to run the collection of scripts. A single simulation runs on `Octave`, but `Sampler.m` needs to be changed to accommodate for a difference in the random number generator command.

## Pipeline

### Running single simulation

`MainFit.m` runs simulation for 25 years, as used for fitting. `MainLong.m` runs simulation for 250 years, as used to investigate long-term behaviour. `RunOne.m` runs a single simulation for 25 years, using the best-fitting parameter sets.

### MCMC parameter sweep

`Sampler.m` performs a sweep, outputting summary statistics into a text file.

### Rerun best 200 parameters

`RunLongTerms.m` runs a subset of the best 200 parameters and stores detailed information in a binary m-file.

## Authors

* **Veronika Siska** - *Technical work* - [siskavera](https://github.com/siskavera)
* **Andrea Manica** - *Supervisor* 
* **Bernhard Mehlig** - *Supervisor*
* **Anders Eriksson** - *Supervisor* 
