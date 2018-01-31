# Predicting the fate of the Tasmanian devil with a spatial epidemiological model

Spatially explicit, individual-based framework to model the spread of DFTD in Tasmania, implemented in Matlab.

## Metapopulation model

Our spatially explicit model consisted of an individual-based (thus stochastic) metapopulation of local populations, arranged on a square grid of 10km by 10km cells representing mainland Tasmania. The local populations on each grid were assumed to be well-mixed and only nearest-neighbour interactions (contact and migration) in the von Neumann neighbourhood (4 immediately neighbouring cells) were allowed. The carrying capacity of each cell was set using published data from [Hawkins et al., 2006](https://doi.org/10.1016/j.biocon.2006.04.010). The structure of the model is shown below.

![alt text](https://github.com/siskavera/tasmanian-devil/tree/master/images/model_structure.png "Model structure")

The local dynamics consisted of a compartmental epidemiological model built on demography following logistic growth. There were four disease-related compartments: susceptible (S), exposed (E, no symptoms or infectivity), diagnosable (D, show symptoms but not infective) and infective. The cancer was assumed to lead to certain death, since very few recoveries have ever been observed4. We used a frequency-dependent transmission rate, since it was shown to fit the data better than a density-dependent alternative4. All delays (average healthy lifetime, transfer time between age classes and epidemiological compartments) were assumed to be exponentially distributed to make the stochastic simulation computationally efficient.

All demographic parameters (birth and death rates) were taken from the [Beeton and McCallum, 2011](https://doi.org/10.1111/j.1365-2664.2011.02060.x) and all disease-related parameters were estimated through the fitting process, except for the scaling matrix of the infection rate between age classes.

The local populations were coupled in two different ways: by migration and contact. In a migration event, individuals permanently moved to a neighbouring patch, representing a change in their home-range. With contact, individuals living on neighbouring populations could infect each other with a rate lower than that within the cell, but stayed at their current population. This corresponded to overlapping home ranges that can cross the (arbitrary) cell borders.

All transition rates are shown below.

![alt text](https://github.com/siskavera/tasmanian-devil/tree/master/images/rates.png "Transition rates")

The model was simulated using the continuous-time [Gillespie algorithm](https://doi.org/10.1021/j100540a008), implemented in Matlab R2015b  was used both throughout the fitting procedure and to record the long-term behaviour of the best 200 parameter sets. The algorithm consisted of an initialisation of the state of the system and the rates of each event, and then performing Monte Carlo steps until the predefined ending time of the total simulation was reached. In each step, first the next event was chosen (e.g. a devil on a given cell is infected), with a probability proportional to its rates. Then, a time interval until the event was drawn from an exponential distribution with its parameter equal to the sum of all rates. Finally, both state and time was updated before performing the next Monte Carlo step. A flowchart of the algorithm is shown below.

![alt text](https://github.com/siskavera/tasmanian-devil/tree/master/images/gillespie.png "Gillespie algorithm")

## Field data for demographic parameters
To set the local demography, we used the birth and death rates in an age-structured model with five yearly age classes, as estimated from current knowledge of devil life history and some basic modelling by [Beeton and McCallum, 2011](https://doi.org/10.1111/j.1365-2664.2011.02060.x)

|                    |0-1 year-old   |	1-2 year-old	| 2-3 year-old | 3-4 year-old | 4< year-old |
|--------------------|---------------|---------------|---------------|--------------|-------------|
|Birth rate [1/year] |	0            |	0.13         |	1.3	         | 1.65	         | 1.21       |
|Death rate [1/year] | 0.1925703	   | 0.1904623     | 0.1610143     | 0.1491002	   | 0.3277331  |

We also used published data to set some spatial characteristics of the model. We used the probability of occurrence of Tasmanian devils in 1 km2 blocks, presented by [Hawkins et al., 2006](https://doi.org/10.1016/j.biocon.2006.04.010), as a proxy for the pre-disease population density. This density map was predicted by species environmental domain analysis on the basis of results from a combination of presence survey techniques, and adapted from [Jones & Rose](https://trove.nla.gov.au/work/22862526?selectedversion=NBD13173214). We then scaled these densities such that the overall pre-disease population size of devils in Tasmania was 60,000 (Nicolas Beeton, PhD thesis).

Last, to obtain estimates of the ratio of disease transmission rates between different age classes, we used published scaling factors from [Beeton and McCallum, 2011](https://doi.org/10.1111/j.1365-2664.2011.02060.x). Animals below 1 year-old were not susceptible or infectious, the rate for those in the 1-2 year-old age class were scaled by a factor of 0.602, while the rest of the age classes had the same rate.

