# Predicting the fate of the Tasmanian devil with a spatial epidemiological model

Spatially explicit, individual-based framework to model the spread of DFTD in Tasmania, implemented in Matlab.

## Metapopulation model

Our spatially explicit model consisted of an individual-based (thus stochastic) metapopulation of local populations, arranged on a square grid of 10km by 10km cells representing mainland Tasmania. The local populations on each grid were assumed to be well-mixed and only nearest-neighbour interactions (contact and migration) in the von Neumann neighbourhood (4 immediately neighbouring cells) were allowed. The carrying capacity of each cell was set using published data from [Hawkins et al., 2006](https://doi.org/10.1016/j.biocon.2006.04.010). The structure of the model is shown below.
![model structure](https://github.com/siskavera/tasmanian-devil/tree/master/images/model_structure.png)

The local dynamics consisted of a compartmental epidemiological model built on demography following logistic growth. There were four disease-related compartments: susceptible (S), exposed (E, no symptoms or infectivity), diagnosable (D, show symptoms but not infective) and infective. The cancer was assumed to lead to certain death, since very few recoveries have ever been observed4. We used a frequency-dependent transmission rate, since it was shown to fit the data better than a density-dependent alternative4. All delays (average healthy lifetime, transfer time between age classes and epidemiological compartments) were assumed to be exponentially distributed to make the stochastic simulation computationally efficient.

All demographic parameters (birth and death rates) were taken from the [Beeton et al., 2011](https://doi.org/10.1111/j.1365-2664.2011.02060.x) and all disease-related parameters were estimated through the fitting process, except for the scaling matrix of the infection rate between age classes.

The local populations were coupled in two different ways: by migration and contact. In a migration event, individuals permanently moved to a neighbouring patch, representing a change in their home-range. With contact, individuals living on neighbouring populations could infect each other with a rate lower than that within the cell, but stayed at their current population. This corresponded to overlapping home ranges that can cross the (arbitrary) cell borders.

All transition rates are shown below.
![model structure](https://github.com/siskavera/tasmanian-devil/tree/master/images/rates.png)

The model was simulated using the continuous-time [Gillespie algorithm](https://doi.org/10.1021/j100540a008), implemented in Matlab R2015b  was used both throughout the fitting procedure and to record the long-term behaviour of the best 200 parameter sets. The algorithm consisted of an initialisation of the state of the system and the rates of each event, and then performing Monte Carlo steps until the predefined ending time of the total simulation was reached. In each step, first the next event was chosen (e.g. a devil on a given cell is infected), with a probability proportional to its rates. Then, a time interval until the event was drawn from an exponential distribution with its parameter equal to the sum of all rates. Finally, both state and time was updated before performing the next Monte Carlo step. A flowchart of the algorithm is shown below.
![model structure](https://github.com/siskavera/tasmanian-devil/tree/master/images/gillespie.png)

