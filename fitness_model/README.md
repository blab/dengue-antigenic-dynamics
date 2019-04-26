# About the fitness model

In meaningfully antigenically diverse viral populations, antigenic novelty (relative to standing population immunity) contributes to viral fitness: as a given virus _i_ circulates in a population, the proportion of the population that is susceptible to infection with _i_--and other viruses antigenically similar to _i_--decreases over time as more people acquire immunity.
Antigenically novel viruses that are able to escape this population immunity are better able to infect hosts and sustain transmission chains, making them fitter than the previously circulating viruses.
Thus, if antigenic novelty constitutes a fitness advantage for DENV, then we would expect greater antigenic distance from recently circulating viruses to correlate with higher growth rates.
We use a simple model, adapted from Luksza and Lassig (_Nature_, 2014) and illustrated below, to estimate population immunity and viral fitness. We then use viral fitness to predict clade frequencies and growth rates.

# Running the fitness model

### Prepare input files

_1 - Estimate empirical clade frequencies over time by [running augur](../augur/)._  
_2 - [Parse clade frequencies](../data_wrangling_scripts/parse-augur-frequencies-output.ipynb) to identify which clades correspond to serotypes and genotypes of interest._  
_3 - Parse the titer model output (also output from augur) to estimate the antigenic distance between each pair of clades._  
**N.B.: For this dataset, steps 1-3 have already been completed; prepared input files for the fitness model can be found [here](../data/frequencies/)**

### Fit model parameters

_4 - Explore parameter space_  
We fit parameters to minimize the root mean squared error between predicted and actual clade frequencies using the Nelder-Mead algorithm as implemented in SciPy. To run parameter fitting:
`python2 antigenic_fitness.py --mode fit`

**beta** Slope of linear relationship between population immunity and viral fitness  
**gamma** Slope of linear relationship between titers and probability of protection  
**sigma** Proportion of titers waning each year since primary infection  
**f\_{s0}** Relative intrinsic fitness of each serotype (f_0 = 0 for DENV4)  
**N** Number of years of previous immunity that contribute to antigenic fitness  
**dt** Number of years in the future to predict clade frequencies

### Run the fitness model

6 - Run `antigenic_fitness.py -h` to see all available options. Simply running with the default options will run the serotype model with parameters fit to optimize `RMSE`.

7 - Explore model output and visualize results with [this notebook](../figures/fitness-frequencies.ipynb)
