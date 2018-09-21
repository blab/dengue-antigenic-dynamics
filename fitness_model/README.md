# About the fitness model  
In meaningfully antigenically diverse viral populations, antigenic novelty (relative to standing population immunity) contributes to viral fitness: as a given virus _i_ circulates in a population, the proportion of the population that is susceptible to infection with _i_--and other viruses antigenically similar to _i_--decreases over time as more people acquire immunity.
Antigenically novel viruses that are able to escape this population immunity are better able to infect hosts and sustain transmission chains, making them fitter than the previously circulating viruses.
Thus, if antigenic novelty constitutes a fitness advantage for DENV, then we would expect greater antigenic distance from recently circulating viruses to correlate with higher growth rates.
We use a simple model, adapted from Luksza and Lassig (_Nature_, 2014) and illustrated below, to estimate population immunity and viral fitness. We then use viral fitness to predict clade frequencies and growth rates.

![Fitness model overview](https://raw.githubusercontent.com/blab/dengue-antigenic-dynamics/master/figures/png/serotype_fitness_model.png)
**A** The relative frequency of each serotype in Southeast Asia is estimated every three months based on available sequence data.  
**B** We calculate antigenic fitness for each serotype over time as its frequency-weighted antigenic distance from recently circulating viruses.  
**C** illustrates how the model predicts clade growth rates. At each timepoint _t_, we blind the model to all empirical data from timepoints later than _t_ and predict each serotype's future trajectory based on its initial frequency, time-invariant intrinsic fitness, and antigenic fitness at time _t_. We predict forward in three-month increments for a total prediction period of _dt = 5_ years.
At each increment, we use the predicted stepwise frequency change to adjust our estimates of antigenic fitness on a rolling basis.
Predicted growth rates are calculated as the predicted final frequency over the actual initial frequency. These predicted growth rates are  compared to empirically observed growth rates to assess model performance. The example illustrated in **C** is also shown in **D** as the blue point.

# Running the fitness model  
### Prepare input files
_1 - Estimate empirical clade frequencies over time by [running augur](../augur/)._   
_2 - [Parse clade frequencies](./helper_scripts/parse-frequencies.ipynb) to identify which clades correspond to serotypes and genotypes of interest._  
_3 - Parse the titer tree (also output from augur) to estimate the antigenic distance between each pair of clades (sum values of `dTiter` for each branch that lies between the two clades on the tree)._  
**N.B.: For this dataset, steps 1-3 have already been completed; prepared input files for the fitness model can be found [here](../data/frequencies/)**  

### Fit model parameters
_4 - Explore parameter space_  
_(Assuming you are on a standard scientific computing cluster using the slurm queue manager)_
_`helpers_scripts$ python sbatch_param_wrapper.py clade_resolution antigenic_resolution`_   
_`clade_resolution` should be either `serotype` or `genotype`._  
_`antigenic_resolution` should be either `interserotype_model` or `fulltree_model`_  

_5 - Examine model performance_  
_Load the out file, `model_performance.csv` into the visualization notebook [here](../helper_scripts/profile-likelihoods.ipynb) to visualize model performance as defined by a variety of optimization metrics._  
  
**N.B.: For this dataset, steps 4-5 have already been completed; model performance for the entire parameter space can be found [here](./profile_likelihoods/). Optimized parameter values are shown below:**
  
**beta** Slope of linear relationship between population immunity and viral fitness   
**gamma** Slope of linear relationship between titers and probability of protection   
**sigma** Proportion of titers waning each year since primary infection   
**f_{s0}** Relative intrinsic fitness of each serotype (f_0 = 0 for DENV4)   
**N** Number of years of previous immunity that contribute to antigenic fitness   
**dt** Number of years in the future to predict clade frequencies   
  
|Genetic resolution | Antigenic resolution | Metric | Metric value | beta | gamma | sigma | DENV1 f_0 | DENV2 f_0 | DENV3 f_0 | DENV4 f_0| 
|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|-------|
|Serotype | Interserotype |    Delta  SSE | 15.02 | 2.57 | 0.57 | 0.86 | 4.57 | 3.43 | 2.14 | 0.00|   
|Serotype | Interserotype | Pearson   R^2  | 0.63 | 2.57 | 0.57 | 0.86 | 3.43 | 2.29 | 0.71 | 0.00|   
|Genotype | Interserotype |    Delta  SSE | 14.83 | 2.57 | 0.57 | 0.86 | 5.71 | 4.57 | 3.57 | 0.00|   
|Genotype | Interserotype | Pearson   R^2  | 0.36 | 2.57 | 0.57 | 0.86 | 5.71 | 5.71 | 2.86 | 0.00|  
|Genotype | Full tree |    Delta  SSE | 14.22 | 1.71 | 0.57 | 0.43 | 1.40 | 0.80 | 0.40 | 0.00|  
|Genotype | Full tree | Pearson   R^2  | 0.33 | 1.29 | 0.57 | 0.43 | 1.40 | 1.60 | 0.40 | 0.00|    
  
### Run the fitness model  
6 - Run `antigenic_fitness.py -h` to see all available options. Simply running with the default options will run the serotype model with parameters fit to optimize `delta_sse`.  
 
7 - Explore model output and visualize results with [this notebook](../figures/fitness-frequencies.ipynb)
