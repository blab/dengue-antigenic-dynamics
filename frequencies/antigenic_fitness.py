import argparse
import pandas as pd
import numpy as np
import baltic as bt
from scipy import stats
from random import choice
from collections import defaultdict
from pprint import pprint
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.optimize import minimize
from math import ceil
from itertools import product
from copy import deepcopy

def normalize_frequencies_by_timepoint(frequencies):
    ''' Normalize each row so that the sum of all frequencies at a single timepoint = 1'''
    def normalize(row):
        total = row.sum()
        if np.isnan(total) or total == 0:
            return row
        else:
            return row.map( lambda x: x / total)

    if isinstance(frequencies, dict):
        frequencies = pd.DataFrame(frequencies)
        normalized_frequencies = frequencies.apply(normalize, axis=1)
        return normalized_frequencies.to_dict()

    else:
        normalized_frequencies = frequencies.apply(normalize, axis=1)
        return normalized_frequencies

def sum_over_j(cls, i, timepoint, proportion_remaining):
    '''
    Look at all a single time point.
    At that timepoint, look at all clades, j, that are cocirculating with i.
    Pull the titers between i and j (D_ij) and adjust for waning immunity (proportion_remaining).
    Use this to calculate a frequency-weighted estimate sum of the probability of protection against i given prior exposure to j.
    '''

    frequency_weighted_protection = 0.

    for j in cls.clades:
        if i==j:
            init_titers = 0.
        else:
            init_titers = cls.titers[tuple(sorted([i,j]))]
        remaining_titers = proportion_remaining * init_titers
        probability_protected = max(cls.sigma*remaining_titers + 1., 0.)
        j_frequency = cls.frequencies[j][timepoint]

        frequency_weighted_protection += probability_protected*j_frequency

    return frequency_weighted_protection

def sum_over_past_t(cls, i, timepoint_of_interest):
    '''
    For a given time point of interest, look at what the population has acquired protection to
    over the past `tp_back` timepoints.
    For each previous timepoint, sum_over_j to calculate the accumulated immunity.
    '''

    def waning(gamma, n):
        ''' Assume immunity wanes linearly with slope gamma per year (n)'''
        return max(gamma*n + 1., 0.)
        # return max( np.exp(gamma*n), 0. )

    tp_idx = cls.timepoints.index(timepoint_of_interest) # index of timepoint of interest
    t_to_sum = cls.timepoints[tp_idx - cls.tp_back : tp_idx] # previous timepoints to sum immunity over

    # proportion of titers acquired in each interval expected to remain by timepoint of interest
    waning_over_time = [ waning(cls.gamma, t - timepoint_of_interest) for t in t_to_sum]

    # proportion of the population that acquired protection in each interval
    protection_over_time = [ sum_over_j(cls, i, t, p_remaining)
                          for (t, p_remaining) in zip(t_to_sum, waning_over_time) ]

    accumulated_protection = sum(protection_over_time)/float(len(protection_over_time))
    return accumulated_protection

def population_exposure(cls, i):
    ''' estimate the proportion of the population that is immune to i at the beginning of season t
    '''

    valid_timepoints = cls.timepoints[cls.tp_back:]

    # if antigenic_resolution == 'null':
    #     return { t: 0. for t in valid_timepoints }
    # else:
    population_exposure = { t: sum_over_past_t(cls, i, t) for t in valid_timepoints}

    return population_exposure

def predict_timepoint(initial_frequency, initial_fitness, years_forward, beta):
    # if initial_frequency < 0.1:
    #     return np.nan
    # else:
    return initial_frequency*np.exp(beta*initial_fitness*years_forward)

def clade_rolling_prediction(cls, i):
    '''
    For each timepoint t, predict the frequency of i based on
    its fitness and initial frequency at time t-years_forward
    '''

    initial_fitnesses = cls.fitness[i]
    initial_frequencies = cls.frequencies[i]
    initial_timepoints = cls.timepoints[cls.tp_back: -1*cls.tp_forward]
    predicted_timepoints = cls.timepoints[cls.tp_forward+cls.tp_back:]

    predicted_frequencies = { pred_t : predict_timepoint(initial_frequencies[init_t],
                                                         initial_fitnesses[init_t],
                                                         cls.years_forward, cls.beta)
                            for (init_t, pred_t)
                            in zip(initial_timepoints, predicted_timepoints) }

    return predicted_frequencies

def clade_rolling_growth_rate(cls,i,predicted=True):

    initial_timepoints = cls.timepoints[cls.tp_back: -1*cls.tp_forward]
    initial_frequencies = cls.frequencies[i][initial_timepoints]

    final_timepoints = cls.timepoints[cls.tp_forward+cls.tp_back:]
    if predicted==True:
        final_frequencies = cls.predicted_rolling_frequencies[i][final_timepoints]
    else:
        final_frequencies = cls.frequencies[i][final_timepoints]

    time_intervals = [ str(f)+'/'+str(i) for (i, f) in zip(initial_timepoints, final_timepoints)]
    growth_rates = [ f / i for (i,f) in zip(initial_frequencies, final_frequencies)]

    return pd.Series(growth_rates, index=time_intervals)

def predict_trajectory(cls, i, initial_timepoint):
    '''
    Predict the frequency of clade i at each time interval between t and t+years_forward,
    based on its fitness and frequency at time t
    '''

    dt_values = [ (1./cls.tppy)*dt for dt in range(1, cls.tp_forward+1)] # fraction of year per timepoint * number of timepoints forward
    predicted_timepoints = [ initial_timepoint + dt for dt in dt_values ]
    initial_frequency = cls.frequencies[i][initial_timepoint]
    initial_fitness = cls.fitness[i][initial_timepoint]

    predicted_trajectory = [ predict_timepoint(initial_frequency, initial_fitness, dt, cls.beta) for dt in dt_values ]

    return pd.Series(predicted_trajectory, index=predicted_timepoints, name=i)

class AntigenicFitness():
    def __init__(self, args):

        for k,v in vars(args).items():
            setattr(self, k, v) # copy over cmd line args

        # actual (observed) frequencies
        self.frequencies = pd.read_csv(self.frequency_path, index_col=0) # pd.DataFrame(index=timepoints, columns=clades, values=relative frequencies)
            self.clades = self.frequencies.columns.values
        self.frequencies = normalize_frequencies_by_timepoint(self.frequencies[self.clades]) # restrict to clades of interest, normalize
        self.frequencies = self.frequencies.loc[(self.frequencies.index >= self.date_range[0]) & (self.frequencies.index <= self.date_range[1])] # restrict to timepoints of interest
        self.timepoints = self.frequencies.index.tolist()
        n_years = int(self.timepoints[-1]) - int(self.timepoints[0]) # number of years in the frequencies dataset
        self.tppy = int(len(self.timepoints)/n_years) # timepoints per year
        self.tp_forward = self.tppy*self.years_forward # number of timepoints forward
        self.tp_back = self.tppy*self.years_back # number of timepoints back

        self.noisy_predictions_mask = self.frequencies < 0.1 # log which initial values were low
        # keep track of which predictions will be made based on low initial values
        self.noisy_predictions_mask.index = self.noisy_predictions_mask.index.map(lambda x: x+self.years_forward)

        if self.fitness_path:
            if self.fitness_path == 'null': # negative control: fitnesses all = 0.
                self.fitness = pd.DataFrame(index=self.timepoints, columns=self.clades)
                self.fitness.fillna(0., inplace=True)
            else: # load from file if provided
                self.fitness = pd.read_csv(self.fitness_path, index_col=0)
        else:
            self.fitness = None

        self.trajectories = {}

        # load pre-computed antigenic distances between clades
        self.titers = {(str(k1), str(k2)):v for (k1,k2),v in pd.Series.from_csv(args.titer_path, header=None,index_col=[0,1]).to_dict().items()}

        if self.save == True:
            assert self.name, 'ERROR: Please provide an analysis name if you wish to save output'


    def calculate_fitness(self):
        ''' fitness = 1.-population exposure'''
        self.fitness = 1. - pd.DataFrame({i: population_exposure(self, i) for i in self.clades})
        if self.save:
            self.fitness.to_csv(self.out_path+self.name+'_fitness.csv')

    def predict_rolling_frequencies(self):
        '''
        Making a "rolling prediction" of frequency for each clade.
        Normalize these predicted frequencies so that they sum to 1. at each timepoint,
        then mask out predictions based on noisy initial frequencies.

        Returns pd.DataFrame(index = t+dt, columns=clades,
                            values = predicted frequency at time t+dt, based on frequency & fitness at time t )
        Xi(t+dt) = Xi(t) * e^( Fi(t) * dt), where Xi is frequency and Fi is fitness of i
        '''
        all_predicted_frequencies = pd.DataFrame({ i : clade_rolling_prediction(self, i)
                                    for i in self.clades })

        self.predicted_rolling_frequencies = normalize_frequencies_by_timepoint(all_predicted_frequencies) # normalize based on ALL tips
        self.predicted_rolling_frequencies = self.predicted_rolling_frequencies[~self.noisy_predictions_mask] # keep only predictions based on initial frequencies at >0.1

        if self.save:
            self.predicted_rolling_frequencies.to_csv(self.out_path+self.name+'_predicted_freqs.csv')

    def calc_growth_rates(self):
        self.predicted_growth_rates = pd.DataFrame({i:clade_rolling_growth_rate(self, i, predicted=True) for i in self.clades})
        self.actual_growth_rates = pd.DataFrame({i:clade_rolling_growth_rate(self, i, predicted=False) for i in self.clades})
        if self.save:
            self.predicted_growth_rates.to_csv(self.out_path+self.name+'_predicted_growth_rates.csv')
            self.actual_growth_rates.to_csv(self.out_path+self.name+'_actual_growth_rates.csv')

    def predict_trajectories(self,initial_timepoint):
        '''
        Predict the frequency of all clades at each time interval between t and t+years_forward,
        based on their initial fitnesses and frequencies at time t.

        Normalize these predicted frequencies so that they sum to 1. at each timepoint.
        '''

        all_trajectories = pd.DataFrame({ i : predict_trajectory(self, i, initial_timepoint)
                           for i in self.clades})
        self.trajectories[initial_timepoint] = normalize_frequencies_by_timepoint(all_trajectories)

def run_model(args):
    antigenic_fitness = AntigenicFitness(args)
    if not isinstance(antigenic_fitness.fitness, pd.DataFrame):
        print 'calculating fitness'
        antigenic_fitness.calculate_fitness()
    print 'predicting frequencies'
    antigenic_fitness.predict_rolling_frequencies()
    print 'calculating growth rates'
    antigenic_fitness.calc_growth_rates()

    actual, predicted = antigenic_fitness.actual_growth_rates, antigenic_fitness.predicted_growth_rates
    actual = actual.loc[actual.index.isin(predicted.index.values)]

    assert predicted.columns.tolist() == actual.columns.tolist()
    assert actual.index.tolist() == predicted.index.tolist()

    actual, predicted = actual.values.flatten(), predicted.values.flatten()
    mask = (~np.isnan(actual)) & (~np.isnan(predicted))
    fit = stats.linregress(actual[mask], predicted[mask])

    if antigenic_fitness.plot == True:
        print 'generating plots'
        plot_fitness_v_frequency(antigenic_fitness)
        plot_rolling_frequencies(antigenic_fitness)
        plot_growth_rates(antigenic_fitness)
        plot_trajectory_multiples(antigenic_fitness)
    print '\n\n'
    return fit[2]**2 #r^2 value

def test_parameter_grid(params, args):
    parameters = sorted(params.keys()) # [ 'beta', 'gamma', 'sigma' ]
    value_ranges = [params[p] for p in parameters] #[ [betamin, ... betamax], [gammamin, ..., gammamax], etc.]
    model_parameterizations = product(*value_ranges) # [(b0, g0, s0), (b0, g0, s1), ....]

    model_performance = {}
    for set_values in model_parameterizations:
        test_args = deepcopy(args)
        for param, value in zip(parameters, set_values):
            setattr(test_args, param, value)
        print 'Running model with parameters:\n', zip(parameters, set_values), '\n'
        model_performance[set_values]= run_model(test_args)

    model_performance = pd.Series(model_performance).reset_index() # tuple(param1_value, param2_value) -> pd.Series(index=(param1_value, param2_value)) -> pd.DataFrame()
    model_performance.columns = sorted(params.keys())+['r^2']

    return model_performance

def plot_fitness_v_frequency(cls):
    sns.set_palette('tab20', n_colors=20)

    fig, axes = plt.subplots(2,1,figsize=(8,6))

    fitness=cls.fitness
    frequencies = cls.frequencies

    for clade in fitness.columns.values:
        axes[0].plot(fitness[clade].index.values, fitness[clade], linestyle='--', label='%s fitness'%clade)
        axes[1].plot(frequencies[clade].index.values, frequencies[clade], linestyle='-', label='%s frequency'%clade)

    axes[0].set_title('Clade fitness')
    axes[1].set_title('Observed clade frequency')
    plt.tight_layout()

    if cls.save:
        plt.savefig(cls.out_path+cls.name+'_fitness.png', dpi=300)
    else:
        plt.show()
    plt.clf()
    plt.close()

def plot_rolling_frequencies(cls):
    sns.set_palette('tab20', n_colors=20)

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12,8))
    for clade, predicted_freqs in cls.predicted_rolling_frequencies.iteritems():
        date_min, date_max = predicted_freqs.index.min(), predicted_freqs.index.max()
        axes[0].plot(predicted_freqs.index.values,predicted_freqs.values, linestyle='--', label='Predicted %s frequencies'%clade)

        actual_frequencies = cls.frequencies[clade][date_min:date_max]
        axes[1].plot(actual_frequencies.index.values, actual_frequencies.values, label='Actual %s frequencies'%clade)

    axes[0].set_title('Predicted clade frequencies')
    axes[1].set_title('Observed clade frequencies')

    plt.tight_layout()
    if cls.save:
        plt.savefig(cls.out_path+cls.name+'_frequencies.png', dpi=300)
    else:
        plt.show()
    plt.clf()
    plt.close()

def plot_growth_rates(cls):
    '''
    For the actual and predicted frequencies, find where both values are non-null and > 0.1
    Plot actual vs. predicted
    '''
    sns.set_palette('Set2', n_colors=10)
    actual, predicted = cls.actual_growth_rates, cls.predicted_growth_rates
    actual = actual.loc[actual.index.isin(predicted.index.values)]

    assert predicted.columns.tolist() == actual.columns.tolist()
    assert actual.index.tolist() == predicted.index.tolist()

    actual, predicted = actual.values.flatten(), predicted.values.flatten()
    mask = (~np.isnan(actual)) & (~np.isnan(predicted))
    fit = stats.linregress(actual[mask], predicted[mask])

    ax=sns.regplot(actual[mask], predicted[mask])
    ax.set_xlabel('Actual growth rate')#, X(t+%d)/X(t)'%years_forward)
    ax.set_ylabel('Predicted growth rate')#, X(t+%d)/X(t)'%years_forward)
    ax.text(0,0.2,'r^2 = %.2f'%fit[2]**2, )
    plt.tight_layout()

    if cls.save:
        plt.savefig(cls.out_path+cls.name+'_growth_rates.png', dpi=300)
    else:
        plt.show()
    plt.clf()
    plt.close()

def plot_trajectory(cls, initial_timepoint, clades, ax):
    cmap = sns.color_palette('Set2', len(clades))
    colors = { clade: cmap[i] for i, clade in enumerate(clades)}

    if ax == None:
        fig, ax = plt.subplots(1,1,figsize=(12,4))
    try:
        predicted_trajectory = cls.trajectories[initial_timepoint]
    except KeyError:
        cls.predict_trajectories(initial_timepoint)
        predicted_trajectory = cls.trajectories[initial_timepoint]

    for clade, predicted_trajectory in predicted_trajectory[clades].iteritems():
        actual_vals = cls.frequencies[clade][:initial_timepoint+cls.tp_forward]
        ax.plot(actual_vals.index.values, actual_vals.values, c=colors[clade], label='Actual frequencies')
        ax.plot(predicted_trajectory.index.values, predicted_trajectory.values, c=colors[clade], linestyle='--', label='Predicted frequencies')
        ax.plot(actual_vals.index.values, cls.fitness[clade][actual_vals.index.values], c=colors[clade], linestyle=':', label='Fitness')
    ax.set_xlim(initial_timepoint-cls.years_back, predicted_trajectory.index.values[-1])
    ax.set_ylim(0,1)

def plot_trajectory_multiples(cls, starting_timepoints=None, n_clades_per_plot=2):
    sns.set(style='whitegrid', font_scale=0.8)

    if starting_timepoints == None:
        starting_timepoints = cls.timepoints[cls.tp_forward+cls.tp_back::cls.tp_forward]

    ncols = len(starting_timepoints)
    nrows = int(ceil(len(cls.clades)/n_clades_per_plot))

    fig, axes = plt.subplots(nrows, ncols, figsize=(3*ncols, 2*nrows),sharex=False, sharey=True)
    clade_sets = [cls.clades[i:i + n_clades_per_plot] for i in xrange(0, len(cls.clades), n_clades_per_plot)]

    for clade_set, row in zip(clade_sets, axes):
        for tp, ax in zip(starting_timepoints, row):
            plot_trajectory(cls, tp, clade_set, ax)
        ax.legend()
    plt.tight_layout()
    if cls.save:
        plt.savefig(cls.out_path+cls.name+'_trajectories.png', bbox_inches='tight', dpi=300)
    else:
        plt.show()

def plot_profile_likelihoods(model_performance, args):
    sns.set_palette('Set2', n_colors=10)

    fit_params = [p for p in model_performance.columns.values if p != 'r^2']

    ml_fit = model_performance.ix[model_performance['r^2'].idxmax()]
    print 'MLE: ', ml_fit
    fig, axes = plt.subplots(ncols=len(fit_params), nrows=1, figsize=(3*len(fit_params), 3))
    for param,ax in zip(fit_params, axes):
        p1,p2 = [p for p in fit_params if p != param]
        plot_vals = model_performance.loc[(model_performance[p1]==ml_fit[p1]) & (model_performance[p2]==ml_fit[p2])]

        sns.regplot(param, 'r^2', data=plot_vals, fit_reg=False, ax=ax)
        ax.set_title('Fixed params:\n%s = %.1f,\n%s=%.1f'%(p1, ml_fit[p1], p2, ml_fit[p2]))
        ax.set_xlabel(param)
        ax.set_ylabel('R^2')
    plt.tight_layout()

    # if args.save:
    plt.savefig(args.out_path+args.name+'_profile_likelihoods.png', bbox_inches='tight', dpi=300)
    # else:
    #     plt.show()

def plot_param_performance(model_performance, args, small_multiples_var='sigma', ):
    small_multiples_vals = pd.unique(model_performance[small_multiples_var])
    nplots = len(small_multiples_vals)
    nrows = max(int(ceil(nplots/5)), 1)

    fig, axes = plt.subplots(ncols=min(5, nplots), nrows=nrows, figsize=(3*5, 3*nrows))

    x_var, y_var = sorted([v for v in model_performance.columns.values if v not in ['r^2', small_multiples_var]])
    vmin, vmax = model_performance['r^2'].min(), model_performance['r^2'].max()

    for value, ax in zip(small_multiples_vals, axes.flatten()):
        plot_values = model_performance.loc[model_performance[small_multiples_var] == value]
        plot_values = plot_values.pivot(index=x_var, columns=y_var, values='r^2')
        sns.heatmap(plot_values, vmin=vmin, vmax=vmax, ax=ax,cbar_kws={'label': 'R^2'})
        ax.set_title('%s = %f'%(small_multiples_var, value))
        ax.set_ylabel(x_var)
        ax.set_xlabel(y_var)

    plt.tight_layout()

    # if args.save:
    plt.savefig(args.out_path+args.name+'_param_performance.png', bbox_inches='tight', dpi=300)
    # else:
    #     plt.show()

if __name__=="__main__":
    sns.set(style='whitegrid')#, font_scale=1.5)

    args = argparse.ArgumentParser()
    args.add_argument('--frequency_path', help='frequencies csv', default='../../data/titer-model/frequencies/southeast_asia_clade_frequencies.csv')
    args.add_argument('--titer_path', help='pairwise dTiters csv', default='../../data/titer-model/frequencies/clade_dtiters.csv')
    args.add_argument('--fitness_path', type=str, help='path to precomputed frequencies or \'null\'')
    args.add_argument('--date_range', nargs=2, type=float, help='which dates to look at', default=[1970., 2015.])
    args.add_argument('--years_back', type=int, help='how many years of past immunity to include in fitness estimates', default=3)
    args.add_argument('--years_forward', type=int, help='how many years into the future to predict', default=3)
    args.add_argument('--gamma', nargs='*', type=float, help='Value or value range for -1*proportion of titers that wane per year post-exposure (slope of years vs. p(titers remaining))', default= -0.15)
    args.add_argument('--sigma', nargs='*', type=float, help='Value or value range for -1*probability of protection from i conferred by each log2 titer unit against i', default= -1.)
    args.add_argument('--beta', nargs='*', type=float, help='Value or value range for beta. fitness = 1. - beta*population_exposure', default= 2.)
    args.add_argument('--n_param_vals', type=int, help='Number of values to test for each parameter if fitting model', default=3)
    args.add_argument('--plot', help='make plots?', action='store_true')
    args.add_argument('--save', help='save csv and png files?', action='store_true')
    args.add_argument('--name', type=str, help='analysis name')
    args.add_argument('--out_path', type=str, help='where to save csv and png files', default='./')
    args = args.parse_args()

    # dataset_mle = {'serotype': {'beta': 11, 'gamma': -1.2, 'sigma': -1.05},
    #                'genotype': {'beta': , 'gamma': , 'sigma': },
    #                'flu':      {'beta': , 'gamma': , 'sigma': }}

    # if args.clades[0] in dataset_clades.keys():
    #     setattr(args, 'clades', dataset_clades[args.clades[0]])
    # assert len(args.clades) >= 2, "ERROR: clades must be either a dataset in ['serotype', 'genotype', 'all'] or a list of clade IDs"

    parameter_grid = {}
    for param in ['beta', 'sigma', 'gamma']:
        param_val = vars(args)[param]
        if len(param_val) == 2:
            parameter_grid[param] = np.linspace(param_val[0], param_val[1], args.n_param_vals)
        else:
            setattr(args, param, param_val[0])

    if parameter_grid != {}:
        model_performance = test_parameter_grid(parameter_grid, args)
        setattr(args, 'plot', False)
        setattr(args, 'save', False)
        plot_param_performance(model_performance, args)
        plot_profile_likelihoods(model_performance, args)

        ml_fit = model_performance.ix[model_performance['r^2'].idxmax()]
        setattr(args, 'beta', ml_fit['beta'])
        setattr(args, 'gamma', ml_fit['gamma'])
        setattr(args, 'sigma', ml_fit['sigma'])
        setattr(args, 'save', True)
        setattr(args, 'plot', True)
        run_model(args)

    else:
        run_model(args)
