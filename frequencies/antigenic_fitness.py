import matplotlib as mpl
mpl.use('Agg')
import argparse
import pandas as pd
import numpy as np
from scipy import stats
from random import choice
from collections import defaultdict
import matplotlib.pyplot as plt
import seaborn as sns
from math import ceil
from itertools import product
from copy import deepcopy
from pprint import pprint

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

def clade_population_immunity(cls, i):
    ''' for clade i, estimate the relative population immunity at each timepoint based on
    which clades (j) have circulated previously;
    how antigenically distant i is from j;
    and the relative frequency of j'''

    def sum_over_j(i, past_timepoint):
        ''' Return a frequency-weighted sum of the probability of protection from i given prior exposure to j '''
        antigenic_distance = [ cls.titers[tuple(sorted([i,j]))] if i != j else 0. for j in cls.clades] # Pull precomputed antigenic distance between i and j
        probability_protected = [ max(-1.*cls.sigma*Dij + 1., 0.) for Dij in antigenic_distance ] # Linear transformation from antigenic distance to probability of protection from i given prior exposure to j
        j_frequencies = [ cls.frequencies[j][past_timepoint] for j in cls.clades] # Pull relative frequency of each clade j at the timepoint of interest
        return sum( [ j_frequency * prob_protected for (j_frequency, prob_protected) in zip(j_frequencies, probability_protected)]) # return weighted sum

    def sum_over_past_t(i, current_timepoint):
        ''' For each timepoint, look at the past `tp_back` number of timepoints and add up the relative immunity acquired in each interval.
        Adjust for how long ago the population was exposed by assuming that immunity wanes linearly with slope gamma per year (n)'''

        tp_idx = cls.timepoints.index(current_timepoint) # index of timepoint of interest
        past_timepoints = cls.timepoints[tp_idx - cls.tp_back : tp_idx] # previous timepoints to sum immunity over
        exposure = [ sum_over_j(i, t) for t in past_timepoints ] # total protection acquired at each past timepoint, t: sum over all clades for each past timepoint
        waning = [max(-1.*cls.gamma*(current_timepoint - t) + 1., 0.) for t in past_timepoints] # proportion of protection originally acquired at time t expected to remain by the current_timepoint
        return sum( [ w*e for (w,e) in zip(waning, exposure)] ) # sum up the total waning-adjusted population immunity as of the timepoint_of_interest

    valid_timepoints = cls.timepoints[cls.tp_back:]
    exposure = { t: sum_over_past_t(i, t) for t in valid_timepoints }
    return exposure

def predict_timepoint(initial_frequency, initial_fitness, years_forward):
    # if initial_frequency < 0.1:
    #     return np.nan
    # else:
    return initial_frequency*np.exp(initial_fitness*years_forward)

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
                                                         cls.years_forward)
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

def predict_clade_trajectory(cls, i, initial_timepoint):
    '''
    Predict the frequency of clade i at each time interval between t and t+years_forward,
    based on its fitness and frequency at time t
    '''

    dt_values = [ (1./cls.tppy)*dt for dt in range(1, cls.tp_forward+1)] # fraction of year per timepoint * number of timepoints forward
    predicted_timepoints = [ initial_timepoint + dt for dt in dt_values ]
    initial_frequency = cls.frequencies[i][initial_timepoint]
    initial_fitness = cls.fitness[i][initial_timepoint]

    predicted_trajectory = [ predict_timepoint(initial_frequency, initial_fitness, dt) for dt in dt_values ]

    return pd.Series(predicted_trajectory, index=predicted_timepoints, name=i)

def calc_information_gain(cls):
    ''' How much better were our predictions than the null model for time t+N? '''

    def kl_divergence(cls, valid_clades, starting_timepoint, null):
        ''' Kullback-Leibler divergence '''
        kl_div = 0.
        for clade in valid_clades:
            actual_frequency = cls.frequencies[clade][starting_timepoint + cls.years_forward]
            if null == True:
                predicted_frequency = cls.frequencies[clade][starting_timepoint]
            else:
                predicted_frequency = cls.predicted_rolling_frequencies[clade][starting_timepoint + cls.years_forward]
            kl_div += actual_frequency * np.log(actual_frequency / predicted_frequency)
        return kl_div

    information_gain = 0.
    for starting_timepoint in cls.timepoints[cls.tp_back: -1*cls.tp_forward]:
        valid_clades = [c for c in cls.clades if cls.frequencies[c][starting_timepoint] >= 0.1 ]
        m = float(len(valid_clades))
        model_kl_div = kl_divergence(cls, valid_clades, starting_timepoint, null=False)
        null_kl_div = kl_divergence(cls, valid_clades, starting_timepoint, null=True)
        information_gain +=  m*(-1.*model_kl_div + null_kl_div)
    return information_gain

def calc_accuracy(cls):

    def clade_accuracy(i):
        mask = (~np.isnan(cls.actual_growth_rates[i]) & (~np.isnan(cls.predicted_growth_rates[i])))
        actual, predicted = cls.actual_growth_rates[i][mask], cls.predicted_growth_rates[i][mask]

        correct_predictions = 0.
        for (a,p) in zip(actual, predicted):
            if (a >= 1. and p>=1.) or (a < 1. and p < 1.):
                correct_predictions += 1.
        return {'n_correct': correct_predictions, 'n_total': float(len(actual))}

    n_correct = 0.
    n_total = 0.
    for i in cls.clades:
        accuracy = clade_accuracy(i)
        n_correct += accuracy['n_correct']
        n_total += accuracy['n_total']
    return n_correct / n_total

def calc_model_performance(cls, metric=None):

    def remove_nan(actual, predicted):
        actual = actual.loc[actual.index.isin(predicted.index.values)]
        assert predicted.columns.tolist() == actual.columns.tolist()
        assert actual.index.tolist() == predicted.index.tolist()
        actual, predicted = actual.values.flatten(), predicted.values.flatten()
        mask = (~np.isnan(actual)) & (~np.isnan(predicted))
        return actual[mask], predicted[mask]

    actual_freq, predicted_freq = remove_nan(cls.frequencies, cls.predicted_rolling_frequencies)
    actual_growth, predicted_growth = remove_nan(cls.actual_growth_rates, cls.predicted_growth_rates)

    performance = {
    'pearson_r2': stats.linregress(actual_growth, predicted_growth)[2]**2,
    'spearman_r': stats.spearmanr(actual_growth, predicted_growth)[0],
    'abs_error': sum([abs(a - p) for (a,p) in zip(actual_freq, predicted_freq)]) / float(len(actual_freq)),
    'information_gain': calc_information_gain(cls),
    'accuracy': calc_accuracy(cls)}

    return performance

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
        exposure = pd.DataFrame({i: clade_population_immunity(self, i) for i in self.clades})
        self.fitness = -1.*self.beta*exposure
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

        all_trajectories = pd.DataFrame({ i : predict_clade_trajectory(self, i, initial_timepoint)
                           for i in self.clades})
        self.trajectories[initial_timepoint] = normalize_frequencies_by_timepoint(all_trajectories)

def plot_fitness_v_frequency(cls):
    sns.set_palette('tab20', n_colors=20)

    fig, axes = plt.subplots(2,1,figsize=(8,6), sharex=True)

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

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12,8), sharex=True)
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
    ax.text(0,0.2,'r^2 = %.2f'%fit[2]**2, transform=ax.transAxes)
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
        fig, ax = plt.subplots(1,1,figsize=(12,4), sharey=True)
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
    plt.clf()
    plt.close()

def plot_profile_likelihoods(model_performance, metric, args):
    fit_params = ['beta', 'gamma', 'sigma']

    if metric != 'abs_error':
        best_fit = model_performance.ix[model_performance[metric].idxmax()]
    else:
        best_fit = model_performance.ix[model_performance[metric].idxmin()]

    print 'Best fit: ', best_fit
    fig, axes = plt.subplots(ncols=len(fit_params), nrows=1, figsize=(3*len(fit_params), 3))
    for param,ax in zip(fit_params, axes):
        p1,p2 = [p for p in fit_params if p != param]
        plot_vals = model_performance.loc[(model_performance[p1]==best_fit[p1]) & (model_performance[p2]==best_fit[p2])]

        sns.regplot(param, metric, data=plot_vals, fit_reg=False, ax=ax)
        ax.set_title('Fixed params:\n%s = %.1f,\n%s=%.1f'%(p1, best_fit[p1], p2, best_fit[p2]))
        ax.set_xlabel(param)
        ax.set_ylabel('%s'%metric)
        plt.tight_layout()

    if args.save:
        plt.savefig(args.out_path+args.name+'_profile_likelihoods_%s.png'%metric, bbox_inches='tight', dpi=300)
    else:
        plt.show()

    plt.clf()
    plt.close()

def plot_model_performance(model_performance, metric, args, small_multiples_var='sigma', ):
    small_multiples_vals = pd.unique(model_performance[small_multiples_var])
    nplots = len(small_multiples_vals)
    nrows = max(int(ceil(nplots/5)), 1)

    fig, axes = plt.subplots(ncols=min(5, nplots), nrows=nrows, figsize=(3*5, 3*nrows))

    fit_params = ['beta', 'gamma', 'sigma']
    x_var, y_var = sorted([v for v in fit_params if v not in [metric, small_multiples_var]])
    vmin, vmax = model_performance[metric].min(), model_performance[metric].max()

    for value, ax in zip(small_multiples_vals, axes.flatten()):
        plot_values = model_performance.loc[model_performance[small_multiples_var] == value]
        plot_values = plot_values.pivot(index=x_var, columns=y_var, values=metric)
        sns.heatmap(plot_values, vmin=vmin, vmax=vmax, ax=ax,cbar_kws={'label': '%s'%metric})
        ax.set_title('%s = %f'%(small_multiples_var, value))
        ax.set_ylabel(x_var)
        ax.set_xlabel(y_var)

    plt.tight_layout()

    if args.save:
        plt.savefig(args.out_path+args.name+'_param_performance_%s.png'%metric, bbox_inches='tight', dpi=300)
    else:
        plt.show()

    plt.clf()
    plt.close()

def run_model(args):
    antigenic_fitness = AntigenicFitness(args)
    if not isinstance(antigenic_fitness.fitness, pd.DataFrame):
        print 'calculating fitness'
        antigenic_fitness.calculate_fitness()
    print 'predicting frequencies'
    antigenic_fitness.predict_rolling_frequencies()
    print 'calculating growth rates'
    antigenic_fitness.calc_growth_rates()

    if antigenic_fitness.plot == True:
        print 'generating plots'
        plot_fitness_v_frequency(antigenic_fitness)
        plot_rolling_frequencies(antigenic_fitness)
        plot_growth_rates(antigenic_fitness)
        plot_trajectory_multiples(antigenic_fitness)

    return calc_model_performance(antigenic_fitness, args.metric)

def test_parameter_grid(args):

    def get_range(parameter):
        p = vars(args)[parameter]
        if type(p) == list and len(p) == 2:
            return np.linspace(p[0], p[1], args.n_param_vals)
        elif type(p) == list:
            return p[0]
        else:
            return p

    fit_params = ['beta', 'gamma', 'sigma']
    def run_with_parameterization(parameterization):
        args_copy = deepcopy(args)
        setattr(args_copy, 'save', False)
        setattr(args_copy, 'plot', False)
        for attr, val in zip(fit_params, parameterization):
            setattr(args_copy, attr, val)
        print 'Running model with parameters:\n', zip(fit_params, parameterization), '\n'

        results = run_model(args_copy)
        results.update({k: v for k,v in zip(fit_params, parameterization)})
        print results, '\n\n'
        return results

    beta_vals, gamma_vals, sigma_vals = get_range('beta'), get_range('gamma'), get_range('sigma')
    model_parameterizations = product(beta_vals, gamma_vals, sigma_vals) # [(b0, g0, s0), (b0, g0, s1), ....]

    model_performance = pd.DataFrame([ run_with_parameterization(parameterization)
                          for parameterization in model_parameterizations]).round(3)

    if args.save:
        model_performance.to_csv(args.out_path+args.name+'_model_performance.csv')

    if args.plot:
        metrics_to_plot = [args.metric] if args.metric else model_performance.columns.values
        print metrics_to_plot
        for metric in metrics_to_plot:
            if metric not in fit_params:
                plot_profile_likelihoods(model_performance, metric, args)
                plot_model_performance(model_performance, metric, args)

    return model_performance

if __name__=="__main__":
    sns.set(style='whitegrid')#, font_scale=1.5)

    args = argparse.ArgumentParser()
    args.add_argument('--frequency_path', help='frequencies csv', default='southeast_asia/serotype/southeast_asia_serotype_frequencies.csv')
    args.add_argument('--titer_path', help='pairwise dTiters csv', default='all_effects_Dij.csv')
    args.add_argument('--fitness_path', type=str, help='path to precomputed frequencies or \'null\'')
    args.add_argument('--date_range', nargs=2, type=float, help='which dates to look at', default=[1970., 2015.])
    args.add_argument('--years_back', type=int, help='how many years of past immunity to include in fitness estimates', default=3)
    args.add_argument('--years_forward', type=int, help='how many years into the future to predict', default=3)
    args.add_argument('--gamma', nargs='*', type=float, help='Value or value range for -1*proportion of titers that wane per year post-exposure (slope of years vs. p(titers remaining))', default= 0.15)
    args.add_argument('--sigma', nargs='*', type=float, help='Value or value range for -1*probability of protection from i conferred by each log2 titer unit against i', default= 1.2)
    args.add_argument('--beta', nargs='*', type=float, help='Value or value range for beta. fitness = -1.*beta*population_immunity', default=1.)
    args.add_argument('--n_param_vals', type=int, help='Number of values to test for each parameter if fitting model', default=3)
    args.add_argument('--metric', help='Metric to use when fitting parameters', choices=['pearson_r2', 'spearman_r', 'abs_error', 'information_gain', 'accuracy'], default=None)
    args.add_argument('--plot', help='make plots?', action='store_true')
    args.add_argument('--save', help='save csv and png files?', action='store_true')
    args.add_argument('--name', type=str, help='analysis name')
    args.add_argument('--out_path', type=str, help='where to save csv and png files', default='./')
    args = args.parse_args()


    ## If given a range of parameters, test all combinations.
    if any([ type(args.beta) == list and len(args.beta) == 2,
             type(args.gamma) == list and len(args.gamma) == 2,
             type(args.sigma) == list and len(args.sigma) == 2  ]):

        model_performance = test_parameter_grid(args)

        if args.metric: # If evaluation metric specified, use this to select the best fit and do a clean run of the model
            if args.metric != 'abs_error':
                best_fit = model_performance.ix[model_performance[args.metric].idxmax()]
            else:
                best_fit = model_performance.ix[model_performance[args.metric].idxmin()]
            setattr(args, 'beta', best_fit['beta'])
            setattr(args, 'gamma', best_fit['gamma'])
            setattr(args, 'sigma', best_fit['sigma'])
            run_model(args)

    else: # If given specific values for paramters, just run the model and print performance metrics.
        model_performance = run_model(args)
        print(model_performance)
