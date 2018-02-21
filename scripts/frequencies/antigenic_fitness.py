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

class AntigenicFitness():
    def __init__(self, args):

        self.date_range = args.date_range # which dates to look at

        # actual (observed) frequencies
        self.frequencies = pd.read_csv(args.frequency_path, index_col=0) # pd.DataFrame(index=timepoints, columns=clades, values=relative frequencies)

        if not args.clades:
            self.clades = self.frequencies.columns.values
        else:
            self.clades = args.clades # which non-overlapping clades to look at

        self.frequencies = normalize_frequencies_by_timepoint(self.frequencies[self.clades]) # restrict to clades of interest, normalize

        self.frequencies = self.frequencies.loc[(self.frequencies.index >= self.date_range[0]) & (self.frequencies.index <= self.date_range[1])] # restrict to timepoints of interest

        self.years_forward = args.years_forward # how many years forward to try and predict (rolling)
        self.years_back = args.years_back # how many years of exposure to account for in fitness estimations
        self.timepoints = self.frequencies.index.tolist()
        n_years = int(self.timepoints[-1]) - int(self.timepoints[0]) # number of years in the frequencies dataset
        self.tppy = int(len(self.timepoints)/n_years) # timepoints per year
        self.tp_forward = self.tppy*self.years_forward # number of timepoints forward
        self.tp_back = self.tppy*self.years_back # number of timepoints back

        self.noisy_predictions_mask = self.frequencies < 0.1 # log which initial values were low
        # keep track of which predictions will be made based on low initial values
        self.noisy_predictions_mask.index = self.noisy_predictions_mask.index.map(lambda x: x+self.years_forward)

        if args.fitness:
            if args.fitness == 'null': # negative control: fitnesses all = 0.
                self.fitness = pd.DataFrame(index=self.timepoints, columns=self.clades)
                self.fitness.fillna(0., inplace=True)
            else: # load from file if provided
                self.fitness = pd.read_csv(args.fitness, index_col=0)
        else:
            self.fitness = None

        # load pre-computed antigenic distances between clades
        self.titers = {(str(k1), str(k2)):v for (k1,k2),v in pd.Series.from_csv(args.dTiters_path, header=None,index_col=[0,1]).to_dict().items()}

        self.sigma=args.sigma # slope of C(Dij) (protection vs. titers)
        self.gamma=args.gamma # slope of waning(n) (time vs. titers)

        self.plot=args.plot
        self.save=args.save
        self.name=args.name
        self.out_path=args.out_path
        if self.save:
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

def plot_fitness_v_frequency(cls):

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
    ax.text(0,0.2,'r = %.2f'%fit[2], )
    plt.tight_layout()

    if cls.save:
        plt.savefig(cls.out_path+cls.name+'_growth_rates.png', dpi=300)
    else:
        plt.show()
    plt.clf()
    plt.close()

if __name__=="__main__":
    sns.set(style='whitegrid', font_scale=1.5)
    sns.set_palette('tab20', n_colors=20)

    args = argparse.ArgumentParser()
    args.add_argument('-f', '--frequency_path', help='frequencies csv', default='../../data/titer-model/frequencies/southeast_asia_clade_frequencies.csv')
    args.add_argument('-t', '--dTiters_path', help='pairwise dTiters csv', default='../../data/titer-model/frequencies/clade_dtiters.csv')
    args.add_argument('--fitness', type=str, help='path to precomputed frequencies or \'null\'')
    args.add_argument('-c', '--clades', nargs='*', type=str, help='which clades to look at', default=['2185', '2589', '2238', '2596', '1460', '1393', '1587', '1455', '975', '979', '1089', '33', '497', '117', '543', '4', '638'])
    # args.add_argument('-c', '--clades', nargs='*', type=str, help='which clades to look at', default=None)
    args.add_argument('-d', '--date_range', nargs=2, type=float, help='which dates to look at', default=[1970., 2015.])
    args.add_argument('-yb', '--years_back', type=int, help='how many years of past immunity to include in fitness estimates', default=3)
    args.add_argument('-yf', '--years_forward', type=int, help='how many years into the future to predict', default=5)
    args.add_argument('-gamma', type=float, help='-1*proportion of titers that wane per year post-exposure (slope of years vs. p(titers remaining))', default= -0.15)
    args.add_argument('-sigma', type=float, help='-1*probability of protection from i conferred by each log2 titer unit against i', default=-0.25)
    args.add_argument('--plot', type=bool, help='make plots?', default=True)
    args.add_argument('--save', type=bool, help='save csv and png files?', default=False)
    args.add_argument('--name', type=str, help='analysis name')
    args.add_argument('--out_path', type=str, help='where to save csv and png files', default='./')

    args = args.parse_args()

    antigenic_fitness = AntigenicFitness(args)
    if not isinstance(antigenic_fitness.fitness, pd.DataFrame):
        print 'calculating fitness'
        antigenic_fitness.calculate_fitness()
    print 'predicting frequencies'
    antigenic_fitness.predict_rolling_frequencies()
    print 'calculating growth rates'
    antigenic_fitness.calc_growth_rates()

    if antigenic_fitness.plot:
        print 'generating plots'
        plot_fitness_v_frequency(antigenic_fitness)
        plot_rolling_frequencies(antigenic_fitness)
        plot_growth_rates(antigenic_fitness)
