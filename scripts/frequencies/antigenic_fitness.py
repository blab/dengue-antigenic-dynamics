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

def sum_over_j(antigenic_fitness, i, timepoint, proportion_remaining):
    '''
    Look at all a single time point.
    At that timepoint, look at all clades, j, that are cocirculating with i.
    Pull the titers between i and j (D_ij) and adjust for waning immunity (proportion_remaining).
    Use this to calculate a frequency-weighted estimate sum of the probability of protection against i given prior exposure to j.
    '''

    frequency_weighted_protection = 0.

    for j in antigenic_fitness.clades:
        if i==j:
            init_titers = 0.
        else:
            init_titers = antigenic_fitness.titers[tuple(sorted([i,j]))]
        remaining_titers = proportion_remaining * init_titers
        probability_protected = max(antigenic_fitness.sigma*remaining_titers + 1., 0.)
        j_frequency = antigenic_fitness.frequencies[j][timepoint]

        frequency_weighted_protection += probability_protected*j_frequency

    return frequency_weighted_protection

def sum_over_past_t(antigenic_fitness, i, timepoint_of_interest):
    '''
    For a given time point of interest, look at what the population has acquired protection to
    over the past `tp_back` timepoints.
    For each previous timepoint, sum_over_j to calculate the accumulated immunity.
    '''

    def waning(gamma, n):
        ''' Assume immunity wanes linearly with slope gamma per year (n)'''
        return max(gamma*n + 1., 0.)

    tp_idx = antigenic_fitness.timepoints.index(timepoint_of_interest) # index of timepoint of interest
    t_to_sum = antigenic_fitness.timepoints[tp_idx - antigenic_fitness.tp_back : tp_idx] # previous timepoints to sum immunity over

    # proportion of titers acquired in each interval expected to remain by timepoint of interest
    waning_over_time = [ waning(antigenic_fitness.gamma, t - timepoint_of_interest) for t in t_to_sum]

    # proportion of the population that acquired protection in each interval
    protection_over_time = [ sum_over_j(antigenic_fitness, i, t, p_remaining)
                          for (t, p_remaining) in zip(t_to_sum, waning_over_time) ]

    accumulated_protection = sum(protection_over_time)/float(len(protection_over_time))
    return accumulated_protection

def population_exposure(antigenic_fitness, i):
    ''' estimate the proportion of the population that is immune to i at the beginning of season t
    '''

    valid_timepoints = antigenic_fitness.timepoints[antigenic_fitness.tp_back:]

    # if antigenic_resolution == 'null':
    #     return { t: 0. for t in valid_timepoints }
    # else:
    population_exposure = { t: sum_over_past_t(antigenic_fitness, i, t) for t in valid_timepoints}

    return population_exposure

class AntigenicFitness():
    def __init__(self, args):

        self.clades = args.clades
        self.date_range = args.date_range

        self.frequencies = pd.read_csv(args.frequency_path, index_col=0)
        self.frequencies = normalize_frequencies_by_timepoint(self.frequencies[self.clades]) # normalize the actual frequencies
        self.frequencies = self.frequencies.loc[(self.frequencies.index >= self.date_range[0]) & (self.frequencies.index <= self.date_range[1])]

        self.years_forward = args.years_forward
        self.years_back = args.years_back
        self.timepoints = self.frequencies.index.tolist()
        n_years = int(self.timepoints[-1]) - int(self.timepoints[0]) # number of years in the frequencies dataset
        self.tppy = int(len(self.timepoints)/n_years) # timepoints per year

        self.tp_back = self.years_back*self.tppy # number of timepoints to sum over
        self.tp_forward = self.years_forward*self.tppy # number of timepoints forward to predict

        self.titers = {(str(k1), str(k2)):v for (k1,k2),v in pd.Series.from_csv(args.dTiters_path, header=None,index_col=[0,1]).to_dict().items()}

        self.sigma=args.sigma
        self.gamma=args.gamma

        self.fitness = None

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

def plot_fitness_v_frequency(antigenic_fitness):

    fig, axes = plt.subplots(2,1,figsize=(8,6))

    fitness=antigenic_fitness.fitness
    frequencies = antigenic_fitness.frequencies

    for clade in fitness.columns.values:
        axes[0].plot(fitness[clade].index.values, fitness[clade], linestyle='--', label='%s fitness'%clade)
        axes[1].plot(frequencies[clade].index.values, frequencies[clade], linestyle='-', label='%s frequency'%clade)

    axes[0].set_title('Clade fitness')
    axes[1].set_title('Observed clade frequency')
    # plt.sca(axes[0])
    # plt.legend(loc=(1,1), ncol=2)
    # plt.sca(axes[1])
    # plt.legend(loc=(1,1), ncol=2)
    plt.tight_layout()

    if antigenic_fitness.save:
        plt.savefig(antigenic_fitness.out_path+antigenic_fitness.name+'_fitness.png', dpi=300)
    else:
        plt.show()

if __name__=="__main__":
    sns.set(style='whitegrid', font_scale=1.5)
    sns.set_palette('tab20', n_colors=20)

    args = argparse.ArgumentParser()
    args.add_argument('-f', '--frequency_path', help='frequencies csv', default='../../data/titer-model/frequencies/southeast_asia_clade_frequencies.csv')
    args.add_argument('-t', '--dTiters_path', help='pairwise dTiters csv', default='../../data/titer-model/frequencies/clade_dtiters.csv')
    args.add_argument('-c', '--clades', nargs='*', type=str, help='which clades to look at', default=['2185', '2589', '2238', '2596', '1460', '1393', '1587', '1455', '975', '979', '1089', '33', '497', '117', '543', '4', '638'])
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
    antigenic_fitness.calculate_fitness()

    if antigenic_fitness.plot:
        plot_fitness_v_frequency(antigenic_fitness)
