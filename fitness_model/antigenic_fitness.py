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
import pickle

def safe_ratio(num, denom):
    if denom > 0:
        return float(num) / (denom)
    else:
        return np.nan

def normalize_timepoint(row):
    total = row.sum()
    try:
        assert not np.isnan(total)
        return row.map(lambda x: x / total)
    except: # np.nan or zerodivision
        return row

def normalize_all_timepoints(frequencies):
    ''' Normalize each row so that the sum of all frequencies at a single timepoint = 1'''
    if isinstance(frequencies, dict):
        frequencies = pd.DataFrame(frequencies)
        normalized_frequencies = frequencies.apply(normalize_timepoint, axis=1)
        return normalized_frequencies.to_dict()

    else:
        normalized_frequencies = frequencies.apply(normalize_timepoint, axis=1)
        return normalized_frequencies

def calc_timepoint_exposure(af, current_timepoint, frequencies=None):
    ''' for a given timepoint, estimate the relative population immunity of each clade, i, based on
    which clades (j) have circulated previously;
    how antigenically distant i is from j;
    and the relative frequency of j'''

    if frequencies is None:
        frequencies = af.actual_frequencies

    def get_Dij(j, i):
        serotype_of = lambda genotype: genotype.split('_')[0]
        try:
            return af.titers[(j, i)]
        except KeyError:
            return af.titers[(serotype_of(j), serotype_of(i))]

    def sum_over_j(i, past_timepoint):
        ''' Return a frequency-weighted sum of the probability of protection from i given prior exposure to j '''
        antigenic_distance = [ get_Dij(j, i) if i != j else 0. for j in af.clades] # Pull precomputed antigenic distance between virus i and serum j
        probability_protected = [ max(-1.*af.sigma*Dij + 1., 0.) for Dij in antigenic_distance ] # Linear transformation from antigenic distance to probability of protection from i given prior exposure to j
        # probability_protected = [ np.exp(-1.*af.sigma*Dij) for Dij in antigenic_distance]
        j_frequencies = [ frequencies[j][past_timepoint] for j in af.clades] # Pull relative frequency of each clade j at the timepoint of interest
        cumulative_immunity = sum( [ j_frequency * prob_protected for (j_frequency, prob_protected) in zip(j_frequencies, probability_protected)]) # return weighted sum
        return cumulative_immunity

    def sum_over_past_t(i):
        ''' For each timepoint, look at the past `tp_back` number of timepoints and add up the relative immunity acquired in each interval.
        Adjust for how long ago the population was exposed by assuming that immunity wanes linearly with slope gamma per year (n)'''
        tp_idx = af.timepoints.index(current_timepoint) # index of timepoint of interest
        past_timepoints = af.timepoints[tp_idx - af.tp_back : tp_idx] # previous timepoints to sum immunity over
        exposure = [ sum_over_j(i, t) for t in past_timepoints ] # total protection acquired at each past timepoint, t: sum over all clades for each past timepoint
        # waning = [ np.exp(-1.*af.gamma*(current_timepoint-t)) for t in past_timepoints]
        waning = [max(-1.*af.gamma*(current_timepoint - t) + 1., 0.) for t in past_timepoints] # proportion of protection originally acquired at time t expected to remain by the current_timepoint
        return sum( [ w*e for (w,e) in zip(waning, exposure)] ) # sum up the total waning-adjusted population immunity as of the timepoint_of_interest

    exposure = { i: sum_over_past_t(i) for i in af.clades }
    return exposure

def calc_timepoint_fitness(af, exposure):
    # full model
    fitness = { i: getattr(af, 'DENV%s_f0'%i[4]) - af.beta*exposure for i, exposure in exposure.iteritems() } ## convert from population exposure to fitness

    # antigenic null (intrinsic only)
    # fitness = { i: getattr(af, 'DENV%s_f0'%i[4]) for i, exposure in exposure.iteritems() } ## convert from population exposure to fitness

    # complete null (all equal)
    # fitness = { i: 0. for i, exposure in exposure.iteritems() } ## convert from population exposure to fitness

    return fitness

def predict_single_frequency(initial_frequency, initial_fitness, years_forward):
    return initial_frequency*np.exp(initial_fitness*years_forward)

def predict_timepoint_distant_frequency(af, current_timepoint):
    years_back, years_forward= af.years_back, af.years_forward # total
    tp_back, tp_forward = af.tp_back, af.tp_forward # total
    interval_years, interval_tp = 0.25, 1 # interval length (yrs), number of timepoints forward per interval

    x0 = current_timepoint
    x0_idx = af.timepoints.index(x0) # index of x0

    known_timepoints = af.timepoints[ :x0_idx+1] # we can pull empirical frequency and fitness values for timepoints at or before x0
    fitness = deepcopy(af.fitness.loc[known_timepoints])
    frequencies = deepcopy(af.actual_frequencies.loc[known_timepoints])

    ### Step through each intermediate timepoint, calculating known frequencies --> fitness --> predicted frequencies --> fitness
    interval_timepoints = af.timepoints[x0_idx : x0_idx+tp_forward]

    for tp_idx, numdate in enumerate(interval_timepoints, start=x0_idx):
        if tp_idx == x0_idx: ## we already know fitness for x0; just predict forward
            interval_fitness = fitness.loc[x0]
        else:
            past_timepoints = af.timepoints[tp_idx - tp_back : tp_idx] # previous timepoints to sum immunity over
            interval_exposure =calc_timepoint_exposure(af, current_timepoint=numdate, frequencies=frequencies) ## calculate fitness based on known [before x0] and predicted [after x0] frequencies
            interval_fitness = calc_timepoint_fitness(af, interval_exposure) ## convert from population exposure to fitness
            interval_fitness = pd.Series(interval_fitness, name = numdate)
            fitness = fitness.append(interval_fitness)      ## record calculated fitness vals

        interval_frequencies = predict_timepoint_close_frequency(af, current_timepoint=numdate, ## predict frequencies for each interval
                                                                final_timepoint=numdate+interval_years,
                                                                fitness=fitness, frequencies=frequencies)

        ## calculate model and null SSE contribution and add to af.model_sse / af.null_sse
        model_SE = interval_frequencies - af.actual_frequencies.loc[numdate+interval_years]
        model_SE = model_SE**2
        model_SE = model_SE.sum()


        null_interval_frequencies = predict_timepoint_close_frequency(af, current_timepoint=numdate, ## predict frequencies for each interval
                                                                final_timepoint=numdate+interval_years,
                                                                fitness='null', frequencies=frequencies)

        null_SE = null_interval_frequencies - af.actual_frequencies.loc[numdate+interval_years]
        null_SE = null_SE**2
        null_SE = null_SE.sum()

        af.model_sse += model_SE
        af.null_sse += null_SE

        frequencies = frequencies.append(interval_frequencies) ## record predicted fitness vals

    af.trajectories[current_timepoint] = frequencies.loc[x0:]
    initial_frequencies = af.actual_frequencies.loc[x0]
    interval_weighted_fitness = fitness.loc[interval_timepoints].applymap(lambda f: f*interval_years).sum()
    predicted_final_frequencies = {i: initial_frequencies[i]*np.exp(interval_weighted_fitness[i]) for i in af.clades}
    return predicted_final_frequencies

def predict_timepoint_close_frequency(af, current_timepoint, final_timepoint, fitness=None, frequencies=None):
    '''
    For each clade i, predict the frequency of i based on
    its fitness and initial frequency at time t-years_forward
    '''
    years_forward = final_timepoint - current_timepoint

    if fitness is None:
        fitness = af.fitness
    elif type(fitness) == str and fitness == 'null':
        assert years_forward <= 1.0
    if frequencies is None:
        frequencies = af.actual_frequencies

    if years_forward <= 1.0:
        initial_frequencies = frequencies.loc[current_timepoint]
        if type(fitness) == str and fitness == 'null':
            predicted_frequencies = { i:predict_single_frequency(initial_frequencies[i],
                                                            0.,
                                                            final_timepoint - current_timepoint)
                                        for i in af.clades }
        else:
            initial_fitnesses = fitness.loc[current_timepoint]
            predicted_frequencies = { i : predict_single_frequency(initial_frequencies[i],
                                                            initial_fitnesses[i],
                                                            final_timepoint - current_timepoint)
                                    for i in af.clades }
    else:
        predicted_frequencies = predict_timepoint_distant_frequency(af, current_timepoint)

    predicted_frequencies = pd.Series(predicted_frequencies, name=final_timepoint)
    predicted_frequencies = normalize_timepoint(predicted_frequencies)
    return predicted_frequencies


def calc_timepoint_growth_rates(af,i,predicted=True):

    initial_timepoints = af.timepoints[af.tp_back: -1*af.tp_forward]
    initial_frequencies = af.actual_frequencies[i][initial_timepoints]

    final_timepoints = af.timepoints[af.tp_forward+af.tp_back:]
    if predicted==True:
        final_frequencies = af.predicted_frequencies[i][final_timepoints]
    else:
        final_frequencies = af.actual_frequencies[i][final_timepoints]

    time_intervals = [ str(f)+'/'+str(i) for (i, f) in zip(initial_timepoints, final_timepoints)]
    growth_rates = [ safe_ratio(f, i) for (i,f) in zip(initial_frequencies, final_frequencies)]

    return pd.Series(growth_rates, index=time_intervals)

def predict_trajectories(af, initial_timepoint):
    '''
    return predicted frequencies at all timepoints between t and t+dt
    '''

    if af.years_forward > 1:
        if initial_timepoint not in af.trajectories:
            predict_timepoint_distant_frequency(af, initial_timepoint)

        return af.trajectories[initial_timepoint]
    else:
        raise ValueError, 'Oops. Currently only doing trajectories for dt > 1'
    # dt_values = [ (1./af.tppy)*dt for dt in range(1, af.tp_forward+1)] # fraction of year per timepoint * number of timepoints forward
    # predicted_timepoints = [ initial_timepoint + dt for dt in dt_values ]
    # initial_frequency = af.actual_frequencies[i][initial_timepoint]
    # initial_fitness = af.fitness[i][initial_timepoint]
    #
    # predicted_trajectory = [ predict_single_frequency(initial_frequency, initial_fitness, dt) for dt in dt_values ]
    #
    # return pd.Series(predicted_trajectory, index=predicted_timepoints, name=i)

def calc_delta_sse(af):
    ''' How much better were our predictions than the null model for time t+N? '''

    return af.null_sse - af.model_sse

    # def calc_sse(af, valid_clades, starting_timepoint, null):
    #     sse = 0.
    #     for clade in valid_clades:
    #         actual_frequency = af.actual_frequencies[clade][starting_timepoint + af.years_forward]
    #
    #         if null == True:
    #             null_frequency = af.actual_frequencies[clade][starting_timepoint]
    #             squared_error = (actual_frequency - null_frequency)**2
    #         else:
    #             predicted_frequency = af.predicted_frequencies[clade][starting_timepoint + af.years_forward]
    #             squared_error = (actual_frequency - predicted_frequency)**2
    #
    #         if not np.isnan(squared_error):
    #             sse += squared_error
    #     return sse
    #
    # d_sse = 0.
    # for starting_timepoint in af.timepoints[af.tp_back: -1*af.tp_forward]:
    #     valid_clades = [c for c in af.clades if af.actual_frequencies[c][starting_timepoint] >= 0.1 ]
    #     model_sse = calc_sse(af, valid_clades, starting_timepoint, null=False)
    #     null_sse = calc_sse(af, valid_clades, starting_timepoint, null=True)
    #     d_sse +=  null_sse - model_sse
    # return d_sse

def calc_accuracy(af):

    def clade_accuracy(i):
        mask = (~np.isnan(af.actual_growth_rates[i]) & (~np.isnan(af.predicted_growth_rates[i])))
        actual, predicted = af.actual_growth_rates[i][mask], af.predicted_growth_rates[i][mask]

        correct_predictions = 0.
        for (a,p) in zip(actual, predicted):
            if (a >= 1. and p>=1.) or (a < 1. and p < 1.):
                correct_predictions += 1.
        return {'n_correct': correct_predictions, 'n_total': float(len(actual))}

    n_correct = 0.
    n_total = 0.
    for i in af.clades:
        accuracy = clade_accuracy(i)
        n_correct += accuracy['n_correct']
        n_total += accuracy['n_total']
    return n_correct / n_total

def calc_model_performance(af):

    def remove_nan(actual, predicted):
        actual = actual.loc[actual.index.isin(predicted.index.values)]
        assert predicted.columns.tolist() == actual.columns.tolist()
        assert actual.index.tolist() == predicted.index.tolist()
        actual, predicted = actual.values.flatten(), predicted.values.flatten()
        mask = (~np.isnan(actual)) & (~np.isnan(predicted))
        return actual[mask], predicted[mask]

    actual_freq, predicted_freq = remove_nan(af.actual_frequencies, af.predicted_frequencies)
    actual_growth, predicted_growth = remove_nan(af.actual_growth_rates, af.predicted_growth_rates)

    performance = {
    'pearson_r2': stats.linregress(actual_growth, predicted_growth)[2]**2,
    'spearman_r': stats.spearmanr(actual_growth, predicted_growth)[0],
    'abs_error': sum([abs(a - p) for (a,p) in zip(actual_freq, predicted_freq)]) / float(len(actual_freq)),
    'accuracy': calc_accuracy(af),
    'delta_sse': calc_delta_sse(af)}

    return performance

class AntigenicFitness():
    def __init__(self, args):

        for k,v in vars(args).items():
            setattr(self, k, v) # copy over cmd line args

        # actual (observed) actual_frequencies
        self.actual_frequencies = pd.read_csv(self.frequency_path, index_col=0) # pd.DataFrame(index=timepoints, columns=clades, values=relative actual_frequencies)
        self.clades = self.actual_frequencies.columns.values
        self.actual_frequencies = normalize_all_timepoints(self.actual_frequencies[self.clades]) # restrict to clades of interest, normalize
        self.actual_frequencies = self.actual_frequencies.loc[(self.actual_frequencies.index >= self.date_range[0]) & (self.actual_frequencies.index <= self.date_range[1])] # restrict to timepoints of interest
        self.timepoints = self.actual_frequencies.index.tolist()
        n_years = int(self.timepoints[-1]) - int(self.timepoints[0]) # number of years in the actual_frequencies dataset
        self.tppy = int(len(self.timepoints)/n_years) # timepoints per year
        self.tp_forward = self.tppy*self.years_forward # number of timepoints forward
        self.tp_back = self.tppy*self.years_back # number of timepoints back

        self.noisy_predictions_mask = self.actual_frequencies < 0.05 # log which initial values were low
        # keep track of which predictions will be made based on low initial values
        self.noisy_predictions_mask.index = self.noisy_predictions_mask.index.map(lambda x: x+self.years_forward)

        self.fitness = None

        self.model_sse = 0.
        self.null_sse = 0.

        self.trajectories = {}

        # load pre-computed antigenic distances between clades
        self.titers = {(serum,virus): Dij for (serum,virus), Dij in
                        pd.Series.from_csv(args.titer_path, header=None,index_col=[0,1], sep='\t').to_dict().items()}

        if self.save == True:
            assert self.name, 'ERROR: Please provide an analysis name if you wish to save output'

    def calculate_fitness(self):
        ''' fitness = 1.-population exposure'''
        valid_timepoints = self.timepoints[self.tp_back:]
        exposure = {t: calc_timepoint_exposure(self, t) for t in valid_timepoints}
        fitness = { t: calc_timepoint_fitness(self, e) for t, e in exposure.iteritems()}
        self.fitness = pd.DataFrame.from_dict(fitness, orient='index')
        if self.save:
            self.fitness.to_csv(self.out_path+self.name+'_fitness.csv')

    def predict_frequencies(self):
        '''
        Making a "rolling prediction" of frequency for each clade.
        Normalize these predicted actual_frequencies so that they sum to 1. at each timepoint,
        then mask out predictions based on noisy initial actual_frequencies.

        Returns pd.DataFrame(index = t+dt, columns=clades,
                            values = predicted frequency at time t+dt, based on frequency & fitness at time t )
        Xi(t+dt) = Xi(t) * e^( Fi(t) * dt), where Xi is frequency and Fi is fitness of i
        '''

        initial_timepoints = self.timepoints[self.tp_back: -1*self.tp_forward]
        predicted_timepoints = self.timepoints[self.tp_forward+self.tp_back:]

        all_predicted_frequencies = pd.DataFrame.from_dict({predicted_timepoint : predict_timepoint_close_frequency(self, initial_timepoint, predicted_timepoint)
                                    for (initial_timepoint, predicted_timepoint) in zip(initial_timepoints, predicted_timepoints)}, orient='index')
        self.predicted_frequencies = all_predicted_frequencies[~self.noisy_predictions_mask] # keep only predictions based on initial actual_frequencies at >0.1

        if self.save:
            self.predicted_frequencies.to_csv(self.out_path+self.name+'_predicted_freqs.csv')

    def simulate_forward(self):
        '''
        Making a "rolling prediction" of frequency for each clade.
        Normalize these predicted frequencies so that they sum to 1. at each timepoint,
        *** and then replace the empirical values with the predicted values before moving to the next timepoint ***

        Returns pd.DataFrame(index = t+dt, columns=clades,
                            values = simulated frequency at time t+dt, based on [actual or simulated*] frequency & fitness at time t )
                            * for timepoints before tp_back + tp_forward, based on actual initial frequencies; for timepoints after tp_back + tp_forward, based on simulated initial frequencies
        Xi(t+dt) = Xi(t) * e^( Fi(t) * dt), where Xi is frequency and Fi is fitness of i
        '''

        initial_timepoints = self.timepoints[self.tp_back: -1*self.tp_forward]
        predicted_timepoints = self.timepoints[self.tp_forward+self.tp_back:]

        all_predicted_frequencies = {}
        ## To simulate forward, we just have to replace the empirical data with the predicted data before moving on to the next timepoint
        for initial_timepoint, final_timepoint in zip(initial_timepoints, predicted_timepoints):
            ## Re-estimate fitness based on simulated frequencies.
            ## For timepoints before tp_back + tp_forward, these values should be unchanged. For timepoints after tp_back+tp_forward, these values will be updated.
            t0_exposure = calc_timepoint_exposure(self, initial_timepoint)
            t0_fitness = calc_timepoint_fitness(self, t0_exposure)
            self.fitness.loc[initial_timepoint] = t0_fitness

            ## Predict t+dt based on updated fitnesses
            predicted_frequencies = predict_timepoint_close_frequency(self, initial_timepoint, final_timepoint)
            ## Replace the empirical frequencies at t+dt with the predicted ones
            self.actual_frequencies.loc[final_timepoint] = predicted_frequencies

            ''' I think the predicted frequencies now become basically meaningless / can be discarded....? '''
            # all_predicted_frequencies[final_timepoint] = predicted_frequencies
        # all_predicted_frequencies = pd.DataFrame.from_dict(all_predicted_frequencies, orient='index')

        ''' Re-mask which predictions were based on low initial values? I think this is wrong / needs to be fixed. For now, leave it out. '''
        # self.noisy_predictions_mask = self.actual_frequencies < 0.05 # log which initial values were low
        # self.noisy_predictions_mask.index = self.noisy_predictions_mask.index.map(lambda x: x+self.years_forward)
        # self.predicted_frequencies = all_predicted_frequencies[~self.noisy_predictions_mask] # keep only predictions based on initial actual_frequencies at >0.1

        if self.save:
            self.fitness.to_csv(self.out_path+self.name+'_fitness.csv')
            self.actual_frequencies.to_csv(self.out_path+self.name+'_simulated_freqs.csv')

    def calc_growth_rates(self):
        self.predicted_growth_rates = pd.DataFrame({i:calc_timepoint_growth_rates(self, i, predicted=True) for i in self.clades})
        self.actual_growth_rates = pd.DataFrame({i:calc_timepoint_growth_rates(self, i, predicted=False) for i in self.clades})
        if self.save:
            self.predicted_growth_rates.to_csv(self.out_path+self.name+'_predicted_growth_rates.csv')
            self.actual_growth_rates.to_csv(self.out_path+self.name+'_actual_growth_rates.csv')

    # def predict_trajectories(self,initial_timepoint):
    #     '''
    #     Predict the frequency of all clades at each time interval between t and t+years_forward,
    #     based on their initial fitnesses and actual_frequencies at time t.
    #
    #     Normalize these predicted actual_frequencies so that they sum to 1. at each timepoint.
    #     '''
    #
    #     all_trajectories = pd.DataFrame({ i : predict_clade_trajectory(self, i, initial_timepoint)
    #                        for i in self.clades})
    #     self.trajectories[initial_timepoint] = normalize_all_timepoints(all_trajectories)

def plot_fitness_v_frequency(af):
    sns.set_palette('tab20', n_colors=20)

    fig, axes = plt.subplots(2,1,figsize=(8,6), sharex=True)

    fitness=af.fitness
    actual_frequencies = af.actual_frequencies

    for clade in fitness.columns.values:
        axes[0].plot(fitness[clade].index.values, fitness[clade], linestyle='--', label='%s fitness'%clade)
        axes[1].plot(actual_frequencies[clade].index.values, actual_frequencies[clade], linestyle='-', label='%s frequency'%clade)

    axes[0].set_title('Clade fitness')
    axes[1].set_title('Observed clade frequency')
    plt.tight_layout()

    if af.save:
        plt.savefig(af.out_path+af.name+'_fitness.png', dpi=300)
    else:
        plt.show()
    plt.clf()
    plt.close()

def plot_frequencies(af):
    sns.set_palette('tab20', n_colors=20)

    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(12,8), sharex=True)
    for clade, predicted_freqs in af.predicted_frequencies.iteritems():
        date_min, date_max = predicted_freqs.index.min(), predicted_freqs.index.max()
        axes[0].plot(predicted_freqs.index.values,predicted_freqs.values, linestyle='--', label='Predicted %s actual_frequencies'%clade)

        actual_frequencies = af.actual_frequencies[clade][date_min:date_max]
        axes[1].plot(actual_frequencies.index.values, actual_frequencies.values, label='Actual %s actual_frequencies'%clade)

    axes[0].set_title('Predicted clade actual_frequencies')
    axes[1].set_title('Observed clade actual_frequencies')

    plt.tight_layout()
    if af.save:
        plt.savefig(af.out_path+af.name+'_frequencies.png', dpi=300)
    else:
        plt.show()
    plt.clf()
    plt.close()

def plot_growth_rates(af):
    '''
    For the actual and predicted actual_frequencies, find where both values are non-null and > 0.1
    Plot actual vs. predicted
    '''
    sns.set_palette('Set2', n_colors=10)
    actual, predicted = af.actual_growth_rates, af.predicted_growth_rates
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

    if af.save:
        plt.savefig(af.out_path+af.name+'_growth_rates.png', dpi=300)
    else:
        plt.show()
    plt.clf()
    plt.close()

def plot_trajectory(af, trajectory, ax=None):

    if not ax:
        fig, ax = plt.subplots(figsize=(8,5))
    cmap = sns.color_palette('Set2', len(trajectory.columns.values))
    colors = { clade: cmap[i] for i, clade in enumerate(trajectory.columns.values)}

    initial_timepoint = trajectory.index.values[0]
    initial_tp_idx = af.timepoints.index(initial_timepoint)
    past_timepoints = af.timepoints[initial_tp_idx-af.tp_back : initial_tp_idx+1]

    for clade, clade_trajectory in trajectory.iteritems():
        actual_vals = af.actual_frequencies[clade][past_timepoints]
        ax.plot(actual_vals.index.values, actual_vals.values, c=colors[clade], label=clade)
        ax.plot(clade_trajectory.index.values, clade_trajectory.values, c=colors[clade], linestyle='--')
        # ax.plot(actual_vals.index.values, af.fitness[clade][actual_vals.index.values], c=colors[clade], linestyle=':')
    ax.set_xlim(past_timepoints[0], trajectory.index.values[-1])
    ax.set_ylim(0,1)
    # ax.set_ylim(-0.05,1.05)

    plt.legend()
    plt.tight_layout()
    # if af.save:
    #     plt.savefig(af.out_path + af.name + '%.1f_trajectory.png'%initial_timepoint, dpi=300, bbox_inches='tight')
    # else:
    #     plt.show()
#
def plot_trajectory_multiples(af, starting_timepoints=None, n_clades_per_plot=4):
    sns.set(style='whitegrid', font_scale=0.8)

    if starting_timepoints == None:
        starting_timepoints = af.timepoints[af.tp_forward+af.tp_back::af.tp_forward]
        starting_timepoints.append(1984.)
        starting_timepoints = sorted(starting_timepoints)

    ncols = len(starting_timepoints)

    if len(af.clades) > n_clades_per_plot:
        nrows = int(ceil(len(af.clades)/n_clades_per_plot))
        fig, axes = plt.subplots(nrows, ncols, figsize=(3*ncols, 2*nrows),sharex=False, sharey=True)
        clade_sets = [af.clades[i:i + n_clades_per_plot] for i in xrange(0, len(af.clades), n_clades_per_plot)]
        zipped = zip(clade_sets, axes)
    else:
        fig, axes = plt.subplots(1, ncols, figsize=(3*ncols, 2),sharex=False, sharey=True)
        zipped = [(af.clades, axes)]

    for clade_set, row in zipped:
        for tp, ax in zip(starting_timepoints, row):
            trajectory = predict_trajectories(af, tp)
            plot_trajectory(af, trajectory, ax)

    ax.legend()
    plt.tight_layout()
    if af.save:
        plt.savefig(af.out_path+af.name+'_trajectories.png', bbox_inches='tight', dpi=300)
    else:
        plt.show()
    plt.clf()
    plt.close()

if __name__=="__main__":
    sns.set(style='whitegrid')#, font_scale=1.5)

    args = argparse.ArgumentParser()
    args.add_argument('--frequency_path', help='actual_frequencies csv', default='../data/frequencies/seasia_serotype_frequencies.csv')
    args.add_argument('--titer_path', help='pairwise dTiters csv', default='../data/frequencies/fulltree_Dij.tsv')
    args.add_argument('--date_range', nargs=2, type=float, help='which dates to look at', default=[1970., 2015.])
    args.add_argument('--years_back', type=int, help='how many years of past immunity to include in fitness estimates', default=2)
    args.add_argument('--years_forward', type=int, help='how many years into the future to predict', default=2)
    args.add_argument('--gamma', type=float, help='Value or value range for the proportion of titers that wane per year post-exposure (slope of years vs. p(titers remaining))', default=0.86)
    args.add_argument('--sigma', type=float, help='Value or value range for -1*probability of protection from i conferred by each log2 titer unit against i', default=0.43)
    args.add_argument('--beta', type=float, help='Value or value range for beta. antigenic fitness = -1.*beta*population_immunity', default=2.57)
    args.add_argument('--DENV1_f0', type=float, help='Relative intrinsic fitness value for DENV1', default = 2.0)
    args.add_argument('--DENV2_f0', type=float, help='Relative intrinsic fitness value for DENV2', default = 2.0)
    args.add_argument('--DENV3_f0', type=float, help='Relative intrinsic fitness value for DENV3', default = 1.0)
    args.add_argument('--DENV4_f0', type=float, help='Relative intrinsic fitness value for DENV4', default = 0.)
    args.add_argument('--trajectory', type=float, nargs='*', help='timepoint(s) to compute trajectories for')
    args.add_argument('--plot', help='make plots?', action='store_true')
    args.add_argument('--save', help='save csv and png files?', action='store_true')
    args.add_argument('--name', type=str, help='analysis name')
    args.add_argument('--out_path', type=str, help='where to save csv and png files', default='./')
    args.add_argument('--mode', type=str, choices=['fit', 'run', 'simulate'], help='Fit parameters, simulate, or run model?', default='run')
    args = args.parse_args()

    if args.mode == 'simulate':
        antigenic_fitness = AntigenicFitness(args)
        antigenic_fitness.calculate_fitness()
        antigenic_fitness.simulate_forward()
        # antigenic_fitness.calculate_fitness()
        # antigenic_fitness.predict_frequencies()
        # antigenic_fitness.calc_growth_rates()


    elif args.mode == 'fit':

        d1_vals = np.linspace(0,6,7)
        d2_vals = np.linspace(0,6,7)
        d3_vals = np.linspace(0,6,7)

        output = []
        for (d1,d2,d3) in product(d1_vals, d2_vals, d3_vals):
            args = deepcopy(args)
            setattr(args, 'DENV1_f0', d1)
            setattr(args, 'DENV2_f0', d2)
            setattr(args, 'DENV3_f0', d3)
            antigenic_fitness = AntigenicFitness(args)
            if not isinstance(antigenic_fitness.fitness, pd.DataFrame):
                antigenic_fitness.calculate_fitness()
            antigenic_fitness.predict_frequencies()
            antigenic_fitness.calc_growth_rates()

            model_performance = calc_model_performance(antigenic_fitness)
            model_performance.update({ 'beta': args.beta, 'sigma': args.sigma, 'gamma': args.gamma, 'DENV1_f0': args.DENV1_f0, 'DENV2_f0': args.DENV2_f0, 'DENV3_f0': args.DENV3_f0})
            output.append(model_performance)


        output = pd.DataFrame(output)
        output = output.reindex(columns=sorted(output.columns.values))
        output.to_csv(args.out_path+args.name+'.csv')

    else:
        antigenic_fitness = AntigenicFitness(args)
        antigenic_fitness.calculate_fitness()
        antigenic_fitness.predict_frequencies()
        antigenic_fitness.calc_growth_rates()
        model_performance = calc_model_performance(antigenic_fitness)
        print model_performance
        # plot_trajectory_multiples(antigenic_fitness, n_clades_per_plot=4)

        if args.trajectory:
            assert args.save or args.plot, 'only bother computing trajectories if we are going to save and/or plot them'
            for t in args.trajectory:
                closest_timepoint = sorted(list(antigenic_fitness.timepoints), key = lambda tp: abs(t - tp))[0]
                trajectory = predict_trajectories(antigenic_fitness, closest_timepoint)
                if args.save:
                    trajectory.to_csv(args.out_path + '%s_%.1f_trajectory.csv'%(args.name, t))
                if args.plot:
                    plot_trajectory(antigenic_fitness, trajectory)
