from glob import glob
import pandas as pd
import numpy as np
import os

parameters = ['DENV1_f0', 'DENV2_f0', 'DENV3_f0', 'abs_error', 'accuracy', 'beta', 'delta_sse', 'gamma', 'information_gain', 'pearson_r2', 'sigma', 'spearman_r']
all_values = []

files = glob('./*.csv')
all_dfs = [ pd.read_csv(f) for f in files if 'model_performance' not in f ]

collated = pd.concat(all_dfs)
collated.to_csv('./model_performance.csv')

for f in files:
    if 'model_performance' not in f:
        os.remove(f)
