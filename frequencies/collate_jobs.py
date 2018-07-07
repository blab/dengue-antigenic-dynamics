from glob import glob
import pandas as pd
import numpy as np
import os

parameters = ['DENV2_f0', 'DENV3_f0', 'DENV4_f0', 'abs_error', 'accuracy', 'beta', 'delta_sse', 'gamma', 'information_gain', 'pearson_r2', 'sigma', 'spearman_r']
all_values = []

files = glob('./*.csv')

for file in files:
    if 'model_performance' not in file:
        with open(file, 'r') as f:
            for line in f:
                vals = line.strip().split(',')
                vals = { param : v for (param, v) in zip(parameters, vals)}
                all_values.append(vals)

if len(all_values):
    collated = pd.DataFrame(all_values)
    collated.to_csv('model_performance.csv')

for f in files:
    if 'model_performance' not in f:
        os.remove(f)
