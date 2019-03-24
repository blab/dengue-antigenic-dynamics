import os
import numpy as np
import random
from itertools import product
import sys

frequency_path = './simulation_freqs.csv'
titer_path = '../../data/frequencies/interserotype_Dij.tsv'
out_path = './model_performance/'

if not os.path.isdir(out_path):
	os.mkdir(out_path)

beta_vals = np.linspace(0.05,0.95,7)
gamma_vals = np.linspace(0.0,0.45,4)
sigma_vals = np.linspace(2.55,3.,4)

for (b,g,s) in product(beta_vals, gamma_vals, sigma_vals):
	name = ''.join(random.choice('0123456789abcdef') for n in xrange(30))
	cmd = 'python ../antigenic_fitness.py --mode fit --frequency_path %s --titer_path %s --out_path %s --name %s --years_back 2 --years_forward 2 --beta %f --gamma %f --sigma %f --DENV4_f0 0'%(frequency_path, titer_path, out_path, name, b,g,s)

	print cmd
	os.system# ('sbatch -t 20:00:00 --mem 1024 --wrap="%s"'%cmd)
# 	os.system('sleep 0.05')
