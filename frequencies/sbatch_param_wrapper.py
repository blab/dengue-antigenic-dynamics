import os
import numpy as np
import random
from itertools import product

antigenic_resolutions = ['all_effects', 'interserotype_effects']
clade_resolutions = ['genotype', 'serotype']
beta_vals = np.linspace(0,3,8)
gamma_vals = np.linspace(0,2,8)
sigma_vals = np.linspace(0,3,8)
d2_vals = np.linspace(0,2,8)
d3_vals = np.linspace(0,2,8)
# d4_vals = np.linspace(0,2,8)

for antigenic_resolution in antigenic_resolutions:
	for clade_resolution in clade_resolutions:
		if antigenic_resolution == 'all_effects' and clade_resolution == 'serotype':
			continue

		frequency_path = './source/southeast_asia_%s_frequencies.csv'%(clade_resolution)
		titer_path = './source/%s_Dij.csv'%antigenic_resolution
		out_path = './southeast_asia/%s/%s/'%(clade_resolution, antigenic_resolution)

		# for (b,g,s,d2,d3,d4) in product(beta_vals, gamma_vals, sigma_vals, d2_vals, d3_vals, d4_vals):
		for (b,g,s,d2,d3) in product(beta_vals, gamma_vals, sigma_vals, d2_vals, d3_vals):
			name = ''.join(random.choice('0123456789abcdef') for n in xrange(30))
			cmd = 'python ./antigenic_fitness.py --frequency_path %s --titer_path %s --out_path %s --name %s --years_back 2 --years_forward 5 --beta %f --gamma %f --sigma %f --DENV1_f0 1 --DENV2_f0 %f --DENV3_f0 %f'%(frequency_path, titer_path, out_path, name, b,g,s,d2,d3)
			# cmd = 'python ./antigenic_fitness.py --frequency_path %s --titer_path %s --out_path %s --name %s --years_back 2 --years_forward 5 --beta %f --gamma %f --sigma %f --DENV1_f0 1 --DENV2_f0 %f --DENV3_f0 %f --DENV4_f0 %f'%(frequency_path, titer_path, out_path, name, b,g,s,d2,d3,d4)

			os.system('sbatch -t 20 --wrap="%s"'%cmd)
			import sys
			sys.exit()
