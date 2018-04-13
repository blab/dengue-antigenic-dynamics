import os

antigenic_resolutions = ['all_effects', 'interserotype_effects']
clade_resolutions = ['genotype', 'serotype']

for antigenic_resolution in antigenic_resolutions:
	for clade_resolution in clade_resolutions:
		if antigenic_resolution == 'all_effects' and clade_resolution == 'serotype':
			continue
		frequency_path = './southeast_asia/%s/southeast_asia_%s_frequencies.csv'%(clade_resolution, clade_resolution)
		titer_path = './%s_Dij.csv'%antigenic_resolution
		out_path = './southeast_asia/%s/%s/'%(clade_resolution, antigenic_resolution)
		name = 'seasia_%s_%s'%(clade_resolution, antigenic_resolution)
		cmd = 'python ./antigenic_fitness.py --frequency_path %s --titer_path %s --out_path %s --name %s --years_back 2 --years_forward 1 --beta 0 3 --gamma 0 2 --sigma 0 3 --n_param_vals 12 --metric sse --plot --save'%(frequency_path, titer_path, out_path, name)
		print cmd
		os.system("sbatch --wrap='%s'"%cmd)
