import os
os.system("module load python2")
os.system("pip install --user --upgrade matplotlib") # hack...... :)

# regions = ['southeast_asia', 'south_america']
region = 'southeast_asia'
N = [1,3,5]
# N=[3]
antigenic_resolutions = ['all_effects']#, 'interserotype_effects']
clade_resolutions = ['genotype', 'serotype']

for n in N:
	for antigenic_resolution in antigenic_resolutions:
		for clade_resolution in clade_resolutions:
			frequency_path = './%s/%s/%s_%s_frequencies.csv'%(region, clade_resolution, region, clade_resolution)
			titer_path = './%s_Dij.csv'%antigenic_resolution
			out_path = './%s/%s/%s/'%(region, clade_resolution, antigenic_resolution)
			name = '%s_%s_%s_N%d'%(region, clade_resolution, antigenic_resolution, n)
			cmd = 'python ./antigenic_fitness.py --frequency_path %s --titer_path %s --out_path %s --name %s --years_back %d --years_forward 1 --beta 0 10 --gamma 0 2 --sigma 0 5 --n_param_vals 10 --plot --save'%(frequency_path, titer_path, out_path, name, n)
			print cmd
			os.system("sbatch --wrap='%s'"%cmd)

# python ./antigenic_fitness.py --frequency_path ./flu/flu_frequencies.csv --titer_path ./flu/flu_titers.csv --out_path ./flu/ --name flu --years_back 1 --years_forward 1 --beta 1 15 --gamma -5 0 --sigma -5 0 --n_param_vals 10
