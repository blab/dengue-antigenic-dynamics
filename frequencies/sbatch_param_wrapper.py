import os

antigenic_resolutions = ['all_effects', 'interserotype_effects']
clade_resolutions = ['genotype', 'serotype']
dt_vals = [1,5]

for antigenic_resolution in antigenic_resolutions:
	for clade_resolution in clade_resolutions:
		if antigenic_resolution == 'all_effects' and clade_resolution == 'serotype':
			continue
		for dt in dt_vals:
			frequency_path = './source/southeast_asia_%s_frequencies.csv'%(clade_resolution)
			titer_path = './source/%s_Dij.csv'%antigenic_resolution
			out_path = './southeast_asia/%s/%s/'%(clade_resolution, antigenic_resolution)
			incidence_corrected_name = 'seasia_%s_%s_dt%d_ic'%(clade_resolution, antigenic_resolution, dt)
			uncorrected_name = 'seasia_%s_%s_dt%d'%(clade_resolution, antigenic_resolution, dt)

			cmd = 'python ./antigenic_fitness.py --frequency_path %s --titer_path %s --out_path %s --name %s --years_back 2 --years_forward %d --beta 0 3 --gamma 0 2 --sigma 0 3 --DENV1_f0 1 --DENV2_f0 0 2 --DENV3_f0 0 2 --DENV4_f0 0 2 --n_param_vals 8 --save'%(frequency_path, titer_path, out_path, uncorrected_name, dt)
			ic_cmd = 'python ./antigenic_fitness.py --frequency_path %s --titer_path %s --out_path %s --name %s --years_back 2 --years_forward %d --beta 0 3 --gamma 0 2 --sigma 0 3 --DENV1_f0 1 --DENV2_f0 0 2 --DENV3_f0 0 2 --DENV4_f0 0 2 --n_param_vals 8 --save --incidence_correction'%(frequency_path, titer_path, out_path, incidence_corrected_name, dt)
			print cmd
			print ic_cmd
			os.system("sbatch --wrap='%s'"%cmd)
			os.system("sbatch --wrap='%s'"%ic_cmd)
