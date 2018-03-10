import os

regions = ['southeast_asia', 'south_america']
antigenic_resolutions = ['all_effects', 'interserotype_effects']
clade_resolutions = ['genotype', 'serotype']

for region in regions:
    for antigenic_resolution in antigenic_resolutions:
        for clade_resolution in clade_resolutions:
            frequency_path = './%s/%s/%s_%s_frequencies.csv'%(region, clade_resolution, region, clade_resolution)
            titer_path = './%s_Dij.csv'%antigenic_resolution
            out_path = './%s/%s/'%(region, clade_resolution)
            name = '%s_%s_%s'%(region, clade_resolution, antigenic_resolution)
            cmd = 'python ./antigenic_fitness.py --frequency_path %s --titer_path %s --out_path %s --name %s --years_back 1 --years_forward 1 --beta 1 15 --gamma -2 0 --sigma -5 0 --n_param_vals 10 '%(frequency_path, titer_path, out_path, name)
            print cmd
            os.system(cmd)
