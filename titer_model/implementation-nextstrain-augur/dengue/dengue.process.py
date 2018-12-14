from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../'))
import base.process
from base.process import process
import argparse
import numpy as np
from dengue_titers import titer_model, titer_export ## Set up and parameterize the titer model separately for tidiness

##### Define references and metadata #####
sanofi_vaccine_strains = {
	'denv1': 'DENV1/THAILAND/PUO359/1980',
	'denv2': 'DENV2/THAILAND/PUO218/1980',
	'denv3': 'DENV3/THAILAND/PAH88188/1988',
	'denv4': 'DENV4/INDONESIA/S1228/1978',
	'all': None}

regions = ['africa', 'europe', 'north_america', 'china', 'south_asia',
			'japan_korea', 'south_pacific', 'oceania', 'south_america',
			'southeast_asia', 'west_asia']

##### Parse args, set up config #####
def collect_args():
	"""Returns a dengue-specific argument parser."""
	parser = base.process.collect_args()

	# parser.add_argument('-j', '--json', default=None, nargs='+', type=str, help="Accepts path to prepared JSON(s); overrides -s argument")
	parser.add_argument('-s', '--serotypes', default=["all"], nargs='+', type=str, choices=['denv1', 'denv2', 'denv3', 'denv4', 'all', 'multiple'],
	help="Look for prepared JSON(s) like ./prepared/dengue_SEROTYPE.json; 'multiple' will run all five builds. Default='multiple'")
	parser.add_argument('--no_mut_freqs', default=True, action='store_true', help="skip mutation frequencies")
	parser.add_argument('--no_tree_freqs', default=False, action='store_true', help="skip tree (clade) frequencies")
	parser.add_argument('--no_titers', default=False, action='store_true', help="skip titer models")
	parser.add_argument('--titer_model', default='full_tree', choices=['full_tree', 'interserotype'], type=str)
	parser.add_argument('--output', default=None)
	parser.add_argument('--lam_drop', default=1.0, type=float)
	parser.set_defaults(json = './prepared/dengue_config.json')
	return parser


def make_config (prepared_json, args):
	"""
	Configure your analysis here.
	Parsed as a function to enable running multiple builds with one cmd.
	"""
	if args.output is None:
		output = {"auspice": "../../%s_model_output/"%args.titer_model, "data": "./processed/"}
	else:
		output = {"auspice": args.output, "data": "./processed/"}

	return {
		"dir": "dengue",
		"in": prepared_json,
		"output": output,
		"geo_inference": ['region'], # what traits to perform this on; don't run country (too many demes, too few sequences per deme to be reliable)
		"auspice": { ## settings for auspice JSON export
			"extra_attr": ['serum', 'clade', 'dTiter_sanofi'], # keys from tree.tree.clade['attr'] to include in export
			"color_options": { # which traits to color the tree by in auspice; titer colorbys are added in dengue_titers
				"country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
				"region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
				"gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"},
			},
			"controls": {'authors':['authors']},
			"defaults": {'geoResolution': 'region', 'colorBy': 'region', 'distanceMeasure': 'div', 'mapTriplicate': True}
			},

		"timetree_options": {"Tc": False},
		"fit_titer_model": not args.no_titers,
		"titers": { # regularization parameter values and cross-validation fraction
			"lam_avi":1.2,
			"lam_pot":0.5,
			"lam_drop":args.lam_drop,
			"training_fraction":0.9,
		},
		"estimate_mutation_frequencies": not args.no_mut_freqs,
		"estimate_tree_frequencies": not args.no_tree_freqs,
		"clean": args.clean,
		"pivot_spacing": 1.0/4, # pivots = time points; 1/N timepoints per year
 		}


##### Parse input files/params and run #####
if __name__=="__main__":
	parser = collect_args()
	args = parser.parse_args()

	### Find the right input files ###
	if args.json: # If provided, a specified JSON path overrides serotype argument
		args.json = [args.json]
	else: # Look for serotype-specific JSONs in the ./prepared/ directory
		if 'multiple' in args.serotypes: # "multiple" = run all 5 builds
			args.serotypes = ['denv1', 'denv2', 'denv3', 'denv4', 'all']
		else:
			args.serotypes = args.serotypes
		args.json = ['./prepared/dengue_%s.json'%s for s in args.serotypes] # Look for ./prepared/dengue_SEROTYPE.json if no file paths given

	for j in args.json:            # validate input JSONs exist
		assert os.path.isfile(j)


	### Run analyses ###
	for prepared_json in args.json:
		print("Processing %s"%prepared_json)
		runner = process(make_config(prepared_json, args)) # parse
		runner.align() # run alignment with mafft
		runner.build_tree() # build tree with fasttree -> raxml
		runner.timetree_setup_filter_run() # infer ML ancestral states (geo traits, node dates, mutations)
		runner.run_geo_inference() # run mugration model to infer transmissions

		# estimate mutation frequencies here.
		if runner.config["estimate_mutation_frequencies"]:
			pivots = runner.get_pivots_via_spacing()
			runner.estimate_mutation_frequencies(pivots=pivots, min_freq=0.02, inertia=np.exp(-1.0/12), stiffness=2)

		# estimate tree frequencies here.
		if runner.config["estimate_tree_frequencies"]: # methods @ [ref]
			pivots = runner.get_pivots_via_spacing()
			runner.estimate_tree_frequencies(pivots=pivots, stiffness=2) # stiffness ~= amount of smoothing

			for region in ['southeast_asia']:#, 'south_america']: #regions:
				try:
					runner.estimate_tree_frequencies(region=region, stiffness=2)
				except:
					continue
		# titers
		if runner.config["fit_titer_model"] and runner.config["titers"]: # methods @ Neher et al., PNAS 2016

			if args.titer_model == 'full_tree':
				titer_model_criterium = lambda node: True
			else:
				### Comparison: force dTiter values to be non-zero only on interserotype brances
				def is_interserotype(node):
					descendents = node.get_terminals()
					serotypes = [k.name.split('/')[0] for k in descendents if 'DENV' in k.name]
					serotypes = [s for s in serotypes if s != 'DENV']
					return len(set(serotypes)) > 1

				interserotype_branches = []
				for node in runner.tree.tree.find_clades():
					if is_interserotype(node):
						interserotype_branches.append(node)
						for child in node.clades:
							interserotype_branches.append(child)
				for node in runner.tree.tree.find_clades():
					if node in interserotype_branches:
						node.interserotype = True
					else:
						node.interserotype = False

				titer_model_criterium = lambda node: node.interserotype == True

			titer_model(runner, ## Run 10x with a 90:10 training:test split to estimate model performance / error
						lam_pot = runner.config['titers']['lam_pot'],
						lam_avi = runner.config['titers']['lam_avi'],
					lam_drop = runner.config['titers']['lam_drop'],
					training_fraction = runner.config['titers']['training_fraction'],
					sanofi_strain = sanofi_vaccine_strains[runner.info['lineage']], # vaccine strain for each serotype-specific build
					plot=False,
					criterium = titer_model_criterium, # calculate dTiter for all branches
					cross_validate=100,
					) # calculate dTiter for all branches

			titer_model(runner, ## Run once more with all the data to estimate branch effects for downstream analysis
						lam_pot = runner.config['titers']['lam_pot'],
						lam_avi = runner.config['titers']['lam_avi'],
					lam_drop = runner.config['titers']['lam_drop'],
					training_fraction = 1., # run again with all the data
					sanofi_strain = sanofi_vaccine_strains[runner.info['lineage']], # vaccine strain for each serotype-specific build
					plot=False,
					criterium = titer_model_criterium, # calculate dTiter for all branches
					) # calculate dTiter for all branches

			titer_export(runner)

	### Export for visualization in auspice
		# runner.auspice_export()
