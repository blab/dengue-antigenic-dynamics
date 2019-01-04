from __future__ import print_function
import os, sys
# we assume (and assert) that this script is running from the virus directory, i.e. inside H7N9 or zika
sys.path.append(os.path.join(os.path.dirname(__file__), '../'))
import base.process
from base.process import process
import argparse
import numpy as np
from dengue_titers import titer_model, titer_export ## Set up and parameterize the titer model separately for tidiness

##### Parse args, set up config #####
def collect_args():
	parser = base.process.collect_args()
	# parser.add_argument('-j', '--json', default=None, nargs=1, type=str, help="Accepts path to prepared JSON; overrides -s argument")
	parser.add_argument('-s', '--serotypes', default="all", nargs=1, type=str, choices=['denv1', 'denv2', 'denv3', 'denv4', 'all', 'multiple'],
	help="Look for prepared JSON(s) like ./prepared/dengue_SEROTYPE.json; 'multiple' will run all five builds. Default='multiple'")
	parser.add_argument('--no_tree_freqs', default=False, action='store_true', help="skip tree (clade) frequencies")
	parser.add_argument('--no_titers', default=False, action='store_true', help="skip titer models")
	parser.add_argument('--titer_model', default='substitution', choices=['tree', 'substitution'], type=str)
	parser.add_argument('--titer_res', '--titer_resolution', default='full_tree', choices=['full_tree', 'interserotype'], type=str)
	return parser

parser = collect_args()
args = parser.parse_args()

### Find the right input files ###
if args.json is None:  # Look for ./prepared/dengue_SEROTYPE.json if no file path given
	args.json = './prepared/dengue_%s.json'%args.serotypes
assert os.path.isfile(args.json)

auspice_output_path = "../../%s_%s_output/"%(args.titer_model, args.titer_res)
if not os.path.isdir(auspice_output_path):
	os.mkdir(auspice_output_path)

### Specify analysis config ###
config = {
"dir": "dengue",
"in": args.json,
"output": {"auspice": auspice_output_path,
			"data": "./processed/"},
"timetree_options": {"Tc": False},
"estimate_tree_frequencies": not args.no_tree_freqs,

"fit_titer_model": not args.no_titers,
"titers": { # regularization parameter values and cross-validation fraction
	"lam_avi":0.0,
	"lam_pot":0.5,
	"lam_drop":1.0,
	"training_fraction":0.9,
	},

"clean": args.clean,
"pivot_spacing": 1.0/4, # pivots = time points; 1/N timepoints per year
}
# "geo_inference": ['region'], # what traits to perform this on; don't run country (too many demes, too few sequences per deme to be reliable)
# "auspice": { ## settings for auspice JSON export
# 	"extra_attr": ['serum', 'clade', 'dTiter', 'cTiter'], # keys from tree.tree.clade['attr'] to include in export
# 	"color_options": { # which traits to color the tree by in auspice; titer colorbys are added in dengue_titers
# 		"country":{"key":"country", "legendTitle":"Country", "menuItem":"country", "type":"discrete"},
# 		"region":{"key":"region", "legendTitle":"Region", "menuItem":"region", "type":"discrete"},
# 		"gt":{"key":"genotype", "legendTitle":"Genotype", "menuItem":"genotype", "type":"discrete"},
# 	},
# 	"controls": {'authors':['authors']},
# 	"defaults": {'geoResolution': 'region', 'colorBy': 'region', 'distanceMeasure': 'div', 'mapTriplicate': True}
# 	},


##### Parse input files/params and run #####


### Run analyses ###
print("Processing %s"%args.json)
runner = process(config) # parse
runner.align() # run alignment with mafft
runner.build_tree() # build tree with fasttree -> raxml
runner.timetree_setup_filter_run() # infer ML ancestral states (geo traits, node dates, mutations)
# runner.run_geo_inference() # run mugration model to infer transmissions

# estimate tree frequencies here.
if runner.config["estimate_tree_frequencies"]:
	pivots = runner.get_pivots_via_spacing()
	runner.estimate_tree_frequencies(region='southeast_asia', stiffness=2)

# titers
if runner.config["fit_titer_model"] and runner.config["titers"]: # methods @ Neher et al., PNAS 2016
	if args.titer_model=='tree':
		n=1000

		if args.titer_res == 'full_tree':
			titer_criterium = lambda node: True

		elif args.titer_res == 'interserotype':
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
				node.interserotype = node in interserotype_branches

			titer_criterium = lambda node: node.interserotype == True

	elif args.titer_model == 'substitution':
		n = 100
		titer_criterium = None

		if args.titer_res == 'interserotype':
			serotypes = ['DENV1', 'DENV2', 'DENV3', 'DENV4']
			serotype_tips = {sero: [] for sero in serotypes}
			for k in runner.tree.tree.get_terminals(): # find all tips with a strain name corresponding to each serotype
				sero = k.name.split('/')[0]
				if sero in serotypes:
					serotype_tips[sero].append(k)
			serotype_mrcas = {sero: runner.tree.tree.common_ancestor(serotype_tips[sero]) for sero in serotypes} # find the serotype mrca

			for sero, mrca in serotype_mrcas.items():
				for k in mrca.find_clades(): # pull all descendants (nodes and tips) of the serotype mrca
					k.sequence = mrca.sequence # set sequence of each descendant to the reconstructed ancestral sequence from the serotype mrca
					k.translations = mrca.translations
					k.aa_mutations = {}
			runner.tree.add_translations() # translate and reassign mutations to each branch
			runner.tree.refine()

			for sero, mrca in serotype_mrcas.items():
				for k in mrca.find_clades(): # pull all descendants (nodes and tips) of the serotype mrca
					if k == mrca: continue
					assert ''.join(k.sequence) == ''.join(mrca.sequence)
					assert k.translations['E'] == mrca.translations['E']
					assert len(k.aa_mutations['E']) == 0

	titer_model(runner, ## Run N times with a random 90:10 training:test split to estimate model performance / error
				lam_pot = runner.config['titers']['lam_pot'],
				lam_avi = runner.config['titers']['lam_avi'],
			lam_drop = runner.config['titers']['lam_drop'],
			model_type=args.titer_model,
			training_fraction = runner.config['titers']['training_fraction'],
			plot=False,
			criterium = titer_criterium, # calculate dTiter for all branches
			cross_validate=n,
			) # calculate dTiter for all branches

	titer_model(runner, ## Run once more with all the data to estimate branch effects for downstream analysis
				lam_pot = runner.config['titers']['lam_pot'],
				lam_avi = runner.config['titers']['lam_avi'],
			lam_drop = runner.config['titers']['lam_drop'],
			model_type=args.titer_model,
			training_fraction = 1., # run again with all the data
			plot=False,
			criterium = titer_criterium, # calculate dTiter for all branches
			) # calculate dTiter for all branches

	titer_export(runner, args.titer_model)

### Export for visualization in auspice
runner.auspice_export()
