from fuzzywuzzy import process
from Bio import SeqIO
import re
from collections import defaultdict
from pprint import pprint
import argparse
from glob import glob
import pandas as pd
import math

####    Parse arguments, pull input data    ####
parser = argparse.ArgumentParser()
parser.add_argument('-titers', default='/Users/Sidney/Dropbox/dengue/data/katzelnick2015/raw/agm_1month_titers.csv', type=str, help="Full path of titer table tsv")
parser.add_argument('-sequences', default='/Users/Sidney/Dropbox/dengue/data/katzelnick2015/metadata/Fig1A-aligned-nucleotide-sequences.FASTA', type=str, help="Full path of virus sequence file")
parser.add_argument('-key', default='/Users/Sidney/Dropbox/dengue/data/katzelnick2015/metadata/Fig1A-key-for-tree-names-virus-names.tsv', help="Full path to strain<\t>accession key")
parser.add_argument('-vdb', default='/Users/Sidney/nextstrain/fauna/data/dengue_all.fasta', type=str, help="Full path of vdb sequence file")
args = parser.parse_args()

titers_df = pd.read_csv(args.titers, index_col=0, comment='#', na_values='*')
titers_df.dropna(how='all', inplace=True)
viruses = list(titers_df.index.values)
sera = list(titers_df.columns.values)
fasta = [ s.description for s in SeqIO.parse(args.sequences, 'fasta') ]
key = [ line.split()[2] for line in open(args.key, 'r') if 'fullname' not in line ]
smith_synonyms = defaultdict(list, { s: [] for s in key })

vdb = [ s.description.split('|')[0] for s in SeqIO.parse(open(args.vdb, 'r'), 'fasta') ]
vdb_synonyms_file = glob('vdb_smith_synonyms.tsv')
if vdb_synonyms_file != []:
	vdb_synonyms_file = open(vdb_synonyms_file[0], 'r')
	vdb_synonyms = { line.split('\t')[0]:line.split('\t')[1:] for line in vdb_synonyms_file }
	vdb_synonyms = defaultdict(list, vdb_synonyms)
	vdb_synonyms_file.close()
else:
	vdb_synonyms = defaultdict(list)
#######		functions	########

def subset_strains(strains):
	'''  { (DENV1, 1990): [DENV1/COUNTRY1/STRAIN1/1990, ...] } '''
	subsets = defaultdict(list)
	for s in strains:
		try:
			sero = 'DENV' + re.search('DEN[V]?[1-4]{1}', s, flags=re.IGNORECASE).group(0)[-1] # DEN1 --> DENV1
			year = str(int(re.search('20[\d]{2}|19[\d]{2}', s).group(0)))
			subsets[(sero, year)].append(s)
		except:
			continue
	return subsets

def match_sero_year(strain, subset_data):
	''' return a subset of vdb viruses that have the correct year and serotype '''
	sero = 'DENV' + re.search('DEN[V]?[1-4]{1}', strain, flags=re.IGNORECASE).group(0)[-1] # DEN1 --> DENV1
	year = str(int(re.search('20[\d]{2}|19[\d]{2}', strain).group(0)))
	choices = subset_data[(sero, year)]
	return choices

def get_strain_id(strain):
	''' hacky way of parsing where smith datasets put the strain ID itself '''
	if '/' in strain:
		if strain.count('/') == 3:
			if re.search('20[\d]{2}|19[\d]{2}', strain.split('/')[2]): #'DENV2/Vietnam/2006/32125'
				return strain.split('/')[3]
			else:
				return strain.split('/')[2] #'DENV2/Vietnam/32125/2006'
		else:
			yearstrain = strain.split('/')[-1] #'DENV2/Vietnam/2006-32-135'
			return ''.join(yearstrain.split('-')[1:])
	else: # 'DEN2_Vietnam_2003_DF_670'
		return ''.join(strain.split('_')[3:])

def match_broad_ids(strain1, strain2):
	''' if both strains have broad IDs, they should match exactly; return boolean'''
	BIDV1 = re.search('BID[-_]?[V]?[\d]{1,6}', strain1, flags=re.IGNORECASE)
	BIDV2 = re.search('BID[-_]?[V]?[\d]{1,6}', strain2, flags=re.IGNORECASE)
	if BIDV1 and not BIDV2: # only one has a BID-V ID
		return False
	elif BIDV2 and not BIDV1:
		return False
	elif not BIDV1 and not BIDV2: # neither have BID-V IDs
		return None
	else: # do they match exactly?
		BIDV1 = BIDV1.group(0).upper().replace('-', '').replace('_', '').strip()
		BIDV2 = BIDV2.group(0).upper().replace('-', '').replace('_', '').strip()
		return BIDV1 == BIDV2

def get_closest_match(string, choices, cutoff=87, p=False):
	''' return fuzzywuzzy matches for strain names if they pass score cutoff '''
	closestmatch = process.extractOne(string, choices)
	if closestmatch and closestmatch[1] >= cutoff:
		if p:
			print closestmatch
		return closestmatch[0]
	else:
		return None

def match_strain_names(s, subset_data, p=False):
	''' try and match strain names to existing records in vdb;
	return dictionary of matches and list of orphans'''
	choices = match_sero_year(s, subset_data)   # only compare to names with the right serotype and year
	match = get_closest_match(s, choices, cutoff=87, p=p)	# match whole names first
	if not match:								# try and match just the strain identifier, but require higher specificity
		match = get_closest_match(get_strain_id(s), choices, cutoff=90, p=p)
	if match and match_broad_ids(s, match) != False: # require broad IDs to match exactly (if present)
		return match
	else:
		return None


def make_vdb_name(strain):
	''' return standardized strain name like DENV1/COUNTRY/STRAINID/YEAR '''
	strain = strain.upper().replace('PUERTO_RICO', 'PUERTORICO').replace('PUERTO-RICO', 'PUERTORICO')
	sero = 'DENV' + re.search('DEN[V]?[1-4]{1}', strain, flags=re.IGNORECASE).group(0)[-1] # DEN1 --> DENV1
	year = str(int(re.search('20[\d]{2}|19[\d]{2}', strain).group(0)))
	country = re.search('[a-z,A-Z]{3,20}[_\/]{1}|[_\/]{1}[a-z,A-Z]{3,20}|[\/]{1}[a-zA-z]{4,6}_[a-zA-z]{4,6}[\/]{1}', strain, flags=re.IGNORECASE).group(0)
	country = country.replace('_', '').replace('/', '').replace('PUERTORICO', 'PUERTO_RICO')
	for p in ['DEN[V]?[1-4]{1}', year, country.replace('/', '').replace('_', ''), country, '/', '_', '-', ' ', 'sylvatic']:
		strain = re.sub(p, '', strain, flags=re.IGNORECASE)
	return '/'.join([sero, country, strain, year]).upper()


########	Group smith strain names together		###########
if 'monovalent' in args.titers:
	all_to_match = fasta+viruses
else:
	all_to_match = fasta+viruses+sera
smith_synonyms_subset = subset_strains(smith_synonyms)
fasta_subset = subset_strains(fasta)
viruses_subset = subset_strains(viruses)
sera_subset = subset_strains(sera)

dataset_priorities = [smith_synonyms_subset, fasta_subset, viruses_subset, sera_subset]

for s in all_to_match:
	if s in smith_synonyms:
		continue
	found_match = None
	priority = 0
	while not found_match and priority < len(dataset_priorities):
		found_match = match_strain_names(s, dataset_priorities[priority])
		priority += 1
	else:
		if found_match == None:
			smith_synonyms[s] = []
			smith_synonyms_subset = subset_strains(smith_synonyms)
		else:
			smith_synonyms[found_match].append(s)
smith_synonyms = dict(smith_synonyms)

######	Match to VDB names		#######
vdb_synonyms_subset = subset_strains(vdb_synonyms.keys()) # Previously identified vdb synonyms
vdb_subset = subset_strains(vdb)	# Other strains in vdb

dataset_priorities = [ vdb_synonyms_subset, vdb_subset ]
for s in smith_synonyms.keys():
	found_match = None
	priority = 0
	while not found_match and priority < len(dataset_priorities):
		found_match = match_strain_names(s, dataset_priorities[priority])
		priority += 1
	else:
		if found_match:
			vdb_synonyms[found_match].append(s)
			vdb_synonyms[found_match] += smith_synonyms[s]
			vdb_synonyms_subset = subset_strains(vdb_synonyms)
		else:
			newname = make_vdb_name(s)
			vdb_synonyms[newname+'*'] = [s]
			vdb_synonyms[newname+'*'] += smith_synonyms[s]

#####	Write to file	####
vdb_synonyms_file = open('vdb_smith_synonyms.tsv', 'w')

for v, synonyms in sorted(vdb_synonyms.items()):
	# print v+'\t'+'\t'.join(list(set(synonyms)))+'\n'
	vdb_synonyms_file.write(v+'\t'+'\t'.join(list(set(synonyms)))+'\n')

vdb_synonyms_file.close()
