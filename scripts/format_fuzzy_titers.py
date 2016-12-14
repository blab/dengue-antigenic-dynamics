from fuzzywuzzy import process
from Bio import SeqIO
import re
from collections import defaultdict
from pprint import pprint

####     Pull data       ####
smith_strains = sorted([ line.split()[0] for line in open('/Users/Sidney/Dropbox/dengue/data/smith2015/key1A_strain_synonymns.tsv', 'r') if not line.startswith('fullname') ])
smith_sequences = { s.description:str(s.seq) for s in SeqIO.parse('/Users/Sidney/Dropbox/dengue/data/smith2015/Fig1A-aligned-nucleotide-sequences.FASTA', 'fasta')}

vdb = [ s.description for s in SeqIO.parse(open('/Users/Sidney/nextstrain/fauna/data/dengue.fasta', 'r'), 'fasta') ]
vdb_strains = [s.split('|')[0] for s in vdb]
vdb_accessions = {s.split('|')[1] : s.split('|')[0] for s in vdb }

####    Define functions        ####
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
    BIDV1 = re.search('BID[-]?[V]?[\d]{1,6}', strain1, flags=re.IGNORECASE)
    BIDV2 = re.search('BID[-]?[V]?[\d]{1,6}', strain2, flags=re.IGNORECASE)
    if BIDV1 and not BIDV2: # only one has a BID-V ID
        return False
    elif BIDV2 and not BIDV1:
        return False
    elif not BIDV1 and not BIDV2: # neither have BID-V IDs
        return None
    else: # do they match exactly?
        return BIDV1.group(0).upper().replace('-', '').replace('_', '').strip() == BIDV2.group(0).upper().replace('-', '').replace('_', '').strip()

def get_closest_match(string, choices, cutoff=87):
    ''' return fuzzywuzzy matches for strain names if they pass score cutoff '''
    closestmatch = process.extractOne(string, choices)
    if closestmatch and closestmatch[1] >= cutoff:
        return closestmatch[0]
    else:
        return None

def match_strain_names(titerstrains, subset_data):
    ''' try and match strain names to existing records in vdb;
    return dictionary of matches and list of orphans'''
    rename = {}
    unmatched = []
    for t in sorted(list(set(titerstrains))):   # first try and match whole name
        choices = match_sero_year(t, subset_data)   # only compare to names with the right serotype and year
        match = get_closest_match(t, choices, 87)
        if match and match_broad_ids(t, match) != False: # require broad IDs to match exactly (if present)
            rename[t] = match
            continue
        else:
            unmatched.append(t)

    for t in unmatched: # try and match just the strain identifier, but require higher specificity
        choices = match_sero_year(t, subset_data)
        match = get_closest_match(get_strain_id(t), choices, 90)
        if match and match_broad_ids(t, match) != False: # still enforce identical broad IDs
            rename[t] = match
            continue

    unmatched = set(unmatched).difference(set(rename.keys()))
    print 'Found matches in vdb for these strains:'
    pprint(sorted(rename.items()), width=130)
    return rename, unmatched

def find_sequences(missing_strains):
    ''' find the closest strain name in the fasta file, return corresponding sequence '''
    unmatched = []
    for s in missing_strains:
        choices = match_sero_year(s, smith_sequences_subsets) # check against sequences with the right serotype and year
        match = get_closest_match(get_strain_id(s), choices, 80, p=True) # match on strain identifier
        if match and match_broad_ids(t, match) != False: # require broad IDs to be identical if they exist
            return smith_sequences[match]
        else:
            return None

####    Find strains with existing VDB records      ####

vdb_strain_subsets = subset_strains(vdb_strains)
found, missing = match_strain_names(smith_strains, vdb_strain_subsets)

####    Find sequence records for strains not in VDB    ####
smith_sequences_subsets = subset_strains(smith_sequences.keys())
