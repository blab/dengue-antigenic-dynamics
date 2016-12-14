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
parser.add_argument('-table', default='/Users/Sidney/Dropbox/dengue/data/smith2015/agm_1month_titers.csv', type=str, help="Full path of titer table tsv")
parser.add_argument('-sequences', default='/Users/Sidney/Dropbox/dengue/data/smith2015/Fig1A-aligned-nucleotide-sequences.FASTA', type=str, help="Full path of virus sequence file")
parser.add_argument('-key', default='/Users/Sidney/Dropbox/dengue/data/smith2015/Fig1A-key-for-tree-names-virus-names.tsv', help="Full path to strain<\t>accession key")
parser.add_argument('-vdb', default='/Users/Sidney/nextstrain/fauna/data/dengue.fasta', type=str, help="Full path of vdb sequence file")
parser.add_argument('-src', default = 'agm_1mo', type=str,  help="source designation")
args = parser.parse_args()

titers_df = pd.read_csv(args.table, index_col=0, comment='#', na_values='*')
titers_df.dropna(how='all', inplace=True)

strain_acc = { line.split()[2]:line.split()[3] for line in open(args.key, 'r')}
smith_sequences = { s.description:str(s.seq) for s in SeqIO.parse(args.sequences, 'fasta')}

vdb = [ s.description for s in SeqIO.parse(open(args.vdb, 'r'), 'fasta') ]
vdb_strains = [s.split('|')[0] for s in vdb]
# vdb_accessions = {s.split('|')[1] : s.split('|')[0] for s in vdb }

####    Pull existing titer measurements (if any)   #####

existing_titers = glob('dengue_strains.tsv')
if existing_titers != []:
    existing_titers = { line.split('\t')[0]:int(line.split('\t')[1]) for line in open(existing_titers[0], 'r') }

try:
    seen_smith_viruses = glob('smith_viruses.tsv')
    seen_smith_viruses = pd.read_csv(seen_smith_viruses[0], dtype='str', delimiter='\t', index_col=0, header=0, skiprows=1)
    seen_smith_viruses = list(seen_smith_viruses['Name'])
except:
    seen_smith_viruses = []

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
    # print 'Found matches in vdb for these strains:'
    # pprint(sorted(rename.items()), width=130)
    return rename, unmatched

def find_sequence(missing_strain):
    ''' find the closest strain name in the fasta file, return corresponding sequence '''
    choices = match_sero_year(missing_strain, smith_sequences_subsets) # check against sequences with the right serotype and year
    match = get_closest_match(get_strain_id(missing_strain), choices, 80) # match on strain identifier
    if match and match_broad_ids(missing_strain, match) != False: # require broad IDs to be identical if they exist
        return match
    else:
        return None

def make_vdb_name(strain):
    ''' return standardized strain name like DENV1/COUNTRY/STRAINID/YEAR '''
    sero = 'DENV' + re.search('DEN[V]?[1-4]{1}', strain, flags=re.IGNORECASE).group(0)[-1] # DEN1 --> DENV1
    year = str(int(re.search('20[\d]{2}|19[\d]{2}', strain).group(0)))
    country = re.search('[a-z,A-Z]{3,20}[_\/]{1}|[_\/]{1}[a-z,A-Z]{3,20}|[\/]{1}[a-zA-z]{4,6}_[a-zA-z]{4,6}[\/]{1}', strain, flags=re.IGNORECASE).group(0)
    for p in [sero, year, country, '/', '_', '-']:
        strain = re.sub(p, '', strain, flags=re.IGNORECASE)
    country = country.replace('_', '').replace('/', '')
    return '/'.join([sero, country, strain, year]).upper()

def format_smith_upload(missing_strains, no_append=True):
    '''
    Given strain list, pull metadata from strain name and put into LANL format
    to make uploading missing records to vdb easy.
    Write this to tsv and print suggested upload command.
    Return list of any strains that weren't found in the fasta and/or key.
    '''
    records = []
    not_in_key = []
    for s in missing_strains:
        newstrain = make_vdb_name(s)
        if newstrain in seen_smith_viruses:
            continue
        try:
            fasta_match = find_sequence(s)
            accession = strain_acc[fasta_match]
        except:
            not_in_key.append(s)
            continue
        row = {}

        NA_fields = ['Species', 'Isolate Name', 'Georegion', 'Author', 'Sampling City', 'Pubmed ID']
        for n in NA_fields:
            row[n] = None
        row['Serotype'], row['Country'], row['Strain'], row['Sampling Year'] = newstrain.split('/')
        row['Accession'] = accession
        row['Name'] = newstrain
        row['Start'] = 935
        row['Stop'] = 2413
        row['Segment'] = 'E'
        row['Organism'] = 'dengue virus '+row['Serotype'][-1]
        row['Species'] = 'dengue virus'
        row['Sequence'] = smith_sequences[fasta_match]
        records.append(row)

    if no_append:
        open('smith_viruses.tsv', 'a').write('# viruses from smith data not in vdb\n')
    else:
        pass
    pd.DataFrame(records).to_csv('smith_viruses.tsv', sep='\t', mode='a', header=no_append)
    return not_in_key

def table_to_tsv(titerdf):
    titer_file = open('dengue_titers.tsv', 'a')
    strain_file = open('dengue_strains.tsv', 'w')
    value_counts = defaultdict(int, existing_titers)
    for virus, seraseries in titerdf.iterrows():
        for sera, value in seraseries.iteritems():
            if type(value)== float and math.isnan(value):
                continue
            elif value == '<10':
                continue
            else:
                titer_file.write('\t'.join([virus, sera, sera, args.src, str(value)])+'\n')
                value_counts[virus] += 1

    for virus in titerdf.index.values:
        strain_file.write('\t'.join([virus, str(value_counts[virus])])+'\n')

    titer_file.close()
    strain_file.close()


####    Find strains with existing VDB records      ####
vdb_strain_subsets = subset_strains(vdb_strains)
smith_sequences_subsets = subset_strains(smith_sequences.keys())

found_viruses, missing_viruses = match_strain_names(titers_df.columns.values, vdb_strain_subsets)
for v in missing_viruses:
    found_viruses[v] = make_vdb_name(v)
found_sera, missing_sera = match_strain_names(titers_df.index.values, vdb_strain_subsets)
for s in missing_sera:
    found_sera[s] = make_vdb_name(s)
titers_df.rename(index = found_viruses, columns = found_sera, inplace=True)
table_to_tsv(titers_df)

####    Deal with records not in VDB    ####
missing_strains = list(set(list(missing_viruses) + list(missing_sera)))

print 'These strains were not found in vdb, but were in provided fasta file.\nTo upload, run `fauna$ python vdb/dengue_upload --fname smith_viruses.tsv --ftype tsv -v dengue -db vdb`\n'
not_in_key = format_smith_upload(missing_strains)
pprint(set(missing_strains).difference(set(not_in_key)))
print 'WARNING: these strains were not found in vdb or provided fasta file:'
pprint(not_in_key)
