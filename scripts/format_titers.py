import pandas as pd
import numpy as np
from Bio import SeqIO
from collections import defaultdict
import argparse
import math
from pprint import pprint
import sys
from glob import glob

'''
Purpose: Match Smith dataset strain names to vdb strain names.
Smith strains --> accessions --> vdb records (if already in db; otherwise write to tsv and upload to vdb)
'''

parser = argparse.ArgumentParser()
parser.add_argument('-table', default='/Users/Sidney/Dropbox/dengue/data/smith2015/agm_1month_titers.csv', type=str, help="Full path of titer table tsv")
parser.add_argument('-sequences', default='/Users/Sidney/Dropbox/dengue/data/smith2015/Fig1A-aligned-nucleotide-sequences.FASTA', type=str, help="Full path of virus sequence file")
parser.add_argument('-key', default='/Users/Sidney/Dropbox/dengue/data/smith2015/Fig1A-key-for-tree-names-virus-names.txt', help="Full path to strain<\t>accession key")
parser.add_argument('-vdb', default='/Users/Sidney/nextstrain/fauna/data/dengue.fasta', type=str, help="Full path of vdb sequence file")
parser.add_argument('-outf', default = 'dengue', type=str,  help="file stem for outfiles")
parser.add_argument('-src', default = 'agm_1mo', type=str,  help="source designation")
args = parser.parse_args()

seen_titers = glob('%s_strains.tsv'%args.outf)
if seen_titers != []:
    existing_titers = { line.split()[0]:int(line.split()[1]) for line in open(seen_titers[0], 'r') }
else:
    existing_titers = {}

seen_not_in_vdb = glob('smith_viruses.tsv')
if seen_not_in_vdb != []:
    existing_not_in_vdb = pd.read_csv(seen_not_in_vdb[0], dtype='str', delimiter='\t', index_col=0, header=0, skiprows=1)
    existing_not_in_vdb = list(existing_not_in_vdb['Accession'])
else:
    existing_not_in_vdb = []
############    Parse input data      ################

titers_df = pd.read_csv(args.table, index_col=0, comment='#', na_values='*')
titers_df.dropna(how='all', inplace=True)
smith_key_df = pd.read_csv(args.key, header=0, sep='\t',index_col=0)

# { vdb_strain: acc }
vdb_acc_strains = { s.description.split('|')[1] : s.description.split('|')[0].split('.')[0] for s in SeqIO.parse(args.vdb, 'fasta')}
vdb_acc_strains['AF326573']='DENV4/DOMINICA/81669/1981'
vdb_acc_strains['AF038403']='DENV2/PAPUANEWGUINEA/NEWGUINEAC/1944'
vdb_acc_strains['L11422']='DENV3/FIJI/29472/1992'

# { acc: vdb_strain } (keep track of already-fixed names)
found_acc_vdb = {'AF326573': 'DENV4/DOMINICA/81669/1981', 'AF038403': 'DENV2/PAPUANEWGUINEA/NEWGUINEAC/1944', 'L11422':'DENV3/FIJI/29472/1992'}

# { smith_strain: accession }, manually add missing metadata
smith_strain_acc = { s['fullname'] : s['genbank'] for i,s in smith_key_df.iterrows()}
smith_strain_acc['DEN1_Cambodia_2003_GenBankGQ868619'] = 'GQ868619'
smith_strain_acc['DEN1_Burma_2005_61117'] = 'KT452791'
smith_strain_acc['DEN2_Vietnam_2003_DF_670']='KT452797'
smith_strain_acc['DEN3_Burma_2008_80931'] = 'KT452792'
smith_strain_acc['DEN1_Bolivia_2010_FSB_3363'] = 'KT382187'
smith_strain_acc['DEN4_Burma_2008_81087'] = 'KT452793'
smith_strain_acc['DENV2/VIETNAM/AC212/2003'] = 'KT452796'
smith_strain_acc['DENV3/Fiji/1992-29472'] = 'L11422'
smith_strain_acc['DEN3_Fiji_1992__']='L11422'
smith_strain_acc['DEN2_Senegal_2003_Sendak_HD_0674']='EF105384'
smith_strain_acc['DENV2/Tonga/1974-Tonga-74']='AY744147'
smith_strain_acc['DENV3/Fiji/1992-29472'] = 'L11422'
smith_strain_acc['DENV2/Senegal/1970/Sendak_H D_0674'] = smith_strain_acc['DENV2/Senegal/1970/Sendak_H']
smith_strain_acc['DENV2/SENEGAL/SENDAKHD0674/1970'] = 'EF105384'
smith_strain_acc['DENV2/Senegal/1970/Sendak_HD_0674'] = 'EF105384'
smith_strain_acc['DENV2/SENEGAL/1970/SENDAKHD0674'] = 'EF105384'

# { accession: sequence object }
smith_sequences = { smith_strain_acc[s.description] : s for s in SeqIO.parse(args.sequences, 'fasta')}


###############     Convert virus and sera strain names to existing or created vdb strain names ###############

def pull_virus_smithmetadata(strain):
    serotype, country, yearstrain = strain.split('/', 2)
    yearstrain = yearstrain.replace('/', '-')
    year = yearstrain.split('-')[0]
    strain = ''.join(yearstrain.split('-')[1:]).split(' ')[0]
    return serotype, country, strain, year

def pull_sera_smithmetadata(strain):
    serotype, country, year, strain = strain.split('_', 3)
    serotype = 'DENV'+serotype[-1]
    return serotype, country, strain, year

def fix_smith_strain(strain, type='virus'):
    '''
    Make new strain names like
    DENV1234/country/ID/year
    '''

    if type == 'virus':
        sero, country, strain_id, year = pull_virus_smithmetadata(strain) # Pull metadata from pre-processed annotations
    else:
        sero, country, strain_id, year = pull_sera_smithmetadata(strain)
    strain = '%s/%s/%s/%s'%(sero, country, strain_id, year)
    strain = strain.replace('-', '').replace('_', '').upper().strip()
    return strain

# Difficult to match from Smith virus strains -> Smith sera strains -> acc.
# So, make everything into a canonical name first.

for s in smith_strain_acc.keys():
    try:
        smith_strain_acc[fix_smith_strain(s)] = smith_strain_acc[s]
    except:
        smith_strain_acc[fix_smith_strain(s, type='sera')] = smith_strain_acc[s]

def convert_smith_vaccine_strains(titerdf):
    vaccine_strains = {'DEN1': 'DENV1/NA/WESTERNPACIFICDELTA30/NA', 'DEN2': 'DENV2/NA/NEWGUINEAC/NA', 'DEN3': 'DENV3/NA/SLEMANDELTA30/NA', 'DEN4': 'DENV4/NA/81669DELTA30/NA'}
    titer_file = open(args.outf+'_titers.tsv', 'a')
    strain_file = open(args.outf+'_strains.tsv', 'w')
    value_counts = defaultdict(int, existing_titers)
    for virus, seraseries in titerdf.iterrows():
        for sera, value in seraseries.iteritems():
            smith_fields = sera.split('_')
            sera = vaccine_strains[smith_fields[0]]
            ID = smith_fields[-1]
            if type(value)== float and math.isnan(value):
                continue
            elif value == '<10':
                continue
            else:
                titer_file.write('\t'.join([virus, sera, ID, args.src, str(value)])+'\n')
                value_counts[virus] += 1

    for virus in titerdf.index.values:
        strain_file.write('\t'.join([virus, str(value_counts[virus])])+'\n')
    titer_file.close()
    strain_file.close()

def convert_smith_vdb_strains(smith_strains, type='virus'):
    smith_vdb_strains = {}
    not_in_key = []
    not_in_vdb = []

    def try_acc(attempt):
        if attempt in smith_strain_acc: # Do we have an accession for this strain name?
            return smith_strain_acc[attempt]
        else:
            return None

    def try_vdb(acc, smith, fixed): # Is the accession in vdb?
        if acc in vdb_acc_strains:
            smith_vdb_strains[smith] = vdb_acc_strains[acc]
        else:
            not_in_vdb.append(smith)
            smith_vdb_strains[smith] = fixed

    for smith_strain in smith_strains:
        fixed_smith_strain = fix_smith_strain(smith_strain, type=type) # Try both the original and canonicalized name

        for name in [smith_strain, fixed_smith_strain]:
            acc = try_acc(name) # Can we find an accession for this strain?
            if acc:
                if acc in found_acc_vdb:
                    smith_vdb_strains[smith_strain] = found_acc_vdb[acc]
                else:
                    try_vdb(acc, smith_strain, fixed_smith_strain) # If we have an accession, is it in the vdb?
                break
            else:
                continue

        if smith_strain not in smith_vdb_strains: # didn't find it in the key? keep track and raise warning.
            not_in_key.append((smith_strain, fixed_smith_strain))

    if not_in_key != [] and type == 'virus':
        print 'These %s strain names not found in key:\n'%type, sorted(not_in_key)
    return smith_vdb_strains, not_in_key, not_in_vdb

def match_virus_sera(smith_vdb_virus_strains, smith_vdb_sera_strains, sera_not_in_key):
    serocountryyear_strain = defaultdict(list) # { concatenated serotypeCountryYear: [list of all strains with this combo]}

    for s in smith_vdb_virus_strains.values():
        fields = s.split('/')
        serocountryyear = fields[0]+fields[1]+fields[3]
        serocountryyear = serocountryyear.strip()
        serocountryyear_strain[serocountryyear].append(s)
    serocountryyear_strain = dict(serocountryyear_strain) # disable default dict behavior
    for k,v in serocountryyear_strain.items():
        if len(v) == 1:
            serocountryyear_strain[k] = v[0] # single strain meets criteria? great, no chance of mixups here.
        else:
            del serocountryyear_strain[k] # otherwise, don't use this as a proxy.

    for s in sera_not_in_key:
        s, fixed_s = s
        fields = fixed_s.split('/')
        serocountryyear = fields[0]+fields[1]+fields[3]
        serocountryyear = serocountryyear.strip()
        try:
            smith_vdb_sera_strains[s] = serocountryyear_strain[serocountryyear] # try to match to previous
        except:
            continue

    sera_not_in_key = [ s for s in sera_not_in_key if s[0] not in smith_vdb_sera_strains ] # remove the ones we found matches to
    if sera_not_in_key != []:
        pprint('still no match for these sera strains:' )
        pprint(sorted([s[0] for s in sera_not_in_key])) # warn and output still unmatched strains
        print '\n\nserocountryyear keys:\n', sorted(serocountryyear_strain.keys())
    return smith_vdb_sera_strains, sera_not_in_key

############### Write missing viruses and metadata to file ###############

def format_smith_upload(acc_list):
    '''
    Given acc list, pull metadata from strain name and put into LANL format
    to make uploading missing records to vdb easy.
    Write this to tsv and print suggested upload command.
    '''
    if existing_not_in_vdb == []:
        open('smith_viruses.tsv', 'a').write('# %d viruses from smith data not in vdb\n'%len(acc_list))
        header = True
    else:
    	header = False
    missing_records = []
    for a in acc_list:
        if a in existing_not_in_vdb:
            continue
        row = {}
        try:
            s = found_acc_vdb[a]
        except:
            s = fix_smith_strain(smith_sequences[a].description)

        NA_fields = ['Species', 'Isolate Name', 'Georegion', 'Author', 'Sampling City', 'Pubmed ID']
        for n in NA_fields:
            row[n] = None
            serotype, country, strain, year = s.split('/')
        row['Accession'] = smith_strain_acc[s]
        row['Name'] = s
        row['Start'] = 935
        row['Stop'] = 2413
        row['Segment'] = 'E'
        row['Organism'] = 'dengue virus '+serotype[-1]
        row['Country'] = country
        row['Sampling Year'] = year
        row['Species'] = 'dengue virus'
        row['Sequence'] = str(smith_sequences[a].seq)
        missing_records.append(row)
    pd.DataFrame(missing_records).to_csv('smith_viruses.tsv', sep='\t', mode='a', header=header)

def table_to_tsv(titerdf, ofile_stem):
    titer_file = open(ofile_stem+'_titers.tsv', 'a')
    strain_file = open(ofile_stem+'_strains.tsv', 'w')
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

########## Run ###############

# Parse virus strain names
smith_vdb_virus_strains, missing_key_viruses, missing_vdb_viruses = convert_smith_vdb_strains(titers_df.index.values)
acc_vdb_strains = { smith_strain_acc[smith_strain]:smith_vdb_virus_strains[smith_strain]
                    for smith_strain in smith_vdb_virus_strains
                          if smith_strain in smith_strain_acc }

if 'monovalent' in args.src:
    # Name known vaccine strains by serotype; handle writing separately (these come with ferret_ids, unlike all the other data)
    titers_df.rename(index = smith_vdb_virus_strains, inplace=True)
    convert_smith_vaccine_strains(titers_df)
    # Accessions not in vdb (found in viruses)
    not_in_vdb = set([smith_strain_acc[fix_smith_strain(v)] for v in missing_vdb_viruses])
else:
    # Parse sera strain names
    smith_vdb_sera_strains, missing_key_sera, missing_vdb_sera = convert_smith_vdb_strains(titers_df.columns.values, type='sera')
    # Try and match uncaptured sera names to parsed virus names
    smith_vdb_sera_strains, missing_key_sera = match_virus_sera(smith_vdb_virus_strains, smith_vdb_sera_strains, missing_key_sera)
    # Rename dataset with standardardized strain names (that match vdb)
    titers_df.rename(index = smith_vdb_virus_strains, columns = smith_vdb_sera_strains, inplace=True)
    table_to_tsv(titers_df, args.outf)
    # Accessions not in vdb (found in either viruses or sera)
    not_in_vdb = set([smith_strain_acc[fix_smith_strain(v)] for v in missing_vdb_viruses] + [smith_strain_acc[fix_smith_strain(s, type='sera')] for s in missing_vdb_sera])
    table_to_tsv(titers_df, args.outf)

print 'These accessions were not found in vdb.\nTo upload, run `fauna$ python vdb/dengue_upload --fname smith_viruses.tsv --ftype tsv -v dengue -db vdb`\n', not_in_vdb
format_smith_upload(not_in_vdb)

print '\n\nWrote all {strain name synonymns : accessions} to tsv: key1A_strain_synonymns.tsv'
key_synonymns = open('key1A_strain_synonymns.tsv', 'a')
key_synonymns.write('fullname\tgenbank\n')
for s,a in sorted(smith_strain_acc.items(), key=lambda i: i[1]):
    key_synonymns.write('%s\t%s\n'%(s,a))
key_synonymns.close()
