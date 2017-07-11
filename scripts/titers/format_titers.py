from Bio import SeqIO
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
parser.add_argument('-synonyms', default='./vdb_smith_synonyms.tsv', type=str, help="Full path of vdbName<\t>all<\t>smith<\t>synonyms file")
parser.add_argument('-vdb', default='/Users/Sidney/nextstrain/fauna/data/dengue_all.fasta', type=str, help="Full path of vdb dengue fasta file")
parser.add_argument('-src', default = 'agm_1mo', type=str,  help="source designation")
args = parser.parse_args()

titers_df = pd.read_csv(args.titers, index_col=0, comment='#', na_values='*')
titers_df.dropna(how='all', inplace=True)

vdb_synonyms = { line.strip().split('\t')[0]:line.strip().split('\t')[1:] for line in open(args.synonyms, 'r') } # { vdbname: [smith1, smith2]}
smith_synonyms = {}	# { smith1: vdbname, smith2:vdbname, ...}
for v, slist in vdb_synonyms.items():
	for s in slist:
		smith_synonyms[s] = v

sequences = { smith_synonyms[s.description]: str(s.seq) for s in SeqIO.parse(open(args.sequences, 'r'), 'fasta') }
accessions = { smith_synonyms[line.split()[2]]: line.split()[3] for line in open(args.key, 'r') if 'fullname' not in line }
extra_accessions = [('DEN1_Cambodia_2003_GenBankGQ868619', 'GQ868619'),('DEN1_Burma_2005_61117', 'KT452791'),('DEN2_Vietnam_2003_DF_670','KT452797'),('DEN3_Burma_2008_80931', 'KT452792'),('DEN1_Bolivia_2010_FSB_3363', 'KT382187'),('DEN4_Burma_2008_81087', 'KT452793'),('DENV2/VIETNAM/AC212/2003', 'KT452796'),('DENV3/Fiji/1992-29472', 'L11422'),('DEN3_Fiji_1992__','L11422'),('DEN2_Senegal_2003_Sendak_HD_0674','EF105384'),('DENV2/Tonga/1974-Tonga-74','AY744147'),('DENV3/Fiji/1992-29472', 'L11422'),('DENV2/SENEGAL/SENDAKHD0674/1970', 'EF105384'),('DENV2/Senegal/1970/Sendak_HD_0674', 'EF105384'),('DENV2/SENEGAL/1970/SENDAKHD0674', 'EF105384')]
for strain, a in extra_accessions:
	try:
		newstrain = smith_synonyms[strain]
	except:
		print 'no strain synonym', strain
		continue
	try:
		accessions[newstrain] = a
	except:
		print 'no accession', newstrain

vdb = [ s.description.split('|')[0] for s in SeqIO.parse(open(args.vdb, 'r'), 'fasta') ]
vaccine_names = {'DEN4':'DENV4/DOMINICA/814669DELTA30/1981','DEN1':'DENV1/NAURU/WESTERNPACIFICDELTA30/1974','DEN3':'DENV3/INDONESIA/SLEMANDELTA30/1978','DEN2':'DENV2/TONGA/DELTA30/1974'}

####    Pull existing titer measurements (if any)   #####

existing_titers = glob('dengue_strains.tsv')
if existing_titers != []:
	existing_titers = { line.split('\t')[0]:int(line.split('\t')[1]) for line in open(existing_titers[0], 'r') }

try:
	seen_smith_viruses = glob('smith_viruses.tsv')
	seen_smith_viruses = pd.read_csv(seen_smith_viruses[0], dtype='str', delimiter='\t', index_col=0, header=0, skiprows=1)
	seen_smith_viruses = list(seen_smith_viruses['Accession'])
except:
	seen_smith_viruses = []

####    Define functions        ####

def format_smith_upload(missing_strains, no_append=True):
	'''
	Given strain list, pull metadata from strain name and put into LANL format
	to make uploading missing records to vdb easy.
	Write this to tsv and print suggested upload command.
	Return list of any strains that weren't found in the fasta and/or key.
	'''
	records = []
	for s in missing_strains:
		accession = accessions[s]
		if accession in seen_smith_viruses:
			continue
		else:
			sequence = sequences[s]
			seen_smith_viruses.append(accession)

		row = {}
		NA_fields = ['Species', 'Isolate Name', 'Georegion', 'Author', 'Sampling City', 'Pubmed ID']
		for n in NA_fields:
			row[n] = None
		row['Serotype'], row['Country'], row['Strain'], row['Sampling Year'] = s.split('/')
		row['Accession'] = accession
		row['Name'] = s
		row['Start'] = 935
		row['Stop'] = 2413
		row['Segment'] = 'E'
		row['Organism'] = 'dengue virus '+row['Serotype'][-1]
		row['Species'] = 'dengue virus'
		row['Sequence'] = sequence
		records.append(row)
	if no_append:
		open('smith_viruses.tsv', 'a').write('# viruses from smith data not in vdb\n')
	pd.DataFrame(records).to_csv('smith_viruses.tsv', sep='\t', mode='a', header=no_append)

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

def convert_smith_vaccine_strains(titerdf):
	titer_file = open('dengue_titers.tsv', 'a')
	strain_file = open('dengue_strains.tsv', 'w')
	value_counts = defaultdict(int, existing_titers)
	for virus, seraseries in titerdf.iterrows():
		for sera, value in seraseries.iteritems():
			smith_fields = sera.split('_')
			sera = vaccine_names[smith_fields[0]]
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

####    Find strains with existing VDB records      ####
if 'monovalent' in args.src:
	titers_df = titers_df.rename(index = smith_synonyms)
	convert_smith_vaccine_strains(titers_df)
	missing_strains = list(set([ v for v in list(titers_df.index.values)+vaccine_names.values() if v not in vdb ]))
else:
	titers_df = titers_df.rename(index = smith_synonyms, columns = smith_synonyms)   # Relabel titer data
	table_to_tsv(titers_df)                                                     # Write to augur-friendly files
	missing_strains = list(set([ v for v in list(titers_df.index.values)+list(titers_df.columns.values) if v not in vdb ]))

####    Deal with records not in VDB    ####
print 'These strains were not found in vdb, but were in provided fasta file.\nTo upload, run `fauna$ python vdb/dengue_upload --fname smith_viruses.tsv --ftype tsv -v dengue -db vdb`\n'
format_smith_upload(missing_strains, no_append=(seen_smith_viruses==[]))
pprint(sorted(missing_strains))

#####	Record how we formatted strains per accession numbers		######
fauna_key = { line.split()[1]: line.split()[0] for line in open('/Users/Sidney/nextstrain/fauna/source-data/dengue_strain_name_fix.tsv', 'r')}
for s, a in accessions.iteritems():
	fauna_key[s] = a
fauna_key_file = open('/Users/Sidney/nextstrain/fauna/source-data/dengue_strain_name_fix.tsv', 'w')
for s, a in sorted(fauna_key.items()):
	fauna_key_file.write(a+'\t'+s+'\n')
fauna_key_file.close()
