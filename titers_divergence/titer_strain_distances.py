from Bio import SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from glob import glob
import os
import pandas as pd

with open('/Users/Sidney/nextstrain/fauna/data/dengue_titers.tsv', 'r') as f: # Pull titers from fauna/data directory
	titerstrains = set([ line.split()[0] for line in f ])
with open('/Users/Sidney/nextstrain/fauna/data/dengue_titers.tsv', 'r') as f:
	serastrains = set([ line.split()[1] for line in f ])

autologous = titerstrains.intersection(serastrains) # Find list of strains that have autologous titers

# Pull the sequence for all the strains with autologous titers
strains_with_titers = [s for s in SeqIO.parse(open('/Users/Sidney/nextstrain/fauna/data/dengue.fasta', 'r'), 'fasta') if s.description.split('|')[0] in autologous ]

# Write each sequence set to file
for sero in ['1', '2', '3', '4']:
	SeqIO.write([s for s in strains_with_titers if s.description.startswith('DENV%s'%sero)], './control_S%s.fasta'%sero, 'fasta')
SeqIO.write(strains_with_titers, './control_any.fasta', 'fasta')

# Align
for s in glob('./control*.fasta'):
	outfile = s.split('/')[-1].split('.')[0]+'_aln.fasta'
	os.system('mafft %s > %s'%(s, outfile)) # Default alignment in mafft
	os.system('rm %s'%s)

	aln = AlignIO.read(open(outfile, 'r'), 'fasta')
	distance_matrix = DistanceCalculator('identity').get_distance(aln) # Basic identity scoring
	strains = [s.split('|')[0] for s in distance_matrix.names] # Truncate to strain names
	distance_matrix = pd.DataFrame(distance_matrix.matrix, columns=strains, index=strains)
	outfile = outfile.split('_aln.')[0]+'_distMatrix.csv'
	distance_matrix.to_csv(outfile) # Write distance matrix to file
