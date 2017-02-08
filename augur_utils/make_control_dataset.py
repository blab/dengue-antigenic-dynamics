from Bio import SeqIO
from pprint import pprint

records = open('../../../fauna/data/dengue_titers.tsv', 'r').readlines()
titerstrains = set(sorted([ l.split()[0] for l in records ]))
serastrains = set(sorted([ l.split()[1] for l in records ]))
wantstrains = titerstrains.union(serastrains)

fasta = { s.description.split('|')[0] :s  for s in SeqIO.parse(open('../../../fauna/data/dengue.fasta', 'r'), 'fasta')}

want_seqs = [ fasta[s] for s in wantstrains ]

SeqIO.write(want_seqs, 'control.fasta', 'fasta')
