from Bio import SeqIO
from glob import glob

def remove_no_date(infile):
	fstem = infile.split('/')[-1].split('.')[0]
	sequences = []
	for seq in SeqIO.parse(infile, 'fasta'):
		if len(seq.description.split('|')) < 2:
			print seq.description
		elif seq.description.split('|')[2].startswith('19') or seq.description.split('|')[2].startswith('20'):
			seq.description = seq.description.replace('-XX-', '').replace('XX|', '|').replace('-|', '|')		
			sequences.append(seq)
	print len(sequences)
	SeqIO.write(sequences, '/Users/Sidney/Dropbox/dengue/data/subsampled_alignments/%s_dated.fasta'%fstem, 'fasta')

files = glob('/Users/Sidney/Dropbox/dengue/data/subsampled_alignments/*.fasta')
print files
for f in files:
	remove_no_date(f)
