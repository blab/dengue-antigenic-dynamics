from glob import glob
from Bio import SeqIO

infiles = glob('/Users/Sidney/Dropbox/dengue/data/subsampled_alignments/denv*.fasta')

for f in infiles:
    ofile = open(f.split('/')[-1].split('.')[0]+'_gard.fasta', 'w')
    seqs = []
    for s in SeqIO.parse(f, 'fasta'):
        if len(s.description.split(' ')) < 2 or len(s.description.split('|')) < 2:
            print s.description
            continue
        else:
            strain = s.description.split(' ')[0]
            acc = s.description.split('|')[1]
            ofile.write('>'+strain+'/'+ acc + '\n' +str(s.seq)+'\n')
    #         seqs.append(s)
    # SeqIO.write(seqs, ofile_name, 'fasta')
