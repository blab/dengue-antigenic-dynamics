from Bio import AlignIO
import numpy as np
import pandas as pd
from collections import defaultdict
from glob import glob
from matplotlib import pyplot as plt
import seaborn as sns


def plot_r_squared(filename):
	filestem = filename.split('/')[-1].split('.')[0]
	align = AlignIO.read(filename, 'fasta')
	align = np.array([list(rec) for rec in align], np.character, order="F") # Read in alignment, make a dataframe
	nsites = len(align.T)
	nseqs = len(align)
	align = pd.DataFrame(align, columns=range(nsites), index=range(nseqs)) # Use pandas indexing to keep site numbers in tact
	print 'alignment ready'

	nt = ['a', 'c', 't', 'g', 'A', 'C', 'T', 'G'] # Allowed nucleotide characters

	def get_majorminor(site):
		site = site[site.isin(nt)]           # Remove the gaps and other invalid characters
		counts = site.value_counts()
		total = float(counts.sum())
		if total < 20 or len(counts) < 2: # If not polymorphic, invalid.
			return np.nan                   # If < 10 characters left, invalid.
		counts = counts[ counts/total >= 0.01] # Remove rare variants
		if len(counts) == 2:
			return tuple(counts.index)    # If that leaves us with two alleles (major and minor), return just the alleles
		else:
			return np.nan                   # Otherwise, the site is not biallelic; return None.

	print 'finding major/minor alleles'
	majorminor = align.apply(lambda x: get_majorminor(x)) # Get the (major,minor) alleles for each site if possible
	majorminor.dropna(inplace=True) # Drop the non-biallelic sites
	validsites = majorminor.index.values # Valid sites are the biallelic ones
	results = pd.DataFrame(columns=validsites, index=validsites, dtype='float') # Store results

	print 'found %d valid sites'%len(validsites)
	print 'starting comparisons'

	def r_squared(haplotypes):
		hapObs = pd.Series(haplotypes).value_counts() # Tally observed counts of haplotypes and alleles
		if len(haplotypes) < 20: # Check, yet again, that we have enough data for this site
			return np.nan
		else:
			pass
	
		allele_i_obs = pd.Series([ h[0] for h in haplotypes]).value_counts()
		allele_j_obs = pd.Series([ h[1] for h in haplotypes ]).value_counts()
	
		if len(allele_i_obs) != 2 or len(allele_j_obs) != 2: # Check again that we have two alleles per site
			return np.nan
		else:
			pass
	

		hapExp = {}

		N = float(hapObs.sum()) # Total N of observed haplotypes
	
		for i in allele_i_obs.index.values:
			for j in allele_j_obs.index.values:
				# Exp counts = row (allele0|site0) total*column (allele1|site1) total / table (allele) total
				hapExp[i+j]=(float(allele_i_obs[i])*float(allele_j_obs[j]))/N # If this haplotype is never observed, use 0
			
		chisq = 0.0
		for hap in hapExp.keys(): # Calculate chi-squared
								  # chisq = sum( (E-O)**2/E for each possible haplotype )
			if hap in hapObs:     # deal with unobserved haplotypes
				chisq += ((float(hapObs[hap]) - hapExp[hap])**2)/hapExp[hap] 
			else:
				chisq += ((0.0 - hapExp[hap])**2)/hapExp[hap] 
	
		r_sq = chisq / (N)  # r^2 = chisq / (Nsamples * df) ; because all sites biallelic, df == (2-1)(2-1) == 1
		return r_sq
	print 'ready'

	for x in range(len(validsites)): # Make upper-triangle comparisons
		i = validsites[x]
		column_i = align.loc[:,i] # Pull the alignment column
		alleles_i = majorminor.loc[i] # And the (major,minor) alleles we found earlier
	
		for y in range(x+1, len(validsites)):
			j = validsites[y]
			column_j = align[j]
			alleles_j = majorminor.loc[j] # Make a list of all the observed haplotypes that consist of the major or minor allele for each site

			haplotypes = [ (h[0]+h[1]).upper() for h in zip(column_i, column_j) if h[0] in alleles_i and h[1] in alleles_j]

			if len(haplotypes) < 10: # Do we still have enough data? 
				continue
			else:
				results.at[i,j] = r_squared(haplotypes) # Save results in the df
	print 'done'


	results.dropna(axis=(0,1), how='all', inplace=True) # Drop any sites where we didn't have enough information to calculate chisqdf
	results.to_csv('r_squared_%s.csv'%filestem)

	plt.figure(figsize=(13,9))
	heatmap = sns.heatmap(results, cmap = plt.cm.GnBu, square = True, robust=True, annot=False)
	heatmap.axes.set_title('Pairwise R^2')
	plt.yticks(rotation=0)
	plt.xticks(rotation=90) 
	plt.tight_layout(pad=1)
	plt.savefig('r_squared_%s.png'%filestem)


file_list = glob('/Users/Sidney/nextstrain/dengue/recombination/*.fasta')
for f in file_list:
	plot_r_squared(f)

