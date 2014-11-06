import re
import subprocess

chrs = [['chr1', 'chr1A', 'chr1B'], ['chr2', 'chr3'], ['chr4', 'chr5'], ['chr6', 'chr7'],
	['chr8', 'chr4A', 'chr9'], ['chr10', 'chr11', 'chr12', 'chr13', 'chr14'], ['chr15',
	'chr16', 'chr17', 'chr18', 'chr19'], ['chr20', 'chr21', 'chr22', 'chr23', 'chr24',
	'chr25'], ['chr26', 'chr27', 'chr28', 'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ']]
bpen = 5

rhos = {'chrLG2': 0.5, 'chrLG5': 0.5, 'chrLGE22': 0.5, 'chr22': 0.5, 'chr1': 0.021, 'chr10': 0.099, 'chr11': 0.058, 'chr12': 0.069, 'chr13': 0.030, 'chr14': 0.116, 'chr15': 0.127, 'chr16': 1.0, 'chr17': 0.182, 'chr18': 0.146, 'chr19': 0.181, 'chr1A': 0.042, 'chr1B': 1.0, 'chr2': 0.012, 'chr20': 0.129, 'chr21': 0.269, 'chr23': 0.330, 'chr24': 0.046,  'chr25': 1.129, 'chr26': 0.413, 'chr27': 0.096, 'chr28': 0.000, 'chr3': 0.021, 'chr4': 0.009, 'chr4A': 0.097, 'chr5': 0.033, 'chr6': 0.060, 'chr7': 0.039, 'chr8': 0.060, 'chr9': 0.071, 'chrZ': 0.014}

for species in ['ZF', 'LTF']:
	pade = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/%s.pade' % (species, species)
	lk = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/%s.lk' % (species, species)
	mm = '/mnt/gluster/home/sonal.singhal1/%s/analysis/mutation_matrix.txt' % (species)

	for ix, chr_list in enumerate(chrs):
		outsh = '%s_%s_%s.sh' % (species, bpen, ix)
		o = open(outsh, 'w')
		for chr in chr_list:
			aa = '/mnt/gluster/home/sonal.singhal1/%s/ancestral_allele/ancestral_allele.%s.ldhelmet.txt' % (species, chr)
			haplo = '/mnt/gluster/home/sonal.singhal1/%s/phasing/PIR_approach/%s_haplotypes.fasta' % (species, chr)
			out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/%s_recombination_bpen%s' % (species, chr, bpen)
			prior = rhos[chr]
			max_lk = prior * 2	

			call = '/mnt/lustre/home/sonal.singhal1/bin/LDhelmet_v1.6/ldhelmet rjmcmc --num_threads 12 -o %s -n 1000000 --burn_in 100000 -b %s -s %s -l %s -p %s -a %s -m %s -w 50 --max_lk_end %s --prior_rate %s' % (out, bpen, haplo, lk, pade, aa, mm, max_lk, prior)
			o.write(call + '\n')
		o.close()
		subprocess.call("chmod a+x %s" % outsh, shell=True) 

