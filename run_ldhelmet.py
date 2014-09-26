import re
import subprocess

chrs = [['chr1', 'chr1A', 'chr1B'], ['chr2', 'chr3'], ['chr4', 'chr5'], ['chr6', 'chr7'],
	['chr8', 'chr4A', 'chr9'], ['chr10', 'chr11', 'chr12', 'chr13', 'chr14'], ['chr15',
	'chr16', 'chr17', 'chr18', 'chr19'], ['chr20', 'chr21', 'chr22', 'chr23', 'chr24',
	'chr25'], ['chr26', 'chr27', 'chr28', 'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ']]
bpen = 100


for species in ['ZF', 'LTF']:
	pade = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/%s.pade' % (species, species)
	lk = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/%s.lk' % (species, species)
	mm = '/mnt/gluster/home/sonal.singhal1/%s/analysis/mutation_matrix.txt' % (species)

	for ix, chr_list in enumerate(chrs):
		outsh = '%s_%s.sh' % (species, ix)
		o = open(outsh, 'w')
		for chr in chr_list:
			aa = '/mnt/gluster/home/sonal.singhal1/%s/ancestral_allele/ancestral_allele.%s.ldhelmet.txt' % (species, chr)
			haplo = '/mnt/gluster/home/sonal.singhal1/%s/phasing/PIR_approach/%s_haplotypes.fasta' % (species, chr)
			out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/%s_recombination_bpen%s' % (species, chr, bpen)
	
			call = '/mnt/lustre/home/sonal.singhal1/bin/LDhelmet_v1.6/ldhelmet rjmcmc --num_threads 12 -o %s -n 3000000 --burn_in 300000 -b %s -s %s -l %s -p %s -a %s -m %s -w 50' % (out, bpen, haplo, lk, pade, aa, mm)
			o.write(call + '\n')
		o.close()
		subprocess.call("chmod a+x %s" % outsh, shell=True) 

