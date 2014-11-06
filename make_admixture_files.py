import glob
import re
import gzip
import random
from itertools import izip

species = {	'DBF': '/mnt/gluster/home/sonal.singhal1/DBF/after_vqsr/by_chr/gatk.ug.dbf.%s.filtered.coverage.biallelicSNPs.vqsr.vcf.gz',
		'ZF': '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.filtered.recoded_biallelicSNPs.nomendel.vcf.gz',
		'LTF': '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.filtered.coverage.recoded_biallelicSNPs.vqsr.vcf.gz' }
sampling_freq = 400
chrs = ['chr1', 'chr1A', 'chr1B', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr4A', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28'] 
ind = {'DBF': 1, 'LTF': 20, 'ZF': 19}
# count = {'DBF': 10, 'LTF': 10, 'ZF': 10}
# chrs = ['chrLGE22', 'chrLG5']

def get_masked_genome(sp):
	gen_file = '/mnt/gluster/home/sonal.singhal1/%s/masked_genome/%s.masked_genome.fa' % (sp, sp)
	f = open(gen_file, 'r')
	return f

snps = {}
genome = []
for sp in species:
	genome.append(get_masked_genome(sp))
id = ''
pos = 0 
count = 0
for l1, l2, l3 in izip(genome[0], genome[1], genome[2]):
	l1 = l1.rstrip()
	l2 = l2.rstrip()
	l3 = l3.rstrip()

	if re.search('>', l1):
		id = re.search('>(\S+)', l1).group(1)
		snps[id] = {}
		pos = 0
		count = 0
	else:
		for bp1, bp2, bp3 in izip(list(l1), list(l2), list(l3)):
			pos += 1
			if bp1 == '1' or bp2 == '1' or bp3 == '1':
				if not count % sampling_freq:
					snps[id][pos] = 1
				count += 1
for gen in genome:
	gen.close()					

snps_id = {}
out_file = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/all_species.geno'
out = open(out_file, 'w')

for chr in chrs:
	snps_id[chr] = {}
	for sp in species:
		vcf = species[sp] % chr
		f = gzip.open(vcf)
		for l in f:
			if not re.search('^#', l):
				d = re.split('\t', l.rstrip())
				pos = int(d[1])
				if pos in snps[chr]:
					if pos not in snps_id[chr]:
						snps_id[chr][pos] = {'ZF': [], 'LTF': [], 'DBF': []}
					for geno in d[9:]:
						if re.match('0/0', geno):
							snps_id[chr][pos][sp].append(0)
						elif re.match('0/1', geno):
							snps_id[chr][pos][sp].append(1)
						elif re.match('1/1', geno):
							snps_id[chr][pos][sp].append(2)
						else:
							snps_id[chr][pos][sp].append(9)
		f.close()

for chr in snps_id:
	for pos in snps_id[chr]:
		for sp in snps_id[chr][pos]:
			if len(snps_id[chr][pos][sp]) == 0:
				snps_id[chr][pos][sp] = [0] * ind[sp]

for chr in snps_id:
	for pos in snps_id[chr]:
		snps = ''
		for species in ['DBF', 'ZF', 'LTF']:
			for bp in snps_id[chr][pos][species]:
				snps += str(bp)
		out.write(snps + '\n')
out.close()
