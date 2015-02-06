import re
import gzip
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chr for which to run analysis")
args = parser.parse_args()
chr = args.chr

zf_vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/gatk.ug.all_zf.%s.coverage.repeatmasked.filtered.nomendel.shared.noswitch.vqsr2.vcf.gz' % chr
ltf_vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.coverage.repeatmasked.filtered.vqsr2.vcf.gz' % chr

if chr == 'chrZ':
	zf_vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/gatk.ug.all_zf.chrZ.coverage.repeatmasked.filtered.nomendel.shared.recodedsex.vqsr2.vcf.gz'
	ltf_vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chrZ.coverage.repeatmasked.filtered.recodedsex.vqsr2.vcf.gz'

zf_genome = '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.repeat_masked.switch_masked.fa'
ltf_genome = '/mnt/gluster/home/sonal.singhal1/LTF/masked_genome/LTF.masked_genome.repeat_masked.fa'
out_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/shared_polymorphism/'
file1 = '%s%s_sites.csv' % (out_dir, chr)
sites_f = open(file1, 'w')
file2 = '%s%s_counts.csv' % (out_dir, chr)
counts_f = open(file2, 'w')

def get_chromosome(chr, genome):
	seq = subprocess.Popen('~/bin/samtools-0.1.19/samtools faidx %s %s' % (genome, chr), shell=True, stdout=subprocess.PIPE)
	chr = ''
	for l in seq.stdout:
		if not re.match('>', l):
			chr += l.rstrip()
	chr = [int(x) for x in list(chr)]
	return chr

def get_vcf(vcf, num):
	f = gzip.open(vcf, 'r')
	var = {}
	for l in f:
		if not re.match('^#', l):
			d = re.split('\t', l)
			if int(d[1]) not in var:
				alleles = [d[3]] + re.split(',', d[4])
				genos = []
				for geno in d[9:(num+9)]:
					geno = re.search('^([^:]+)', geno).group(1)
					genos += re.split('/', geno)
				genos = filter(lambda x: x != '.', genos)
				af = {}
				if len(genos) > 0:
					for ix, allele in enumerate(alleles):
						freq = genos.count(str(ix)) / float(len(genos))
						if freq > 0 and freq < 1:
							af[allele] = freq
					if len(af) > 1:
						var[int(d[1])] = dict()
						for allele, af in af.items():
							var[int(d[1])][allele] = af
	f.close()
	return var
					
zf_chr = get_chromosome(chr, zf_genome)
ltf_chr = get_chromosome(chr, ltf_genome)
zf_var = get_vcf(zf_vcf, 19)
ltf_var = get_vcf(ltf_vcf, 20)

sites_f.write('chr,pos,shared_data\n')
denom = 0
counts = {'indel': {'LTFunique': 0, 'ZFunique': 0, 'shared': 0}, 'snp':  {'LTFunique': 0, 'ZFunique': 0, 'shared': 0}}
for pos, (zfbp, ltfbp) in enumerate(zip(zf_chr, ltf_chr)):
	pos = pos + 1
	if zfbp < 1 and ltfbp < 1:
		denom += 1
		type = None
		shared = None
		if pos in zf_var:
			if pos in ltf_var:
				shared_alleles = list(set(ltf_var[pos].keys()) & set(zf_var[pos].keys()))
				if len(shared_alleles) > 1:
					type = 'snp'
					shared = 'shared'
					printalleles = []
					for var in shared_alleles:
						printalleles.append('%s|%.3f|%.3f' % (var, zf_var[pos][var], ltf_var[pos][var]))
						if len(var) > 1:
							type = 'indel'
					counts[type][shared] += 1
					sites_f.write('%s,%s,%s\n' % (chr, pos, ';'.join(printalleles)))
				else:
					# variant in both species, but not shared
					type = 'snp'
                                	shared = 'ZFunique'
                                	for var in zf_var[pos]:
                                        	if len(var) > 1:
                                        	        type = 'indel'
                                	counts[type][shared] += 1
	
					type = 'snp'
                                        shared = 'LTFunique'
                                        for var in ltf_var[pos]:
                                                if len(var) > 1:
                                                        type = 'indel'
                                        counts[type][shared] += 1	
			else:
				type = 'snp'
				shared = 'ZFunique'
				for var in zf_var[pos]:
                                        if len(var) > 1:
                                                type = 'indel'
				counts[type][shared] += 1
		else:
			if pos in ltf_var:
				# unique to ltf
				type = 'snp'
				shared = 'LTFunique'
				for var in ltf_var[pos]:
					if len(var) > 1:
						type = 'indel'
				counts[type][shared] += 1

counts_f.write('chr,denominator,snp_indel,type,count\n')
for type1 in counts:
	for type2 in counts[type1]:
		counts_f.write('%s,%s,%s,%s,%s\n' % (chr, denom, type1, type2, counts[type1][type2]))
counts_f.close()
sites_f.close()
