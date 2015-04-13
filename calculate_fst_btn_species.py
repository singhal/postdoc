import re
import gzip
import argparse
import subprocess
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
chr = args.chr

'''
ZF Files
'''
vcf1 = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.shared.noswitch.phased.vqsr2.vcf.gz' % chr
if chr == 'chrZ':
	vcf1 = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.chrZ.coverage.repeatmasked.filtered.nomendel.shared.recodedsex.phased.vqsr2.vcf.gz'
zf_genome = '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.repeat_masked.switch_masked.fa'

'''
LTF Files
'''
vcf2 = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.coverage.repeatmasked.filtered.phased.vqsr2.vcf.gz'  % chr
if chr == 'chrZ':
	vcf2 = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chrZ.coverage.repeatmasked.filtered.recodedsex.phased.vqsr2.vcf.gz'
ltf_genome = '/mnt/gluster/home/sonal.singhal1/LTF/masked_genome/LTF.masked_genome.repeat_masked.fa'

out = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/fst/%s.Wrights_Fst.ZF_LTF.csv.gz' % chr


def get_sequence(genome1, genome2, chr):
        seq_call = subprocess.Popen('samtools faidx %s %s' % (genome1, chr), shell=True, stdout=subprocess.PIPE)
        seq1 = ''
        for l in seq_call.stdout:
                if not re.match('>', l):
                        seq1 += l.rstrip().upper()
        seq1 = list(seq1)

        seq_call = subprocess.Popen('samtools faidx %s %s' % (genome2, chr), shell=True, stdout=subprocess.PIPE)
        seq2 = ''
        for l in seq_call.stdout:
                if not re.match('>', l):
                        seq2 += l.rstrip().upper()
        seq2 = list(seq2)

	ditch_sites = {}
	for ix, (a, b) in enumerate(zip(seq1, seq2)):
		pos = ix + 1
		if a == '8' or b == '8':
			ditch_sites[pos] = 1

	return(ditch_sites)


def get_vcf(vcf, ditch_sites):
	var = {}
	f = gzip.open(vcf, 'r')
	for l in f:
		if not re.search('^#', l):
			d = re.split('\t', l.rstrip())
			
			alleles = [d[3]] + re.split(',', d[4])
			pos = int(d[1])
			indel = False
			
			for allele in alleles:
				if len(allele) > 1:
					indel = True

		
			if not indel and pos not in ditch_sites:
				genos = []
				for geno in d[9:]:
					geno = re.search('^([^:]+)', geno).group(1)
                                        genos += re.split('[|/]', geno)
				genos = [x for x in genos if not re.search('\.', x)]

				sites = {}
				for ix, allele in enumerate(alleles):
					count = genos.count(str(ix))
					if count > 0:
						sites[allele] = count
				if len(sites) > 0:
					var[pos] = sites

	return var


def calc_fst(pop1_p, pop2_p, tot_p):
	hs_1 = 2 * pop1_p * (1 - pop1_p)
	hs_2 = 2 * pop2_p * (1 - pop2_p)
	ht = 2 * tot_p * (1 - tot_p)

	if ht != 0:
		fst = (ht - ((hs_1 + hs_2) / 2.0))/ ht
		return fst
	else:
		return None


def get_fst(chr, var1, var2, nsam1, nsam2, o):
	sites = sorted(set(var1.keys() + var2.keys()))
	fst = None
	for site in sites:
		if site in var1 and site in var2:
			alleles = [allele for allele in var1[site] if allele in var2[site]]
			if len(alleles) > 0:
				pop1_p = var1[site][alleles[0]] / float(np.sum(var1[site].values()))
				pop2_p = var2[site][alleles[0]] / float(np.sum(var2[site].values()))
				tot_p = (var1[site][alleles[0]] + var2[site][alleles[0]]) / float(np.sum(var1[site].values()) + np.sum(var2[site].values()))
				fst = calc_fst(pop1_p, pop2_p, tot_p)
		else:
			if site in var1:
				pop1_p = np.max(var1[site].values()) / float(np.sum(var1[site].values()))
				pop2_p = 0				
				tot_p = np.max(var1[site].values()) / float(np.sum(var1[site].values()) + nsam2)
				fst = calc_fst(pop1_p, pop2_p, tot_p)
			if site in var2:
				pop2_p = np.max(var2[site].values()) / float(np.sum(var2[site].values()))
                                pop1_p = 0
                                tot_p = np.max(var2[site].values()) / float(np.sum(var2[site].values()) + nsam1)
				fst = calc_fst(pop1_p, pop2_p, tot_p)

		if fst:
			o.write('%s,%s,%s\n' % (chr, site, fst))


ditch_sites = get_sequence(zf_genome, ltf_genome, chr)
var1 = get_vcf(vcf1, ditch_sites)
var2 = get_vcf(vcf2, ditch_sites)
o = gzip.open(out, 'w')
o.write('chr,pos,fst\n')
get_fst(chr, var1, var2, 38, 40, o)
o.close()
