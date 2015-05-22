import gzip
import re
import pandas as pd
import glob
import argparse
import subprocess

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chr")
args = parser.parse_args()
chr = args.chr

# raw snps
raw = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/finch.gatk18May.vcf.gz'

# /mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/finch.gatk18May.vcf.gz
file = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.shared.noswitch.phased.vqsr2.vcf.gz'  % chr
discord = '/mnt/lustre/home/emleffler/finches/zebrafinch/individual_stats/compare_duplicate/auts_diff_genotypes.txt'

if chr == 'chrZ':
	file = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.chrZ.coverage.repeatmasked.filtered.nomendel.shared.recodedsex.phased.vqsr2.vcf.gz'
	discord = '/mnt/lustre/home/emleffler/finches/zebrafinch/individual_stats/compare_duplicate/chrZ_diff_genotypes.txt'

f = gzip.open(file, 'r')
snps = {}
for l in f:
	if not re.search('^#', l):
		d = re.split('\t', l.rstrip())
		alleles = [d[3]] + re.split(',', d[4])
		indel = False
		for allele in alleles:
			if len(allele) > 1:
				indel = True
		if not indel:
			snps[int(d[1])] = 1
f.close()

data = {'tot_final_snps': len(snps), 'tot_raw_snps': 0, 'tot_discard': 0, 'tot_raw_match': 0, 'tot_discard_match': 0, 'tot_discard_missing': 0}

f = open(discord, 'r')
for l in f:
	d = re.split('\s+', l.rstrip())
	if d[0] == chr:
		data['tot_discard'] += 1
		if int(d[1]) in snps:
			if d[2] != './.' and d[3] != './.':
				data['tot_discard_match'] += 1
			else:
				data['tot_discard_missing'] += 1
f.close()

out = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/raw.%s.vcf' % chr
subprocess.call('tabix %s %s > %s' % (raw, chr, out), shell=True)
f = open(out, 'r')
for l in f:
        if not re.search('^#', l):
                d = re.split('\t', l.rstrip())
                alleles = [d[3]] + re.split(',', d[4])
                indel = False
                for allele in alleles:
                        if len(allele) > 1:
                                indel = True
                if not indel:
			data['tot_raw_snps'] += 1
			if int(d[1]) in snps:
				data['tot_raw_match'] += 1
f.close()

out = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/duplicate/%s.duplicate_genotypes.csv' % chr
o = open(out, 'w')
for key in data:
	o.write('%s,%s,%s\n' % (chr, key, data[key]))
o.close()



