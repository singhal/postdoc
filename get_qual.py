import gzip
import numpy as np
import scipy as sci
import re
import pandas as pd
import argparse

chr_lengths = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351}

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
chr = args.chr

hot_file = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/spot2kb_flank40kb.seqldhot_validate_hotspots.csv'
d = pd.read_csv(hot_file)
d = d[d.zlk >= 10]
d = d[d.chr == chr]

out_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/hotspot_checks/data_quality/snp_qual/'
hot = open('%shot_%s.out' % (out_dir, chr), 'w')
out = open('%sout_%s.out' % (out_dir, chr), 'w')

vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/for_shapeit/gatk.ug.finch19.%s.allfilters.recoded_biallelicSNPs.nomissing.vcf.gz' % (chr)
if chr == 'chrZ':
	vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/for_shapeit/gatk.ug.finch19.chrZ.allfilters.recodedsex.recoded_biallelicSNPs.vcf.gz'

sites = {}
f = gzip.open(vcf, 'r')
for l in f:
	if not re.search('#', l):
		a = re.split('\t', l)
		alleles = [a[3]] + re.split(',', a[4])
		indel = False
		for allele in alleles:
			if len(allele) > 1:
				indel = True
		if not indel:
			if re.search('VQSLOD=[\d|\.|\-]+', l):
				vqslod = re.search('VQSLOD=([\d|\.|\-]+)', l).group(1)
				sites[ int(a[1]) ] = float(vqslod)
f.close()

for start in range(1, chr_lengths[chr], 50000):
	end = start + 50000
	tmp_sites = filter(lambda x: x >= start, sites.keys())
	tmp_sites = filter(lambda x: x < end, tmp_sites)

	values = [sites[x] for x in tmp_sites]
	if len(values) > 0:
		out.write('%s\n' % np.mean(values))
	
for hot_start in d.zstart:
	start = hot_start - 25000
	if start < 1:
		start = 1
	end = hot_start + 25000
	if end > chr_lengths[chr]:
		end = chr_lengths[chr]

	tmp_sites = filter(lambda x: x >= start, sites.keys())
	tmp_sites = filter(lambda x: x < end, tmp_sites)

	values = [sites[x] for x in tmp_sites]
	if len(values) > 0:
		hot.write('%s\n' % np.mean(values))
hot.close()
out.close()
