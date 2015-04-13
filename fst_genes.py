import re
import numpy as np
from itertools import izip
import pandas as pd
import gzip
import os

gff = '/mnt/gluster/home/sonal.singhal1/reference/Taeniopygia_guttata.gff'
gff = pd.read_csv(gff, sep='\t', header=None, names=['chr', 'type', 'cds_mrna', 
                                                'start', 'stop', 'score', 
                                                'orientation', 'codon_pos', 'id'])
gff['id'] = [x.replace('ID=', '') for x in gff.id]
gff['id'] = [x.replace('Parent=', '') for x in gff.id]
gff['id'] = [x.replace(';', '') for x in gff.id]

genes = {}
for chr, group in gff[gff.cds_mrna == 'mRNA'].groupby('chr'):
	genes[chr] = {}
	for id, start, stop in izip(group.id, group.start, group.stop):
		genes[chr][id] = {'gene': [start, stop], 'exons': []}

for (chr, id), group in gff[gff.cds_mrna == 'CDS'].groupby(['chr', 'id']):
	for start, stop in izip(group.start, group.stop):
		genes[chr][id]['exons'].append([start,stop])

o = open('/mnt/gluster/home/sonal.singhal1/LTF/analysis/fst/ZF_LTF.gene_fst.csv', 'w')
o.write('chr,gene,gene_length,num_snps,gene_fst,exon_fst\n')
for chr in genes:
	fst_file = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/fst/chr%s.Wrights_Fst.ZF_LTF.csv.gz' % chr
	if os.path.isfile(fst_file):
		# get the snps
		snps = {}
		f = gzip.open(fst_file, 'r')
		heaader = f.next()
		for l in f:
			d = re.split(',', l.rstrip())
			fst = d[2]
			if fst != 'nan':
				if fst < 0:
					fst = 0
				snps[int(d[1])] = float(fst)
		f.close()

		snp_list = sorted(snps.keys())

		for gene in genes[chr]:
			gene_fst = 0
			exon_fst = 0
			num_exon_snps = 0

			gene_snps = [x for x in snp_list if x >= genes[chr][gene]['gene'][0]]
			gene_snps = [x for x in gene_snps if x <= genes[chr][gene]['gene'][1]]

			for snp in gene_snps:
				gene_fst += snps[snp]
			if len(gene_snps) > 0:
				gene_fst = gene_fst / float(len(gene_snps))

			for (exon_start, exon_end) in genes[chr][gene]['exons']:
				for snp in gene_snps:
					if snp >= exon_start and snp <= exon_end:
						exon_fst += snps[snp]
						num_exon_snps += 1
			if num_exon_snps > 0:
				exon_fst = exon_fst / float(num_exon_snps)

			gene_length = genes[chr][gene]['gene'][1] - genes[chr][gene]['gene'][0]

			o.write('%s,%s,%s,%s,%s,%s\n' % (chr, gene, gene_length, len(gene_snps), gene_fst, exon_fst))

o.close()	
