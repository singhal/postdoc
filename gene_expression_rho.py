import pandas as pd
import numpy as np
from itertools import izip

recom_file = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/maps/chr%s_recombination_bpen100.txt'
gff = '/mnt/gluster/home/sonal.singhal1/reference/Taeniopygia_guttata.gff'
# gexp = '/mnt/gluster/home/sonal.singhal1/reference/gexp/Balakrishnan2015testes_gexp.csv'
gexp = '/mnt/gluster/home/sonal.singhal1/reference/gexp/Balakrishnan2015_gexp.csv'

chrs = ['1', '1A', '2', '3', '4', '4A', '5', '6', \
	'7', '8', '9', '10', '11', 'chr12', '13', '14', '15']
allchrs =  [	'1', '1A', '1B', '2', '3', '4', '4A', '5', '6', \
        	'7', '8', '9', '10', '11', '12', '13', '14', '15',
		'17', '18', '19', '20', '21', '22', '23', '24', '25', '26',
		'27', '28', 'LG2', 'LGE22', 'Z']


gexp_f = pd.read_csv(gexp)
gexp = {}
for id, fpkm in izip(gexp_f.gene, gexp_f.fpkm):
	gexp[id] = fpkm 
# for id, hfpkm, tfpkm, afpkm in izip(gexp_f.gene, gexp_f.hybrid_fpkm, gexp_f.timor_fpkm, gexp_f.all_fpkm):
#	gexp[id] = {'hybrid': hfpkm, 'timor': tfpkm, 'all': afpkm}
gexp_f = ''

# need to identify location of first gene and span
gff = pd.read_csv(gff, sep='\t', header=None, names=['chr', 'type', 'cds_mrna', 
                                                'start', 'stop', 'score', 
                                                'orientation', 'codon_pos', 'id'])
gff['id'] = [x.replace('ID=', '') for x in gff.id]
gff['id'] = [x.replace('Parent=', '') for x in gff.id]
gff['id'] = [x.replace(';', '') for x in gff.id]
#gff = gff[gff.chr.isin(chrs)]

# chr -> gene -> first exon & entire length & orient
genes = {}
mrna = gff[gff.cds_mrna == 'mRNA']
for chr, gene, s, e, orient in izip(mrna.chr, mrna.id, mrna.start, mrna.stop, mrna.orientation):
	if chr not in genes:
		genes[chr] = {}
	genes[chr][gene] = {'ends': [s, e], 'exon': None, 'orient': orient}
cds = gff[gff.cds_mrna == 'CDS'].groupby(['chr', 'id'])

for (chr, gene), group in cds:
	if genes[chr][gene]['orient'] == '+':
		genes[chr][gene]['exon'] = [group.start.tolist()[0] , group.stop.tolist()[0]]
	else:
		genes[chr][gene]['exon'] = [group.start.tolist()[-1] , group.stop.tolist()[-1]]
gff = ''
cds = ''
mrna = ''

def get_rhos(rho, chr, gene, genes, type):
	chunk_start = genes[chr][gene][type][0]
	chunk_end = genes[chr][gene][type][1]
	tmp = rho[rho.right_snp >= chunk_start]
	tmp = tmp[tmp.left_snp <= chunk_end]

	rate_mult = 0
	diff = 0

	for start, end, rate in izip(tmp.left_snp, tmp.right_snp, tmp.meanrho):
		if start >= chunk_start and end <= chunk_end:
			rate_mult += (end - start) * rate
			diff += (end - start)
		elif start < chunk_start and (end <= chunk_end and end >= chunk_start):
			rate_mult += (end - chunk_start) * rate
			diff += (end - chunk_start)
		elif (start >= chunk_start and start <= chunk_end) and end > chunk_end:
			rate_mult += (chunk_end - start) * rate
			diff += (chunk_end - start)
		elif start <= chunk_start and end >= chunk_end:
			diff += (chunk_end - chunk_start)
			rate_mult += (chunk_end - chunk_start) * rate

	rho = 'NA'
	if diff > 0:
		rho = rate_mult / float(diff)

	return rho

out = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/TSS/gene_expression_rho.Balakrishnan2015.all.csv'
# out = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/TSS/gene_expression_rho.Balakrishnan2015_testes.all.csv'
o = open(out, 'w')
o.write('chr,gene,all_fpkm,zfexonrho,ltfexonrho,zfgenerho,ltfgenerho\n')
# o.write('chr,gene,hybrid_fpkm,timor_fpkm,all_fpkm,zfexonrho,ltfexonrho,zfgenerho,ltfgenerho\n')
for chr in allchrs:
	ltfrho = pd.read_csv(recom_file % ('LTF', chr), sep=" ", skiprows=3, header=None,
			names=['left_snp', 'right_snp', 'meanrho', 'p0.025', 'p0.975'])
	zfrho = pd.read_csv(recom_file % ('ZF', chr), sep=" ", skiprows=3, header=None,
			names=['left_snp', 'right_snp', 'meanrho', 'p0.025', 'p0.975'])

	for gene in genes[chr]:
		if gene in gexp:
			ltf_exon = get_rhos(ltfrho, chr, gene, genes, 'exon')
			zf_exon = get_rhos(zfrho, chr, gene, genes, 'exon')

			ltf_gene = get_rhos(ltfrho, chr, gene, genes, 'ends')
			zf_gene = get_rhos(zfrho, chr, gene, genes, 'ends')

			o.write('%s,%s,%s,%s,%s,%s,%s\n' % (chr, gene, gexp[gene], zf_exon, ltf_exon, zf_gene, ltf_gene))
			# o.write('%s,%s,%s,%s,%s,%s,%s,%s,%s\n' % (chr, gene, gexp[gene]['hybrid'], gexp[gene]['timor'], gexp[gene]['all'], zf_exon, ltf_exon, zf_gene, ltf_gene))
o.close()
