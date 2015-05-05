import re
import gzip
import pandas as pd
import argparse
from itertools import izip
import subprocess
import os

def get_genes(gff, chr):
        gff = pd.read_csv(gff, sep='\t', header=None, names=['chr', 'type', 'cds_mrna', 
                                                'start', 'stop', 'score', 
                                                'orientation', 'codon_pos', 'id'])
        gff['chr'] = 'chr' + gff.chr
        gff = gff[gff.chr == chr]
        gff['id'] = [x.replace('ID=', '') for x in gff.id]
        gff['id'] = [x.replace('Parent=', '') for x in gff.id]
        gff['id'] = [x.replace(';', '') for x in gff.id]

	features = {}

        genes = {}
        mrna = gff[gff.cds_mrna == 'mRNA']
        for id, orient in izip(mrna.id, mrna.orientation):
                genes[id] = {'orient': orient, 'exons': []}

        for id, group in gff[gff.cds_mrna == 'CDS'].groupby(['id']):
                for start, stop in izip(group.start, group.stop):
                        genes[id]['exons'].append([start,stop])
	
	for id in genes:
		if genes[id]['orient'] == '+':
			for ix, (exon_start, exon_end) in enumerate(genes[id]['exons']):
				length = exon_end - exon_start
				for pos_ix, pos in enumerate(range(exon_start, exon_end)):
					rel_pos = '%.2f' % (pos_ix / float(length))
					features[pos] = {'rel_pos': rel_pos, 'type': 'exon', 'num': ix, 'tot_num': len(genes[id]['exons'])}
			for ix, (first, second) in enumerate(zip(genes[id]['exons'], genes[id]['exons'][1:])):
				intron_start = first[1]
				intron_end = second[0]
				length = intron_end - intron_start
				for pos_ix, pos in enumerate(range(intron_start, intron_end)):
                                        rel_pos = '%.2f' % (pos_ix / float(length))
                                        features[pos] = {'rel_pos': rel_pos, 'type': 'intron', 'num': ix, 'tot_num': len(genes[id]['exons']) - 1}
		if genes[id]['orient'] == '-':
                        for ix, (exon_start, exon_end) in enumerate(genes[id]['exons']):
                                length = exon_end - exon_start
				rel_ix = len(genes[id]['exons']) - ix - 1
                                for pos_ix, pos in enumerate(range(exon_start, exon_end)):
                                        rel_pos = '%.2f' % (1 - (pos_ix / float(length)))
                                        features[pos] = {'rel_pos': rel_pos, 'type': 'exon', 'num': rel_ix, 'tot_num': len(genes[id]['exons'])}
                        for ix, (first, second) in enumerate(zip(genes[id]['exons'], genes[id]['exons'][1:])):
                                intron_start = first[1]
                                intron_end = second[0]
                                length = intron_end - intron_start
				rel_ix = (len(genes[id]['exons']) - 1) - ix - 1
                                for pos_ix, pos in enumerate(range(intron_start, intron_end)):
                                        rel_pos = '%.2f' % (1 - (pos_ix / float(length)))
                                        features[pos] = {'rel_pos': rel_pos, 'type': 'intron', 'num': rel_ix, 'tot_num': len(genes[id]['exons']) - 1}	
	genes = {}

        return features


def get_rhos(out_file, rho_file, features):
	f = open(rho_file, 'r')
	out = open(out_file, 'w')

	out.write('position,relative_position,exon_intron,number,total_number,rho\n')

	for i in range(3):
		junk = f.next()
	for l in f:
		d = re.split('\s+', l.rstrip())
		start = int(d[0])
		end = int(d[1])
		rho = float(d[2])

		for bp in range(start, end):
			if bp in features:
				out.write('%s,%s,%s,%s,%s,%s\n' % (bp, features[bp]['rel_pos'], features[bp]['type'], features[bp]['num'], \
						features[bp]['tot_num'], rho))
	out.close()


def main():
        parser = argparse.ArgumentParser()
        parser.add_argument("--chr", help="chromosome for which to run analysis")
	parser.add_argument("--sp", help="chromosome for which to run analysis")
        args = parser.parse_args()
        chr = args.chr
	sp = args.sp

	# gff file
        gff = '/mnt/gluster/home/sonal.singhal1/reference/Taeniopygia_guttata.gff'

	# out dir
	out_file = '/mnt/gluster/home/sonal.singhal1/%s/analysis/gene_rho/gene_rho.%s.csv' % (sp, chr)

	# rho file
	rho_file = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/maps/%s_recombination_bpen5.txt' % (sp, chr)

	features = get_genes(gff, chr)
	get_rhos(out_file, rho_file, features)
	print 'yay, finished'

if __name__ == "__main__":
    main()
