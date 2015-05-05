import re
import gzip
import pandas as pd
import argparse
from itertools import izip
import subprocess
import os
import glob

def parse_files(out_files, out):
	vals = {}
	for type1 in ['exon', 'intron']:
		vals[type1] = {'0': {}, '1': {}, '2': {}, '3': {}, '4': {}, 'greater': {}}
		for type2 in vals[type1]:
			for i in range(0,101):
				vals[type1][type2][i / 100.] = {'rho': 0, 'num': 0}
	for file in out_files:
		f = open(file, 'r')
		header = f.next()
		for l in f:
			d = re.split(',', l.rstrip())
			rel_pos = float(d[1])
			exon_intron = d[2]
			num = int(d[3])
			tot_num = int(d[4])
			if tot_num >= 4:
				type = str(num)
				if num > 4:
					type = 'greater'
				vals[exon_intron][type][rel_pos]['rho'] += float(d[5])
				vals[exon_intron][type][rel_pos]['num'] += 1
		f.close()
	
	o = open(out, 'w')
	o.write('exon_intron,location,relative_position,rho,num_sites\n')
	for type1 in vals:
		for type2 in vals[type1]:
			for pos in vals[type1][type2]:
				rho = vals[type1][type2][pos]['rho']
				num = vals[type1][type2][pos]['num']
				rho = rho / float(num)
				o.write('%s,%s,%s,%s,%s\n' % (type1, type2, pos, rho, num))
	o.close()

def main():
        parser = argparse.ArgumentParser()
        parser.add_argument("--sp", help="chromosome for which to run analysis")
        args = parser.parse_args()
        sp = args.sp

	# out files
        out_files = glob.glob('/mnt/gluster/home/sonal.singhal1/%s/analysis/gene_rho/gene_rho.*.csv' % (sp))
	out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/gene_rho/gene_rho_long_genes_5_cats.csv' % sp
	parse_files(out_files, out)
	print 'yay, finished'

if __name__ == "__main__":
    main()
