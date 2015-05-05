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
		vals[type1] = {'last': {}, 'middle': {}, 'first': {}}
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
			type = None
			num = int(d[3])
			tot_num = int(d[4])
			if (num + 1) == tot_num:
				type = 'last'
			elif num > 0:
				type = 'middle'
			if num == 0:
				type = 'first'
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
	out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/gene_rho/gene_rho_by_3_cats.csv' % sp
	parse_files(out_files, out)
	print 'yay, finished'

if __name__ == "__main__":
    main()
