import glob
import numpy
import re
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
sp = args.sp

files = glob.glob("/mnt/gluster/home/sonal.singhal1/%s/analysis/mut_skew/chr*" % sp)
out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/mut_skew/summary_mutskew_tss.csv' % sp

types = ['AT_GC', 'GC_AT', 'GC_GC', 'AT_AT']
skew = {}
for file in files:
	f = open(file,'r')
	header = f.next()
	for l in f:
		d = re.split(',', l.rstrip())

		dist = int(d[4])
		rho = float(d[5])
		af = float(d[3])
		type = d[2]

		if abs(dist) <= 100000:
			if af < 1 and af > 0:
				if dist not in skew:
					skew[dist] = {}
					for t in types:
						skew[dist][t] = {'rho': 0, 'num': 0, 'af': 0}
				skew[dist][type]['rho'] += rho
				skew[dist][type]['num'] += 1
				skew[dist][type]['af'] += af
	f.close()

o = open(out, 'w')
o.write('location,type,average_af,average_rho,number_pos\n')
for pos in skew:
	for type in skew[pos]:
		if skew[pos][type]['num'] > 0:
			avg_af = skew[pos][type]['af'] / float( skew[pos][type]['num'] )
			avg_rho = skew[pos][type]['rho'] / float( skew[pos][type]['num'] )
			o.write('%s,%s,%s,%s,%s\n' % (pos, type, avg_af, avg_rho, skew[pos][type]['num'] ))
o.close()
