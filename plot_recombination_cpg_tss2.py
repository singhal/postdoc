import glob
import numpy as np
import re
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
sp = args.sp

dir = '/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/' % sp
files = glob.glob('%s*recombination_tss_cpg.csv.gz' % dir)

out = '%ssummary_tss_vs_cpg.csv' % dir

summary = {}

for file in files:
	f = gzip.open(file, 'r')
	header = f.next()

	for l in f:
		d = re.split(',', l.rstrip())
		rho = float(d[2])
		
		tss_bin	= ((int(d[3]) / 500) + 1) * 500
		
		cpg_bin = ((int(d[4]) / 500) + 1) * 500

		if tss_bin > 1e5:
			tss_bin = 1000000
		if tss_bin < -1e5:
			tss_bin = -1000000
		if cpg_bin > 1e5:
			cpg_bin = 1000000

		if tss_bin not in summary:
			summary[tss_bin] = {}
		if cpg_bin not in summary[tss_bin]:
			summary[tss_bin][cpg_bin] = {'sum_rhos': 0, 'number': 0}
		summary[tss_bin][cpg_bin]['sum_rhos'] += rho
		summary[tss_bin][cpg_bin]['number'] += 1
	f.close()

o = open(out, 'w')
o.write('tss_bin,cpg_bin,rho,number\n')
for tss_bin in summary:
	for cpg_bin in summary[tss_bin]:
		rho = summary[tss_bin][cpg_bin]['sum_rhos'] / float(summary[tss_bin][cpg_bin]['number'])
		o.write('%s,%s,%s,%s\n' % (tss_bin, cpg_bin, rho, summary[tss_bin][cpg_bin]['number']))
o.close()
		
