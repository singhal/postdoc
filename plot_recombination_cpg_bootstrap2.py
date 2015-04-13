import re
import glob
import pandas as pd
import numpy as np
import sys
from itertools import izip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
sp = args.sp

dir = '/mnt/gluster/home/sonal.singhal1/%s/analysis/TSS/cpg/' % sp
out = '%ssummary_bootstrap.csv' % (dir)
boot_files = glob.glob('%sboot*csv' % dir)

values = {}
for file in boot_files:
	print file
	d = pd.read_csv(file)
	d = d.rename(columns={'mean': 'mean_rho'})
	for loc, bin, mean in izip(d.location, d.cpg_bin, d.mean_rho):
		if loc not in values:
			values[loc] = {}
		if bin not in values[loc]:
			values[loc][bin] = {'min': mean, 'max': mean}
		if mean < values[loc][bin]['min']:
			values[loc][bin]['min'] = mean
		if mean > values[loc][bin]['max']:
                        values[loc][bin]['max'] = mean

o = open(out, 'w')
o.write('location,cpg_bin,min,max\n')
for loc in values:
	for bin in values[loc]:
		o.write('%s,%s,%s,%s\n' % (loc, bin, values[loc][bin]['min'], values[loc][bin]['max']))
