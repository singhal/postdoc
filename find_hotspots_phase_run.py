import re
import glob
import numpy as np
import pandas as pd
import random
from itertools import izip
import argparse
import gzip
import os

parser = argparse.ArgumentParser()
parser.add_argument("--sp", help="species for which to run the analysis")
args = parser.parse_args()
sp = args.sp

putative_hotspots = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/LTF_ZF.putative_hotspots.csv'
chrs = ['chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', \
		'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr27', 'chrZ']

# take this much sequence around the putative hotspot
block = 50e3
results_dir = '/mnt/gluster/home/sonal.singhal1/%s/analysis/phase_hotspots/' % sp

# get the most conservative set of hotspots 
d = pd.read_csv(putative_hotspots)
d = d[(d.spot_size == 2000) & (d.flank_size == 40000)]
d = d[d.chr.isin(chrs)]

# now start to create the files
for chr, chrhot in d.groupby('chr'):
	for start, length in zip(chrhot.spot_start, chrhot.length):
		out = '%sputative_hotspot_%s_%s.PHASE.txt' % (results_dir, chr, start)
		out_in = '%sputative_hotspot_%s_%s.PHASE.prior' % (results_dir, chr, start)     
		out_file = out.replace('.txt', '.asphased.out')
		out_hot = out_file + '_hotspot'

		run = True		
		if os.path.isfile(out_hot):
			if os.path.getsize(out_hot) > 0:
				run = False

		if run:
			start_hot = int(block / 2.0) - 50
			stop_hot = start_hot + length + 50
			print '/mnt/lustre/home/sonal.singhal1/bin/phase.2.1.1/PHASE -X100 -k999 -MR2 1 %s %s -r%s %s %s 100 1 100' % (start_hot, stop_hot, out_in, out, out_file)
