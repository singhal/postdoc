import glob
import re
import os

dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/hotspot_checks/phasing_uncertainty/'
files = glob.glob(dir + '*txt')

for file in files:
	out = file + '.sum'
	# assume that things didn't work
	good = False
	if os.path.isfile(out):
		# check if the file ran and is complete
		f = open(file, 'r')
		for l in f:
			if re.search('Positions', l):
				positions = f.next().rstrip()
				sites = [int(x) for x in re.split('\s+', positions)]
		f.close()
		
		last_snp = 0
		f = open(out, 'r')
		for l in f:
			if re.search('^\d', l):
				d = re.split('\s+', l.rstrip())
				if len(d) > 4:
					last_snp = int(d[-1])
		f.close()
		
		if len(sites) > 0:
			if abs(sites[-1] - last_snp) < 1000:
				good = True

	if not good:
		infile = file.replace('txt', 'in')
		print '/mnt/lustre/home/sonal.singhal1/bin/sequenceLDhot/sequenceLDhot %s %s' % (infile, file)
