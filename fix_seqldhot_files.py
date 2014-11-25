import glob
import re
import os

dir = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/hotspots/'
files = glob.glob(dir + '*txt')
# files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/putative_hotspot_chr1_68972000.seqLDhot.*txt')

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
		filein = file.replace('.txt', '.in')
		'''
		f = open(filein, 'r')
		filein2 = filein + '2'
		o = open(filein2, 'w')
		for l in f:
			if re.search('background rho', l):
				brho = float(re.search('= (\S+)', l).group(1))
				if brho == 0:
					o.write('background rho = 0.0001\n')
				else:
					o.write(l)
			else:
				o.write(l)
		f.close()
		o.close()
		os.rename(filein2, filein)
		'''		
	
		print '/mnt/lustre/home/sonal.singhal1/bin/sequenceLDhot/sequenceLDhot %s %s' % (filein, file)
		
