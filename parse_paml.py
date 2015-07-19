import re
import glob

out_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/dnds/results/'
files = glob.glob('%s*out' % out_dir)

print 'gene,lnl,sig,num_sites,p1,w1,p2,w2,p3,w3'
for file in files:
	gene = re.search('([^\/]+).paml.out', file).group(1)
	f = open(file, 'r')

	lnls = []
	grid = False
	sites = []
	p = []
	w = []

	for l in f:
		l = l.rstrip()
		if re.search('^p:', l):
			p = re.split('\s+', l)[1:]
			w = re.split('\s+', f.next().rstrip())[1:]	
		if re.search('^lnL', l):
			d = re.split('\s+', l)
			lnls.append(float(d[-2]))
		if re.search('Bayes Empirical Bayes', l):
			while not grid:
				l = f.next().rstrip()
				if re.search('grid', l):
					grid = True
				d = re.split('\s+', l)
				if len(d) == 7:
					if d[0] == '' and re.search('\*', l):
						sites.append({d[1]: d[3]})
			
		
	if len(lnls) == 2:
		chi_square = 2 * abs(lnls[0] - lnls[1])		
		# with df = 2, chi_square sig = 5.991
		sig = False
		if chi_square >= 5.991:
			sig = True
		joined = ','.join(['%s,%s' % (a,b) for a, b in zip(p, w)])
		print '%s,%.2f,%s,%s,%s' % (gene, chi_square, sig, len(sites), joined)

