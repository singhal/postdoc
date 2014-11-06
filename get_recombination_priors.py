import pandas as pd

file = '/mnt/lustre/home/emleffler/finches/zebrafinch/downloads/linkage_map/linkage_map_rec_rate_in_intervals.txt'

d = pd.read_csv(file, sep='\t')
chrs = d.groupby('chr')

rhos = {}

# start, end, rate
for chr, group in chrs:
	total_length = sum(group.end - group.start)
	total_rho = sum((group.end - group.start) * group.rate)

	rho = total_rho / total_length
	rho = rho * 1e-8 * 3.37e6
	
	rhos[chr] =  '%.3f' % (rho)

string = '{'
for chr in sorted(rhos.keys()):
	string += '\'%s\': %s, ' % (chr, rhos[chr])
string += '}'

print string  
