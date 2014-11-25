import glob
import re

dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/map_simulations/'
seq_size = 1e6
bpen = 5
hotspot_files = glob.glob('%s/hotspot/*.txt' % dir)

out_file = '%sstability_of_uncertainty.csv' % dir
out = open(out_file, 'w')
out.write('rho_strength,simulation_number,bpen,pos,actual_rec,estimated_rec,ratio,within\n')


mean_rho_values = {}
f = open('%smean_rho.txt' % dir, 'r')
for l in f:
	d = re.split('\s+', l.rstrip())
	mean_rho_values[ int(d[0]) ] = float(d[2])
f.close()


for hotspot_file in hotspot_files:

	rho_ix = int(re.search('hotspot(\d+).txt', hotspot_file).group(1))
	rho = mean_rho_values[rho_ix]
		
	rec_rates = {}
	hotspot_f = open(hotspot_file, 'r')
	for l in hotspot_f:
		d = re.split('\s+', l)
		start = int(seq_size * float(d[0]))
		end = int(seq_size * float(d[1]))
		for i in range(start, end):
			rec_rates[i] = float(d[2]) * rho
	hotspot_f.close()

	map_files = glob.glob('%srec_maps/recombination%s_*%s.txt' % (dir, rho_ix, bpen))

	for file in map_files:
		sim_num = re.search('_(\d+)_', file).group(1)

		f = open(file, 'r')
		for i in range(3):
			f.next()
		for l in f:
			d = re.split('\s', l)
			i = int(d[0])

			estimated = float(d[2])
			p025  = float(d[3])
			p975 = float(d[4])

			actual = rho
			if i in rec_rates:
				actual = rec_rates[i]
			within = False
			if actual >= p025:
				if actual <= p975:
					within = True

			out.write('%.3f,%s,%s,%s,%.3f,%.3f,%.3f,%s\n' % (rho, sim_num, bpen, i, 
				actual, float(d[2]), (actual/float(d[2])), within))

		f.close()
out.close()

