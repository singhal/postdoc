import pandas as pd
import re
import numpy as np

putative_hotspots = '/mnt/gluster/home/sonal.singhal1/simulations/shared/seqldhot/ZF_LTF.putative_hotspots.lambda5.csv'

block = 50e3

d = pd.read_csv(putative_hotspots)
hotspots = {}
for file, mid in zip(d.file, d.midpoint):
	if file not in hotspots:
		hotspots[file] = []
	if len(hotspots[file]) > 0:
		min_dist = np.min([abs(x - mid) for x in hotspots[file]])
		if min_dist > 5000:
			hotspots[file].append(mid)
	else:
		hotspots[file].append(mid)

rhos = {'0.001': '0.0005', '0.002': '0.001', '0.01': '0.005', '0.1': '0.05', '0.5': '0.25', '0.8': '0.4'}

def make_seqldhot_file(sp, file, mids, theta, sh):
	hap_file = '/mnt/gluster/home/sonal.singhal1/simulations/shared/%s/haplo/haplo_%s.fa' % (sp, file)
	rho_file = '/mnt/gluster/home/sonal.singhal1/simulations/shared/%s/maps/recombination_%s_5.txt' % (sp, file)

	seq = {}
        f = open(hap_file, 'r')
        for l in f:
		if re.search('>', l):
			id = re.search('>(\S+)', l).group(1)
			seq[id] = ''
		else:
			seq[id] += l.rstrip()
	f.close()

	for id in seq:
		seq[id] = list(seq[id])

	inds = seq.keys()
	haplo = {}
	for ind in inds:
		haplo[ind] = {}
	for i in range(0,1000000):
		snps = [seq[ind][i] for ind in inds]
		# want only varaible sites
		if len(set(snps)) == 2:
			uniq_snps = list(set(snps))
			# want only snps that aren't singletons
			if snps.count(uniq_snps[0]) > 1 and snps.count(uniq_snps[1]) > 1:
				snp_bin = {}
				for ix, snp in enumerate(uniq_snps):
					snp_bin[snp] = str(ix)
				for ind, snp in zip(inds, snps):
					haplo[ind][i] = snp_bin[snp] 
        del seq


	rho_d = pd.read_csv(rho_file, sep=" ", skiprows=3, header=None, 
                		names=['left_snp', 'right_snp', 'meanrho', 'p0.025', 'p0.975']) 

	for mid in mids:
		out = '/mnt/gluster/home/sonal.singhal1/simulations/shared/seqldhot/%s/putative_hotspot_%s_%s.seqLDhot.txt' % (sp, file, mid)

		seq_start = mid - ( 50e3 / 2.0 )
                if seq_start < 1:
                        seq_start = 1
                seq_end = mid + ( 50e3 / 2.0 )

		sites = filter(lambda x: x >= seq_start, haplo[inds[0]].keys())
                sites = filter(lambda x: x <= seq_end, sites)

		back_rho = rho_d[ rho_d.left_snp >= seq_start ]
                back_rho = back_rho[ back_rho.right_snp  <= seq_end ].meanrho.tolist()
                back_rho = filter(lambda x: np.isfinite(x), back_rho)
                back_rho = np.mean(back_rho)
                if np.isfinite(back_rho):
                        if back_rho == 0:
                                back_rho = 0.0001
                else:
                        back_rho = 0.0001

		sorted_sites = sorted(sites)
                # next, let's make the haplotypes
                tmp_haplo = {}
                for ix in sorted(haplo.keys()):
                        hap = ''
                        for pos in sorted_sites:
                                hap += haplo[ix][pos]
                        if hap not in tmp_haplo:
                                tmp_haplo[hap] = 0
                        tmp_haplo[hap] += 1

                # next, let's make the sequenceLDhot
                o = open(out, 'w')
                o.write('Distinct = %s\nGenes = %s\nLoci = %s\n' % (len(tmp_haplo), len(haplo), len(sites)))
                o.write('I=1\nK = -2\nPositions of loci:\n')
                o.write(' '.join([str(x - int(seq_start) + 1) for x in sorted_sites]) + '\n')
                o.write('Haplotypes\n')
                for hap, count in tmp_haplo.items():
                        o.write('\t%s %s\n' % (hap, count))                     
                o.write('#')
                o.close()

		out_in = '/mnt/gluster/home/sonal.singhal1/simulations/shared/seqldhot/%s/putative_hotspot_%s_%s.seqLDhot.in' % (sp, file, mid)
                o = open(out_in, 'w')
                o.write('Number of runs = 5000\nMIN number of iterations per hotspots = 100\ndriving values (for rho) = 2\n')
                o.write('background rho = %.3f\ntheta (per site) = %s\n' % (back_rho * 1000, theta))
                o.write('abs grid for hotspot likelihood\n0.5 40\nrel grid for hotspots likelihood\n5 100\n')
                o.write(' sub-region (number of SNPS; length (bps); frequency (bps))\n')
                o.write(' 10 2000 1000\n#\n')
                o.close()
                sh.write('/mnt/lustre/home/sonal.singhal1/bin/sequenceLDhot/sequenceLDhot %s %s\n' % (out_in, out))		


sh = open('/mnt/gluster/home/sonal.singhal1/simulations/shared/seqldhot/ZF_LTF.run_seqldhot.sh', 'w')
for file in hotspots:
	zf_file = file
	zf_file = zf_file.replace('recombination_', '')
	zf_file = re.sub('_5$', '', zf_file)
	zf_rho = re.search('rho([\d|\-|\.]+)', file).group(1)
	ltf_rho = rhos[zf_rho]
	ltf_file = zf_file.replace(zf_rho, ltf_rho)

	make_seqldhot_file('ZF', zf_file, hotspots[file], 0.0064, sh)
	make_seqldhot_file('LTF', ltf_file, hotspots[file], 0.0045, sh)
sh.close()
