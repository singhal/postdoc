import gzip
import re
import numpy as np
import glob
import pandas as pd

hot_file = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/seqldhot_hotspots/spot2kb_flank40kb.seqldhot_validate_hotspots.csv'
out_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/hotspot_checks/LDhat/'
vcf_base = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.vcf.gz'

d = pd.read_csv(hot_file)
d = d[d.species == 'ZF']
d = d[d.zmatch_heat > 10]
d = d[d.zmatch_lk > 10]
# will take more work before I can check chromosome Z because of the different number of chromosomes
d = d[d.chr != 'chrZ']

chrs = d.groupby('chr')
for chr, groupchr in chrs:
	vcf = vcf_base % chr
	o = gzip.open(vcf, 'r')
	var = {}
	for line in o:
		if re.search('^#CHROM', line):
			d = re.split('\t', line)
			inds = d[9:len(d)]
		elif not re.search("^#", line):
			# this is a possible variant line
			d = re.split("\t", line)

			# don't wantt indels
			indels = False
			alleles = [d[3]] + re.split(',', d[4])
			for allele in alleles:
				if len(allele) > 1:
					indels = True

			if not indels:
				allele_num = 0
				allele_counts = dict()
				for i, allele in enumerate(alleles):
					# figure out how many times each allele was genotyped
					allele_counts[i] = len(re.findall(str(i) + '/', line)) + len(re.findall('/' + str(i), line))
					# figure out how many alleles are actually segregating in the population
					if allele_counts[i] > 0:
						allele_num += 1
					else:
						allele_counts.pop(i, None)

				# if there are exactly two alleles, then yay!
				#	-LDhat cannot use fixed or multiallelic sites
				if allele_num == 2:
					# no singletons to make it comparable
					if min(allele_counts.values()) > 1:
						geno_ids = {}
						geno_ids['./.'] = '?'
						alleles = sorted(allele_counts.keys())
						geno_ids[str(alleles[0]) + '/' + str(alleles[0])] = '0'
						geno_ids[str(alleles[1]) + '/' + str(alleles[1])] = '1'
						geno_ids[str(alleles[0]) + '/' + str(alleles[1])] = '2'

						genos = []
						# should write this as a map function to speed it up
						for ix, ind in enumerate(d[9:len(d)]):
							geno = re.search('(\S\/\S)', ind).group(1)
							genos.append(geno_ids[geno])
						var[int(d[1])] = genos
	o.close()

	for center in groupchr.spot_start:
		start = center - 50000
		if start < 1:
			start = 1
		end = center + 50000

		sites = filter(lambda x: x >= start, var.keys())
		sites = filter(lambda x: x <= end, sites)
		sites = sorted(sites) 

		out_site = out_dir + 'sites_locs/%s_%s.locs' % (chr, center)

		out = open(out_site, 'w')
		# how length is defined is a bit wonky, but i think this should generally work.
		chr_len = max(sites) / float(1000)
		out.write('%s\t%s\tL\n' % (len(sites), chr_len))
		for site in sites:
			out.write('%s\t' % (site / float(1000)))
		out.write('\n')
		out.close()
	
		out_geno = out_dir + 'sites_locs/%s_%s.sites' % (chr, center)
		out = open(out_geno, 'w')
		out.write('%s\t%s\t2\n' % (len(inds), len(sites)))
		for ind_ix, ind in enumerate(inds):
			geno = ''
			for pos in sites:
				geno += str(var[pos][ind_ix])
			out.write('>geno_%s\n' % ind) 
			out.write('%s\n' % geno)
		out.close()
