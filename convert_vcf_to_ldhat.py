import gzip
import re
import numpy as np
import glob

vcfs = glob.glob('/mnt/lustre/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.*filtered.coverage.vqsr.vcf.gz')
out_dir = '/mnt/lustre/home/sonal.singhal1/LTF/phasing/LDhat/'

for vcf in vcfs:
	print vcf
	array = []
	sites = []
	inds = []
	chr = ''
	cur_pos = 0
	if re.search('(chr\S+)\.filtered', vcf):
		chr = re.search('(chr[a-z|A-Z|\d|\__]+)\.filtered', vcf).group(1)
	chr_len = ''
	o = gzip.open(vcf, 'r')

	for line in o:
		if re.search('^#CHROM', line):
			d = re.split('\t', line)
			inds = d[9:len(d)]
		elif re.search('contig=<ID=%s,length=(\d+)' % chr, line):
			chr_len = re.search('contig=<ID=%s,length=(\d+)' % chr, line).group(1)
		elif not re.search("^#", line):
			# this is a possible variant line
			d = re.split("\t", line)
			allele_num = 0
			allele_counts = dict()
			# go through all the alleles (add 1 because of reference!)
			for i in range(len(re.split(',', d[4]))+1):
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
				if int(d[1]) > cur_pos:
					array.append(genos)
					sites.append(int(d[1]))
				cur_pos = int(d[1])
			
	o.close()

	num_sites = len(sites)
	tracker = 1
	genos = np.asarray(array).T.tolist()
	while (3800 * (tracker - 1)) < num_sites:
		out_site = out_dir + 'sites_locs/filtered.coverage.vqsr.%s_chunk%s.locs' % (chr, tracker)

		start = 3800 * (tracker - 1)
		end = start + 4000
		# take care of the boundary condition
		if end > num_sites:
			end = num_sites
		
		chunk_sites = end - start
		
		out = open(out_site, 'w')
		# how length is defined is a bit wonky, but i think this should generally work.
		chr_len = sites[end - 1] / float(1000)
		out.write('%s\t%s\tL\n' % (chunk_sites, chr_len))
		for site in sites[start:end]:
			out.write('%s\t' % (site / float(1000)))
		out.write('\n')
		out.close()
	
		out_geno = out_dir + 'sites_locs/filtered.coverage.vqsr.%s_chunk%s.sites' % (chr, tracker)
		out = open(out_geno, 'w')
		out.write('%s\t%s\t2\n' % (len(inds), chunk_sites)) 
		for id, geno in zip(inds, genos):
			out.write('>geno_%s\n' % id) 
			out.write('%s\n' % ''.join(geno[start:end]))
		out.close()
		tracker += 1
