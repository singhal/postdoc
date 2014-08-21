import re
import copy
import glob
import gzip
import numpy as np
import os.path
import sys

pedfile = '/mnt/gluster/home/sonal.singhal1/ZF/zf_inds_plink.ped'
dir = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/'

chr_dict = {	'chr1': 1, 'chr2': 2, 'chr3': 3, 'chr4': 4, 'chr5': 5, 'chr6': 6, 'chr7': 7, 'chr8': 8, 'chr9': 9, 'chr10': 10, 'chr11': 11, 'chr12': 12,
		'chr13': 13, 'chr14': 14, 'chr15': 15, 'chr16': 16, 'chr17': 17, 'chr18': 18, 'chr19': 19, 'chr20': 20, 'chr21': 21, 'chr22': 22,
		'chr23': 33, 'chr24': 34, 'chr25': 35, 'chr26': 36, 'chr27': 37, 'chr28': 38, 'chr1A': 30, 'chr1B': 31, 'chr4A': 32, 'chrLG2': 39,
		'chrLG5': 40, 'chrLGE22': 41, 'chrZ': 'X' } 

chrs = [ 'chrLGE22', 'chr1A', 'chr1B', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
         'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
         'chrLG2', 'chrLG5', 'chr1', 'chrZ' ]

chrs_long = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr1A']

for chr in [i for i in chrs if i not in chrs_long]:
	merge_vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/gatk.ug.all_zf.%s.coverage.filtered.vqsr.vcf.gz' % chr
	chr_num = str(chr_dict[chr])

	ped_out = dir + 'mendelian_errors/all_zf.%s.ped' % chr
	sites_out = dir + 'mendelian_errors/all_zf.%s.map' % chr

	if not os.path.isfile(sites_out):
		ped = dict()
		inds = list()
		sites = dict()
		genos = list()
		file = gzip.open(merge_vcf, 'r')
	
		ped_open = open(pedfile, 'r')
		for l in ped_open:
			d = re.split('\s+', l)
			ped[d[1]] = d
		ped_open.close()
	
		for l in file:
			if re.search('^#CHROM', l):
				l = l.rstrip()
				d = re.split('\t', l)
				inds = d[9:len(d)]
			if re.search('PASS', l):
				d = re.split('\t', l)
	
				indel = False
				
				allele_counts = dict()
				alleles = [d[3]] + re.split(',', d[4])
				# go through all the alleles (add 1 because of reference!)
				for ix, allele in enumerate(alleles):
					# figure out how many times each allele was genotyped
					allele_counts[allele] = len(re.findall(str(ix) + '/', l)) + len(re.findall('/' + str(ix), l))
					# figure out how many alleles are actually segregating in the population
					if allele_counts[allele] == 0:
						allele_counts.pop(allele, None)
					if len(allele) > 1:
						indel = True
		
				# if there are exactly two alleles, then yay!
				#       -ShapeIt cannot use fixed or multiallelic sites
				if len(allele_counts.keys()) == 2:
					geno = []
					for match in re.finditer('(\S\/\S)', l):
						geno.append(match.group(1))
					if indel:
						tracker = 8
						for ix, allele in enumerate(alleles):
							if allele in allele_counts.keys():
								geno = [g.replace(str(ix) + '/', '%s ' % tracker) for g in geno]
								geno = [g.replace(str(ix), '%s' % tracker) for g in geno]
								tracker += 1
					else:
						for ix, allele in enumerate(alleles):				
							geno = [g.replace(str(ix) + '/', '%s ' %allele) for g in geno]
							geno = [g.replace(str(ix), '%s' % allele) for g in geno]
					geno = [g.replace('./.', '0 0') for g in geno]


					missing = False
					if len(filter(lambda x:'0 0' in x, geno)) == len(geno):
						missing = True
		
					# this logic is sufficient because want to keep SNPs bc LDhelmet can actually use them
					# SNPs always appear before indels in files that have both
					if int(d[1]) not in sites and not missing:
						genos.append(geno)
						sites[int(d[1])] = 1
		file.close()

		p_out = open(ped_out, 'w')
                for ind_ix, ind in enumerate(inds):
                        p_out.write('\t'.join(ped[ind]) + '\t')
                        for snp in range(0, len(genos)):
                                p_out.write(genos[snp][ind_ix] + '\t')
                        p_out.write('\n')
                p_out.close()

		s_out = open(sites_out, 'w')
		for ix, site in enumerate(sorted(sites.keys())):
			s_out.write('%s\t%s\t0\t%s\n' % (chr_num, 'SNP' + str(ix), site))
		s_out.close()
