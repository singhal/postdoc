import re
import glob
import gzip

unrel_vcf = '/mnt/lustre/home/sonal.singhal1/ZF/after_vqsr/by_chr/gatk.ug.zf.chr22.filtered.coverage.vqsr.vcf.gz'
rel_vcf = '/mnt/lustre/home/sonal.singhal1/ZF/phasing/family_approach/family_vcf/gatk.ug.zf_fam.chr22.variable.coverage.vqsr.vcf.gz'

unrelvar = dict()
unrelfile = gzip.open(unrel_vcf, 'r')
for l in unrelfile:
	if not re.search('^#', l):
		d = re.split('\t', l)

                allele_counts = dict()
                alleles = [d[3]] + re.split(',', d[4])
                # go through all the alleles (add 1 because of reference!)
                for ix, allele in enumerate(alleles):
                        # figure out how many times each allele was genotyped
                        allele_counts[allele] = len(re.findall(str(ix) + '/', l)) + len(re.findall('/' + str(ix), l))
                        # figure out how many alleles are actually segregating in the population
                        if allele_counts[allele] == 0:
                                allele_counts.pop(allele, None)

                # if there are exactly two alleles, then yay!
                #       -ShapeIt cannot use fixed or multiallelic sites
                if len(allele_counts.keys()) == 2:
			geno = []
			for match in re.finditer('(\S\/\S)', l):
				geno.append(match.group(1))
			
			for ix, allele in enumerate(alleles):
				geno = [g.replace(str(ix) + '/', '%s ' % allele) for g in geno]
				geno = [g.replace(str(ix), allele) for g in geno]
			geno = [g.replace('./.', '0 0') for g in geno]

			unrelvar[d[1]] = {'all': alleles, 'used': allele_counts.keys(), 'geno': geno}
unrelfile.close()

relfile = gzip.open(rel_vcf, 'r')
for l in relfile:
	if not re.search('^#', l):
		d = re.split('\t', l)

		# only care about positions that are variable in the family
		if d[1] in unrelvar:
			allele_counts = dict()
			alleles = [d[3]] + re.split(',', d[4])
			# go through all the alleles (add 1 because of reference!)
			for ix, allele in enumerate(alleles):
				# figure out how many times each allele was genotyped
				allele_counts[allele] = len(re.findall(str(ix) + '/', l)) + len(re.findall('/' + str(ix), l))
				# figure out how many alleles are actually segregating in the population
				if allele_counts[allele] == 0:
					allele_counts.pop(allele, None)

		# if there are exactly two alleles, then yay!
		#       -ShapeIt cannot use fixed or multiallelic sites
		if len(allele_counts.keys()) == 2:

			geno = []
                        for match in re.finditer('(\S\/\S)', l):
                                geno.append(match.group(1))

			indel = False
			for ix, allele in enumerate(alleles):
				geno = [g.replace(str(ix) + '/', '%s ' % allele) for g in geno]
				geno = [g.replace(str(ix), allele) for g in geno]
				if len(allele) > 1:
					indel = True

			if d[1] in unrelvar:
				tot_alleles = list( set(allele_counts.keys()) | set( relvar[d[1]]['used']) )	
		
				if len(tot_alleles) == 2:
					# this is the winning case
					geno = geno + relvar[d[1]]['geno']
				else:
					# can still keep variation, but missing
					geno = geno + ['0 0'] * 5
				
			else:
				# still can keep the variation, but missing
				geno = geno + ['0 0'] * 5

			
			print geno
"""
