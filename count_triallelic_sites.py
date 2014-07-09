import re

# meant to count the number of sites with more than two alleles
# considers how many of these are actually multi-allelic -- sometimes, no
#	individual has a reference allele, making it actually biallelic

# filtered VCF with only passing sites
file = '/ifs/data/c2b2/mp_lab/ss4776/longtailedfinch/after_vqsr/gatk.ug.ltf.allchrs.snps.indels.vqsr2.filtered.vcf'

# summary outfile
outfile = '/ifs/data/c2b2/mp_lab/ss4776/longtailedfinch/multiallelic_sites.txt'

# the type of sites
sites = {'all_variable': 0, 'alleles>2_reported': 0, 'alleles>2_actual': 0, 'alleles>2_fake': 0}

count = 0 
r = open(file, 'r')
for l in r:
	# don't want header stuff
	if not re.search('^#', l):

		d = re.split("\t", l)
		sites['all_variable'] += 1

		if re.search(',', d[4]):
			sites['alleles>2_reported'] += 1
			
			# initialize the dictionary
			# each allele gets its own number
			allele_counts =  dict()
			allele_counts['0'] = 0
			for i in range(1,len(re.split(",", d[4]))+1):
				allele_counts[str(i)] = 0
			
			# count how many of these alleles got genotyped
			for allele in allele_counts:
				allele_counts[allele] = len(re.findall(allele + "/", l)) + \
							len(re.findall("/" + allele, l))
			
			# how many of these alleles actually segregate in the population?
			num_alleles = 0
			for allele in allele_counts:
				if allele_counts[allele] > 0:
					num_alleles += 1			

			if num_alleles > 2:
				count += 1
				sites['alleles>2_actual'] += 1
			else:
				sites['alleles>2_fake'] += 1
r.close()			

out = open(outfile, 'w')
for site_type in sites:
	out.write('%s\t%s\n' % (site_type, sites[site_type]))
out.close()
