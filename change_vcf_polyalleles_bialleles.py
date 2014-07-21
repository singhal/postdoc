import re
import gzip
import glob

vcfs = glob.glob('/mnt/lustre/home/sonal.singhal1/ZF/after_vqsr/by_chr/*coverage.vqsr*')

for vcf in vcfs:
	infile = gzip.open(vcf, 'r')
	vcfout = vcf.replace('coverage', 'coverage.recoded_biallelic')
	outfile = gzip.open(vcfout, 'w')

	for l in infile:
		if re.search('^#', l):
			outfile.write(l)
		else:
			d = re.split('\t', l)
			if len(d[3]) == 1:
				if len(d[4]) == 1:
					outfile.write(l)
				else:
					if re.search('\,', d[4]):
						alleles = re.split('\,', d[4])
						indel = False
						for allele in alleles:
							if len(allele) > 1:
								indel = True

						if not indel:
							for ix, allele in enumerate(alleles, start=1):
								tmp = l
								tmp = tmp.rstrip()
								for i in range(1,len(alleles)+1):
									if i == ix:
										tmp = tmp.replace('%s/' % i, '%s/' % 1)
										tmp = tmp.replace('/%s' % i, '/%s' % 1)
									else:
										tmp = tmp.replace('%s/' % i, '%s/' % 0)
                                                                        	tmp = tmp.replace('/%s' % i, '/%s' % 0)
								d2 = re.split('\t', tmp)
								outfile.write('\t'.join(d2[0:4]) + '\t' + allele + '\t' + '\t'.join(d2[5:len(d2)]) + '\n')
	infile.close()
	outfile.close()
				 	
