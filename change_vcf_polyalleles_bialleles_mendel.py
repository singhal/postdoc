import re
import gzip
import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
args = parser.parse_args()
chr = args.chr

vcffile = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/gatk.ug.all_zf.%s.coverage.filtered.repeatmasked.vqsr.vcf.gz' % chr

for vcf in [vcffile]:
	infile = gzip.open(vcf, 'r')
	vcfout = vcf.replace('repeatmasked.vqsr', 'repeatmasked.recoded_biallelicSNPs.nomendel')
	outfile = gzip.open(vcfout, 'w')

	errorfile = '/mnt/gluster/home/sonal.singhal1/ZF/mendelian_errors/plink_results/all_zf.me.%s.mendel' % chr
        sites_file = '/mnt/gluster/home/sonal.singhal1/ZF/mendelian_errors/all_zf.%s.map' % chr

        sites_f = open(sites_file, 'r')
        sites = {}
        for l in sites_f:
                d = re.split('\s+', l)
                sites[d[1]] = int(d[3])
        sites_f.close()
        
        errors = {}
        error_f = open(errorfile, 'r')
        junk = error_f.next()
        for l in error_f:
                d = re.split('\s+', l)
                errors[sites[d[4]]] = 1
        error_f.close()

	for l in infile:
		if re.search('^#', l):
			if re.search('^#CHROM', l):
				d = re.split('\t', l.rstrip())
				outfile.write(l)
			else:
				outfile.write(l)
		else:
			d = re.split('\t', l.rstrip())
			if re.search('PASS', l):
				indel = False
				alleles = [d[3]] + re.split(",", d[4])
				for allele in alleles:
					if len(allele) > 1:
						indel = True
	
				if not indel:
					actual_alleles = []
                	                for ix, allele in enumerate(alleles):
                	                        count = len(re.findall('%s\/' % ix, l)) + len(re.findall('\/%s' % ix, l))
						if count > 0:
							actual_alleles.append(allele)
	
					if len(actual_alleles) == 2:
						if int(d[1]) not in errors:
							if len(d[4]) == 1:
								outfile.write(l)
							else:
								tmp = l.rstrip()
								for ix, allele in enumerate(alleles):
									if allele in actual_alleles:
										recode = actual_alleles.index(allele)
										tmp = tmp.replace('%s/' % ix, '%s/' % recode)
										tmp = tmp.replace('/%s' % ix, '/%s' % recode)	
								d2 = re.split('\t', tmp)
	
								outfile.write(  '\t'.join(d2[0:3]) + '\t' + actual_alleles[0] + 
									 	'\t' + actual_alleles[1] + '\t.\t' + d2[6]
										+ '\t.\t' +  '\t'.join(d2[8:]) + '\n')
	infile.close()
	outfile.close()
				 	
