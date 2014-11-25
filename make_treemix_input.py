import re
import glob
import gzip

chrs = ['chr1', 'chr1A', 'chr1B', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11',
	'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24',
	'chr25', 'chr26', 'chr27', 'chr28', 'chrLG2', 'chrLG5', 'chrLGE22']

out_file = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/treemix/treemix_auts_inclDBF.txt'
out_f = open(out_file, 'w')
out_f.write('pop1 pop2 pop3 pop4\n')

def parse_vcf(var, vcf, popname):
	f = gzip.open(vcf, 'r')
	for l in f:
		if not re.search('^#', l):
			d = re.split('\t', l)
			pos = int(d[1])
			# easiest way to get biallelic SNPs
			if len(d[3]) == 1 and len(d[4]) == 1:
				if pos not in var:
					var[pos] = {'ref': d[3], 'alt': d[4]}
				if var[pos]['alt'] == d[4] and var[pos]['ref'] == d[3]:
					ref = 0
					alt = 0
					for geno in d[9:]:
						geno = re.search('^([^:]+)', geno).group(1)
						geno = re.split('[\||/]', geno)
						for allele in geno:
							if allele == '0':
								ref += 1
							else:
								alt += 1	
					var[pos][popname] = [ref, alt]
	f.close()
	return var


def parse_vcf2(var, vcf, popname1, num1, popname2, num2):
        f = gzip.open(vcf, 'r')
        for l in f:
                if not re.search('^#', l):
                        d = re.split('\t', l)
                        pos = int(d[1])
			if len(d[3]) == 1 and len(d[4]) == 1:
                        	if pos not in var:
                                	var[pos] = {'ref': d[3], 'alt': d[4]}
                        	if var[pos]['alt'] == d[4] and var[pos]['ref'] == d[3]:
                        		ref = 0
                                        alt = 0
                                        for geno in d[9:9+num1]:
                                                geno = re.search('^([^:]+)', geno).group(1)
                                                geno = re.split('[\||/]', geno)
                                                for allele in geno:
                                                        if allele == '0':
                                                                ref += 1
                                                        else:
                                                                alt += 1
                                        var[pos][popname1] = [ref, alt]

					ref = 0 
                                        alt = 0 
                                        for geno in d[9+num1:9+num1+num2]:
                                                geno = re.search('^([^:]+)', geno).group(1)
                                                geno = re.split('[\||/]', geno)
                                                for allele in geno:
                                                        if allele == '0':
                                                                ref += 1
                                                        else:
                                                                alt += 1
                                        var[pos][popname2] = [ref, alt]
        f.close()
        return var


for chr in chrs:
	zf_vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.phased.vcf.gz' % chr
	ltf_vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.filtered.coverage.repeatmasked.vqsr.phased.vcf.gz' % chr
	dbf_vcf = '/mnt/gluster/home/sonal.singhal1/DBF/after_vqsr/by_chr/gatk.ug.dbf.%s.filtered.coverage.vqsr.vcf.gz' % chr

	var = {}
	var = parse_vcf(var, zf_vcf, 'pop1')
	var = parse_vcf2(var, ltf_vcf, 'pop2', 10, 'pop3', 10)
	var = parse_vcf(var, dbf_vcf, 'pop4')

	for pos in sorted(var.keys()):
		pop_line = ''
		for pop, count in [('pop1', 38), ('pop2', 20), ('pop3', 20), ('pop4', 2)]:
			if pop in var[pos]:
				pop_line += '%s,%s ' % (var[pos][pop][0], var[pos][pop][1])
			else:
				pop_line += '%s,%s ' % (count, 0)	
		out_f.write(pop_line + '\n')
out_f.close()
