from itertools import izip
import gzip
import re
import sys

def get_vcf(vcffile, type):
        file = gzip.open(vcffile, 'r')
        var = dict()
        for l in file:
                if not re.match("#", l):
                        d = re.split('\t', l)

                        alleles = [d[3]] + re.split(",", d[4])
                        
                        # do not want to include indels in this variant dictionary
                        # cannot use these to polarize
                        indel = False
                        for allele in alleles:
                                if len(allele) > 1:
                                        indel = True
                        
                        if not indel:
                                # do not want any alleles included that are at zero frequency
                                # not useful for anything but to show ref, and we already have that
                                allele_counts = dict()
                                for i in range(len(alleles)):
                                        allele_counts[alleles[i]] = len(re.findall('%s\/' % i, l)) + len(re.findall('\/%s' % i, l))
                                tot_n = sum(allele_counts.values())
                                
				af_real = dict()

				for i in allele_counts:
                                        af = allele_counts[i] / float(tot_n)
                                        if af > 0:
                                        	af_real[i] = '%.2f' % af
				# for the ingroup, do not want to include fixed or polyallelic sites, so deal with that here
				if type == 'in':
					if len(af_real) == 2:
                                		var[int(d[1])] = af_real
				else:
					var[int(d[1])] = af_real
        file.close()
        return var

  
vcf_in = '/mnt/lustre/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chr1.filtered.coverage.vqsr.vcf.test.gz'
vcf_out1 = '/mnt/lustre/home/sonal.singhal1/ZF/after_vqsr/by_chr/gatk.ug.zf.chr1.filtered.coverage.vqsr.vcf.test.gz'
vcf_out2 = '/mnt/lustre/home/sonal.singhal1/DBF/after_vqsr/by_chr/gatk.ug.dbf.chr1.filtered.coverage.vqsr.vcf.test.gz'

# importantly, all these genomes are interleaved at 60chr long
# this prevents me from having to read the entire genome into memory
# should be much faster than the access methods too
#
# masked genome files for each species; genome_out1 and genome_out2 are outgroup
genome_out2 = '/mnt/lustre/home/sonal.singhal1/DBF/masked_genome/DBF.masked_genome.fa'
genome_out1 = '/mnt/lustre/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.fa'
# another outgroup, even further
genome_out3 = '/mnt/lustre/home/sonal.singhal1/Darwin/g_magnirostris/Geospiza_magnirostris.fasta'
# the genome that tells the reference site at each position
genome_ref = '/mnt/lustre/home/sonal.singhal1/reference/taeGut1_60.bamorder.fasta'

var_in = get_vcf(vcf_in, 'in')
var_out1 = get_vcf(vcf_out1, 'out')
var_out2 = get_vcf(vcf_out2, 'out')

f_out1 = open(genome_out1, 'r')
f_out2 = open(genome_out2, 'r')
f_out3 = open(genome_out3, 'r')
f_ref = open(genome_ref, 'r')

def get_history(b, pos, var_out, ref):
	if b < 4:
		if pos in var_out:
			return var_out[pos]
		else:
			return ref
	else:
		return 'N'


def call_ancestral_allele(b_in, b_out1, b_out2, b_out3, b_ref):
	# oh no! there is variation in the outgroup
	if isinstance(b_out1, dict) or isinstance(b_out2, dict):
		if b_out3 == b_ref:
			if b_out3 in b_in:
				return b_ref
	# good, there is no variation in the outgroup
	else:
		if b_out1 != 'N' and b_out2 != 'N':
			if b_out1 == b_out2:
				if b_out1 in b_in:
					# this is the winner. we had good coverage at both outgroups, and they
					#	are the same fixed allele
					return b_out1
		else:
			# trust it if: there are at least two non-N bases, they agree, and they are one of the two variants
			# trust it if: the finch genome base agrees with the reference base, and they are one of the two variants
			bases = []
			for b in [b_out1, b_out2, b_out3]:
				if b != 'N':
					bases.append(b)
			if len(bases) > 1:
				if bases[0] == bases[1]:
					if bases[0] in b_in:
						# reference base is the winner
						return bases[0]
			if b_out3 == b_ref:
				if b_out3 in b_in:
					return b_out3	
	if b_out3 in b_in:
		return b_out3
	return 'N'

line_count = 0
for l_out1, l_out2, l_out3, l_ref in izip(f_out1, f_out2, f_out3, f_ref):
	l_out1 = l_out1.strip()
	l_out2 = l_out2.strip()
	l_out3 = l_out3.strip()
	l_ref = l_ref.strip()

	if re.match(">(\S+)", l_out1):
		chr = re.match(">(\S+)", l_out1).group(1)
		# vcf = dir + 'gatk.ug.ltf.%s.filtered.coverage.biallelicSNPs.vqsr.vcf.gz' % chr
		
		# need to empty out the variation dictionary
		# var = dict()
		#  need to reset the line count
		line_count = 0
	
		# vcfopen = gzip.open(vcf, 'r')
		# for l in vcfopen:
		#	if not re.match("#", l):
		#		d = re.split('\t', l)
		#		var[int(d[1])] = {'ref': d[3], 'alt': d[4]}
		# vcfopen.close()
	else:
		for ix, (b_out1, b_out2, b_out3, b_ref) in enumerate(izip(l_out1, l_out2, l_out3, l_ref)):	
			b_out1 = int(b_out1)
			b_out2 = int(b_out2)
			b_out3 = b_out3.upper()

			# get the current genome position	
			# need to add 1 because python indexes at 0; genomes are at index 1
			pos = line_count * 60 + ix + 1
	
			# this is a variable position
			if pos in var_in:
				b_in = var_in[pos]
				
				b_out1call = get_history(b_out1, pos, var_out1, b_ref)
				b_out2call = get_history(b_out2, pos, var_out2, b_ref)
				
				anc_allele = call_ancestral_allele(b_in, b_out1call, b_out2call, b_out3, b_ref)
	
				if anc_allele != 'N':
					derived_freq = 1 - float(b_in[anc_allele])
				else:
					print b_in
					print b_out1call
					print b_out2call
					print b_out3
					print b_ref
					print '*****'


			if pos > 500000:
				sys.exit()

		line_count += 1

f_in.close()
f_out1.close()
f_out2.close()
f_out3.close()
f_ref.close()
