from itertools import izip
import gzip
import re
import sys
import subprocess
import os
import argparse

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
                                        allele_counts[alleles[i]] = len(re.findall('\s%s[\/|\:]' % i, l)) + len(re.findall('\/%s' % i, l))
                                tot_n = sum(allele_counts.values())
                                
				if tot_n > 0:

					af_real = dict()
	
					for i in allele_counts:

                	                        af = allele_counts[i] / float(tot_n)
                        	                if af > 0:
                                	        	af_real[i] = '%.3f' % af
					# for the ingroup, do not want to include fixed or polyallelic sites, so deal with that here
					if type == 'in':
						if len(af_real) == 2:
                                			var[int(d[1])] = af_real
					else:
						var[int(d[1])] = af_real
        file.close()
        return var


def get_history(b, pos, var_out, ref):
	if b < 4:
		if pos in var_out:
			if len(var_out[pos]) == 1:
				return var_out[pos].keys()[0]	
			else:
				return var_out[pos]
		else:
			return ref
	else:
		return 'N'


def check_equal(iterator):
	return len(set(iterator)) <= 1


def call_ancestral_allele(b_in, b_out1, b_out2, b_out3, b_out4):
	# oh no! there is variation in the outgroup
	if not isinstance(b_out1, dict) and not isinstance(b_out2, dict):
		if b_out1 != 'N' and b_out2 != 'N':
			if b_out1 == b_out2:
				if b_out1 in b_in:
					# this is the winner. we had good coverage at both outgroups, and they
					#	are the same fixed allele
					return b_out1
	far_out = []
	for bp in (b_out3, b_out4):
		if bp in ['A', 'T', 'C', 'G']:
			far_out.append(bp)
	if len(far_out) > 0:
		if check_equal(far_out):
			if far_out[0] in b_in:
				return far_out[0]	
	return 'N'


def get_chromosome(genome, chr):
	outfile = genome + '_' + chr
	subprocess.call('~/bin/samtools-0.1.19/samtools faidx %s %s > %s' % (genome, chr, outfile), shell=True)
	out_f = open(outfile, 'r')
	locus_name = out_f.next()
	return outfile, out_f


def pp_hist(base_info):
	if isinstance(base_info, dict):
		return "|".join("%s:%s" % bp_pair for bp_pair in base_info.items())
	else:
		return base_info

def trawl_genome(out_file, chr, f_out1, f_out2, f_out3, f_out4, f_ref, var_in, var_out1, var_out2):
	out = open(out_file, 'w')
	out.write('chr,position,ingroup,outgroup1,outgroup2,faroutgroup3,faroutgroup4,ancestralallele\n')
	for line_count, (l_out1, l_out2, l_out3, l_out4, l_ref) in enumerate(izip(f_out1, f_out2, f_out3, f_out4, f_ref)):
		l_out1 = list(l_out1.rstrip())
		l_out2 = list(l_out2.rstrip())
		l_out3 = list(l_out3.upper().rstrip())
		l_out4 = list(l_out4.upper().rstrip())
		l_ref =  list(l_ref.upper().rstrip())

		l_out1 = [int(i) for i in l_out1]
		l_out2 = [int(i) for i in l_out2]

		for ix, (b_out1, b_out2, b_out3, b_out4, b_ref) in enumerate(izip(l_out1, l_out2, l_out3, l_out4, l_ref)):
			# get the current genome position	
			# need to add 1 because python indexes at 0; genomes are at index 1
			pos = line_count * 60 + ix + 1
	
			# this is a variable position
			if pos in var_in:
				b_in = var_in[pos]
				
				b_out1call = get_history(b_out1, pos, var_out1, b_ref)
				b_out2call = get_history(b_out2, pos, var_out2, b_ref)
				
				anc_allele = call_ancestral_allele(b_in, b_out1call, b_out2call, b_out3, b_out4)

				out.write('%s,%s,%s,%s,%s,%s,%s,%s\n' % \
					(chr, pos, pp_hist(b_in), pp_hist(b_out1call), pp_hist(b_out2call), b_out3, b_out4, anc_allele))
	out.close()
	return


def clean_up(chr_out, f_out):
	f_out.close()
	os.remove(chr_out)
	

def main():
	parser = argparse.ArgumentParser()
        parser.add_argument("--chr", help="chromosome for which to run analysis")
        args = parser.parse_args()
        chr = args.chr

	# LTF
	# vcf_in = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/for_shapeit/gatk.ug.ltf.%s.allfilters.recoded_biallelicSNPs.vcf.gz' % chr
	# if chr == 'chrZ':
	#	vcf_in = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/for_shapeit/gatk.ug.ltf.chrZ.allfilters.recodedsex.recoded_biallelicSNPs.vcf.gz'
	# vcf_out1 = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/gatk.ug.all_zf.%s.coverage.repeatmasked.filtered.vqsr2.vcf.gz' % chr
	# out_file = '/mnt/gluster/home/sonal.singhal1/LTF/ancestral_allele/ancestral_allele.%s.csv' % chr
	# all genomes must be at 60chr per line!!!
	# genome_out1 = '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.repeat_masked.switch_masked.fa'

	# ZF
	vcf_in = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/for_shapeit/gatk.ug.finch19.%s.allfilters.recoded_biallelicSNPs.nomissing.vcf.gz' % chr
        if chr == 'chrZ':
                vcf_in = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/for_shapeit/gatk.ug.finch19.chrZ.allfilters.recodedsex.recoded_biallelicSNPs.vcf.gz'
	vcf_out1 = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.coverage.repeatmasked.filtered.vqsr2.vcf.gz' % chr
        out_file = '/mnt/gluster/home/sonal.singhal1/ZF/ancestral_allele/ancestral_allele.%s.csv' % chr
        genome_out1 = '/mnt/gluster/home/sonal.singhal1/LTF/masked_genome/LTF.masked_genome.repeat_masked.fa'
       
	# same for both species
	vcf_out2 = '/mnt/gluster/home/sonal.singhal1/DBF/after_vqsr/by_chr/gatk.ug.dbf.%s.filtered.coverage.vqsr.vcf.gz' % chr
	genome_out2 = '/mnt/gluster/home/sonal.singhal1/DBF/masked_genome/DBF.masked_genome.fa'
	# other outgroups, even further
	genome_out3 = '/mnt/gluster/home/sonal.singhal1/Darwin/g_magnirostris/Geospiza_magnirostris.fasta'
	genome_out4 = '/mnt/gluster/home/sonal.singhal1/Darwin/g_fortis/Geospiza_fortis.fasta'
	genome_ref = '/mnt/gluster/home/sonal.singhal1/reference/taeGut1_60.bamorder.fasta'
	
	var_in = get_vcf(vcf_in, 'in')
	var_out1 = get_vcf(vcf_out1, 'out')
	var_out2 = get_vcf(vcf_out2, 'out')

	chr_out1, f_out1 = get_chromosome(genome_out1, chr)
	chr_out2, f_out2 = get_chromosome(genome_out2, chr)
	chr_out3, f_out3 = get_chromosome(genome_out3, chr)
	chr_out4, f_out4 = get_chromosome(genome_out4, chr)
	chr_ref, f_ref = get_chromosome(genome_ref, chr)

	trawl_genome(out_file, chr, f_out1, f_out2, f_out3, f_out4, f_ref, var_in, var_out1, var_out2)
	
	clean_up(chr_out1, f_out1)
	clean_up(chr_out2, f_out2)
	clean_up(chr_out3, f_out3)
	clean_up(chr_out4, f_out4)
	clean_up(chr_ref, f_ref)

if __name__ == "__main__":
    main()
