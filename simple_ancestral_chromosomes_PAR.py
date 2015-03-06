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
				genos = []
                                for geno in d[9:]:
                                        geno = re.search('^([^:]+)', geno).group(1)
                                        genos += re.split('[|/]', geno)

				af = {}
				for ix, allele in enumerate(alleles):
                                        freq = genos.count(str(ix)) / float(len(genos))
                                        if freq > 0:
                                        	af[allele] = '%.3f' % freq
				var[int(d[1])] = af
        file.close()
        return var


def get_history(pos, var, b_ref):
	if pos in var:
		if len(var[pos]) == 1:
			return var[pos].keys()[0]	
		else:
			return var[pos]
	else:
		return b_ref


def check_equal(iterator):
	return len(set(iterator)) <= 1


def most_common(lst):
    return max(set(lst), key=lst.count)


# need to work on this
def call_ancestral_allele(b1, b2, b3, b4, b5):
	fixed_alleles = []
	for b in [b1, b2, b3]:
		if not isinstance(b, dict):
			fixed_alleles.append(b)

	if len(fixed_alleles) > 1:
		if check_equal(fixed_alleles):
			return fixed_alleles[0]

	for b in [b4, b5]:
		if b in ['A', 'T', 'C', 'G']:
			fixed_alleles.append(b)

	if len(fixed_alleles) > 0:
		putative_anc = most_common(fixed_alleles)
		if fixed_alleles.count(putative_anc) >= 2:
			return putative_anc

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


def trawl_genome(out_file, chr, f_out4, f_out5, f_ref, var1, var2, var3):
	new_chr = ''
	for line_count, (l_out4, l_out5, l_ref) in enumerate(izip(f_out4, f_out5, f_ref)):
		l_out4 = list(l_out4.upper().rstrip())
		l_out5 = list(l_out5.upper().rstrip())
		l_ref =  list(l_ref.upper().rstrip())

		for ix, (b4, b5, b_ref) in enumerate(izip(l_out4, l_out5, l_ref)):
			# get the current genome position	
			# need to add 1 because python indexes at 0; genomes are at index 1
			pos = line_count * 60 + ix + 1
		
			# this is a variable position
			if pos in var1 or pos in var2 or pos in var3:	
				b1_call = get_history(pos, var1, b_ref)
				b2_call = get_history(pos, var2, b_ref)
				b3_call = get_history(pos, var3, b_ref)
	
				anc_allele = call_ancestral_allele(b1_call, b2_call, b3_call, b4, b5)					
			else:
				anc_allele = b_ref

			new_chr += anc_allele
	return new_chr


def clean_up(chr_out, f_out):
	f_out.close()
	os.remove(chr_out)

def print_chr(chr, new_chr, outfile):
	f = open(outfile, 'w')
	f.write('>%s\n' % chr)
	new_chr = list(new_chr)
	for i in xrange(0, len(new_chr), 60):
		f.write(''.join(new_chr[i:i+60]) + '\n')
	f.close()

def main():
        chr = 'chrZ_random'

	vcf1 = '/mnt/gluster/home/sonal.singhal1/PAR/zf/gatk.ug.unrel_zf.chrZ_random.snps.filtered.vqsr.vcf.gz'
	vcf2 = '/mnt/gluster/home/sonal.singhal1/PAR/ltf/gatk.ug.ltf.chrZ_random.snps.filtered.vqsr.vcf.gz'
	vcf3 = '/mnt/gluster/home/sonal.singhal1/DBF/after_vqsr/by_chr/gatk.ug.dbf.chrZ_random.filtered.coverage.vqsr.vcf.gz'
	out_file = '/mnt/gluster/home/sonal.singhal1/reference/ancestral_sequence.%s.fa' % chr

	# all genomes must be 60 bp per line
	# reference genomes for these two species that are further out
	genome_out4 = '/mnt/gluster/home/sonal.singhal1/Darwin/g_magnirostris/Geospiza_magnirostris.fasta'
	genome_out5 = '/mnt/gluster/home/sonal.singhal1/Darwin/g_fortis/Geospiza_fortis.fasta'
	genome_ref = '/mnt/gluster/home/sonal.singhal1/reference/Zfinch.fa'

	var1 = get_vcf(vcf1, 'in')
	var2 = get_vcf(vcf2, 'out')
	var3 = get_vcf(vcf3, 'out')

	chr_out4, f_out4 = get_chromosome(genome_out4, chr)
	chr_out5, f_out5 = get_chromosome(genome_out5, chr)
	chr_ref, f_ref = get_chromosome(genome_ref, chr)

	new_chr = trawl_genome(out_file, chr, f_out4, f_out5, f_ref, var1, var2, var3)
	print_chr(chr, new_chr, out_file)	

	clean_up(chr_out4, f_out4)
	clean_up(chr_out5, f_out5)
	clean_up(chr_ref, f_ref)

if __name__ == "__main__":
    main()
