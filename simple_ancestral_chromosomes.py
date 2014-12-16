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

# need to work on this
def get_history(pos, var, b, b_ref):
	if pos in var:
		if len(var[pos]) == 1:
			return var[pos].keys()[0]	
		else:
			return var[pos]
	else:
		if b == 'N':
			return 'N'
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
		if fixed_alleles.count(putative_anc) >= 3:
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

def trawl_genome(out_file, chr, f_out1, f_out2, f_out3, f_out4, f_out5, f_ref, var1, var2, var3):
	new_chr = ''
	for line_count, (l_out1, l_out2, l_out3, l_out4, l_out5, l_ref) in enumerate(izip(f_out1, f_out2, f_out3, f_out4, f_out5, f_ref)):
		l_out1 = list(l_out1.rstrip())
		l_out2 = list(l_out2.rstrip())
		l_out3 = list(l_out3.rstrip())
		l_out4 = list(l_out4.upper().rstrip())
		l_out5 = list(l_out5.upper().rstrip())
		l_ref =  list(l_ref.upper().rstrip())

		l_out1 = [int(i) for i in l_out1]
		l_out2 = [int(i) for i in l_out2]
		l_out3 = [int(i) for i in l_out3]

		for ix, (b1, b2, b3, b4, b5, b_ref) in enumerate(izip(l_out1, l_out2, l_out3, l_out4, l_out5, l_ref)):
			# get the current genome position	
			# need to add 1 because python indexes at 0; genomes are at index 1
			pos = line_count * 60 + ix + 1
			anc_allele = 'N'
	
			# this is a variable position
			if pos in var1 or pos in var2 or pos in var3:	
				b1_call = get_history(pos, var1, b1, b_ref)
				b2_call = get_history(pos, var2, b2, b_ref)
				b3_call = get_history(pos, var3, b3, b_ref)
	
				anc_allele = call_ancestral_allele(b1_call, b2_call, b3_call, b4, b5)			
			else:
				if b1 == 'N' and b2 == 'N' and b3 == 'N':
					anc_allele = 'N'
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
	parser = argparse.ArgumentParser()
        parser.add_argument("--chr", help="chromosome for which to run analysis")
        args = parser.parse_args()
        chr = args.chr

	vcf1 = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.filtered.coverage.repeatmasked.vqsr.phased.vcf.gz' % chr
	vcf2 = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.phased.vcf.gz' % chr
	if chr == 'chrZ':
		vcf1 = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chrZ.filtered.coverage.repeatmasked.recodedsex.vqsr.phased.vcf.gz'
		vcf2 = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.chrZ.coverage.repeatmasked.filtered.recodedsex.vqsr.phased.vcf.gz'
	
	vcf3 = '/mnt/gluster/home/sonal.singhal1/DBF/after_vqsr/by_chr/gatk.ug.dbf.%s.filtered.coverage.vqsr.vcf.gz' % chr
	out_file = '/mnt/gluster/home/sonal.singhal1/reference/ancestral_sequence.%s.fa' % chr

	# all genomes must be at 60chr per line!!!
	# masked genome files for each species
	genome_out2 = '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.repeat_masked.fa'
	genome_out1 = '/mnt/gluster/home/sonal.singhal1/LTF/masked_genome/LTF.masked_genome.repeat_masked.fa'
	genome_out3 = '/mnt/gluster/home/sonal.singhal1/DBF/masked_genome/DBF.masked_genome.fa'
	# reference genomes for these two species that are further out
	genome_out4 = '/mnt/gluster/home/sonal.singhal1/Darwin/g_magnirostris/Geospiza_magnirostris.fasta'
	genome_out5 = '/mnt/gluster/home/sonal.singhal1/Darwin/g_fortis/Geospiza_fortis.fasta'
	genome_ref = '/mnt/gluster/home/sonal.singhal1/reference/taeGut1_60.bamorder.fasta'

	var1 = get_vcf(vcf1, 'in')
	var2 = get_vcf(vcf2, 'out')
	var3 = get_vcf(vcf3, 'out')

	chr_out1, f_out1 = get_chromosome(genome_out1, chr)
	chr_out2, f_out2 = get_chromosome(genome_out2, chr)
	chr_out3, f_out3 = get_chromosome(genome_out3, chr)
	chr_out4, f_out4 = get_chromosome(genome_out4, chr)
	chr_out5, f_out5 = get_chromosome(genome_out5, chr)
	chr_ref, f_ref = get_chromosome(genome_ref, chr)

	new_chr = trawl_genome(out_file, chr, f_out1, f_out2, f_out3, f_out4, f_out5, f_ref, var1, var2, var3)
	print_chr(chr, new_chr, out_file)	

	clean_up(chr_out1, f_out1)
	clean_up(chr_out2, f_out2)
	clean_up(chr_out3, f_out3)
	clean_up(chr_out4, f_out4)
	clean_up(chr_out5, f_out5)
	clean_up(chr_ref, f_ref)

if __name__ == "__main__":
    main()
