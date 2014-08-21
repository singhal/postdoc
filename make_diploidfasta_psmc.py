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
                                        allele_counts[alleles[i]] = len(re.findall('%s\/' % i, l)) + len(re.findall('\/%s' % i, l))
                                tot_n = sum(allele_counts.values())
                                
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


def get_chromosome(genome, chr):
        outfile = genome + '_' + chr
        subprocess.call('~/bin/samtools-0.1.19/samtools faidx %s %s > %s' % (genome, chr, outfile), shell=True)
        out_f = open(outfile, 'r')
        chromosome = ''
        locus_name = out_f.next()
        for l in out_f:
                chromosome = chromosome + l.rstrip()
        out_f.close()
        os.remove(outfile)
        return list(chromosome)


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


def main():
	parser = argparse.ArgumentParser()
        parser.add_argument("--chr", help="chromosome for which to run analysis")
        args = parser.parse_args()
        chr = args.chr

	out_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/PSMC/'
	vcf_in = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.vqsr.vcf.gz' % chr
	masked_genome = '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.fa'
	genome_ref = '/mnt/gluster/home/sonal.singhal1/reference/taeGut1_60.bamorder.fasta'

	var_in = get_vcf(vcf_in)

	chr_out1 = get_chromosome(genome_out1, chr)
	chr_ref = get_chromosome(genome_ref, chr)

	make_seq(var_in, chr_out1, chr_ref, out_dir)

if __name__ == "__main__":
    main()
