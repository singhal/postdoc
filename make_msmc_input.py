import re
import glob
import subprocess
import os

species = 'LTF'
files = sorted(glob.glob('/mnt/gluster/home/sonal.singhal1/%s/phasing/PIR_approach/results/*haps' % species))
# these are the haplotypes to sample
# for ZF, just take the first four haplotypes
haps = [5, 6, 7, 8]
# for LTF, these are the first two and last two haplotypes
haps = [[5, 6,  43, 44], [7, 8, 41, 42], [9, 10, 39, 40], [11, 12, 37, 38], [13, 14, 35, 36]]
out_dir = '/mnt/gluster/home/sonal.singhal1/%s/analysis/msmc/' % species
masked_genome = '/mnt/gluster/home/sonal.singhal1/%s/masked_genome/%s.masked_genome.fa' % (species, species)

def get_chromosome(genome, chr):
        outfile = genome + '_' + chr
        subprocess.call('~/bin/samtools-0.1.19/samtools faidx %s %s > %s' % (genome, chr, outfile), shell=True)
        out_f = open(outfile, 'r')
        chromosome = ''
        locus_name = out_f.next()
        for l in out_f:
                chromosome = chromosome + l.rstrip().upper()
        out_f.close()
        os.remove(outfile)
        return list(chromosome)

def allele(nuc, ref, alt):
	if nuc == '0':
		return ref
	else:
		return alt

for file in files:
	chr = re.search('(chr\S+)_', file).group(1)
	chr_id = chr.replace('chr','')
	if chr != 'chrZ':
		chromosome = get_chromosome(masked_genome, chr)
		for ix, hap in enumerate(haps):
			out_file = '%s%s_%s.%s.msmc' % (out_dir, species, chr, ix)
			out = open(out_file, 'w')
			f = open(file, 'r')
			last_site = 1
			for l in f:
				d = re.split('\s+', l.rstrip())
				alleles = [d[ind] for ind in hap]
				if len(set(alleles)) > 1:
					alleles = [allele(x, d[3], d[4]) for x in alleles]
					cur_site = int(d[2])
					# don't need to change indices despite diff in base 0 and base 1
					# this is because we want the site after last_site up to the cur_site
					# construction gets the number of non-masked bases between two sites
					num_bases_btn = len(filter(lambda x: int(x) < 4, chromosome[last_site:cur_site]))
					out.write('%s\t%s\t%s\t%s\n' % (chr_id, d[2], num_bases_btn, ''.join(alleles)))
					last_site = int(d[2])
			f.close()
			out.close()
