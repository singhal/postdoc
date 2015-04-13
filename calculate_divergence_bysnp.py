import re
import gzip
import glob
import subprocess
import argparse
import random
from itertools import izip

def get_variants(vcf, sp, chr):
	haps = {}
	f = gzip.open(vcf, 'r')	
	nind = None

	for l in f:
		if not re.search('#', l):

			d = re.split('\t', l.rstrip())

			# check not indel
			indel = False
			alleles = [d[3]] + re.split(',', d[4])
			for allele in alleles:
				if len(allele) > 1:
					indel = True

			num_alleles = {}
			for ix, allele in enumerate(alleles):
				num_alleles[str(ix)] = allele
			num_alleles['.'] = 'N'

			if not indel:
                        	genos = []
                                for geno in d[9:]:
                                        geno = re.search('^([^:]+)', geno).group(1)
                                        if re.search('/', geno):
                                                geno = re.split('/', geno)
                                                random.shuffle(geno)
                                                genos += geno
                                        else:
                                                genos += re.split('\|', geno)
				genos = [num_alleles[x] for x in genos]

				if sp == 'LTFh':
					if chr == 'chrZ':
						genos = genos[0:15]
					else:
						genos = genos[0:20]
				if sp == 'LTFa':
					if chr == 'chrZ':
						genos = genos[15:32]
					else:
						genos = genos[20:40]                        

			        haps[int(d[1])] = genos
                                nind = len(genos)

	f.close()

	return haps, nind


def get_sequence(ref_genome, genome1, genome2, chr, start, end):
	seq_call = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome1, chr, start, end), shell=True, stdout=subprocess.PIPE)
	seq1 = ''
	for l in seq_call.stdout:
		if not re.match('>', l):
			seq1 += l.rstrip().upper()
	seq1 = list(seq1)

	seq_call = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome2, chr, start, end), shell=True, stdout=subprocess.PIPE)
        seq2 = ''
        for l in seq_call.stdout:
                if not re.match('>', l):
                        seq2 += l.rstrip().upper()
	seq2 = list(seq2)

	seq_call = subprocess.Popen('samtools faidx %s %s:%s-%s' % (ref_genome, chr, start, end), shell=True, stdout=subprocess.PIPE)
        ref_seq = ''
        for l in seq_call.stdout:
                if not re.match('>', l):
                        ref_seq += l.rstrip().upper()
        ref_seq = list(ref_seq)

	keep_sites = {}
	for ix, (ref, bp1, bp2) in enumerate(izip(ref_seq, seq1, seq2)):
		pos = ix + start
		if bp1 in ['0', '1', '2', '3'] and bp2 in ['0', '1', '2', '3']:
			keep_sites[pos] = ref

	return keep_sites


def get_divergence(start, end, hap1, nind1, hap2, nind2, keep_sites):
	sites = sorted(set(hap1.keys() + hap2.keys()))
	sites = [x for x in sites if x in keep_sites]
	
	pi = 0
	pi_sum = 0

	for site in sites:
		num_diff = 0
		if site in hap1:
			geno1 = hap1[site]
		else:
			geno1 = keep_sites[site]
		if site in hap2:
			geno2 = hap2[site]
		else:
			geno2 = keep_sites[site]

		for i in geno1:
			for j in geno2:
				if i != j and i != 'N' and j != 'N':
					num_diff += 1

		pi_sum += num_diff / float(len(geno1) * len(geno2))		

	if  len(keep_sites) > 0:
		pi = pi_sum / float(len(keep_sites))

	return pi


def get_divergence_within(start, end, hap, nind, keep_sites):
        sites = sorted(hap.keys())
        sites = [x for x in sites if x in keep_sites]

	pi = 0
        pi_sum = 0

        for site in sites:
                num_diff = 0
                geno = hap[site]
             
                for ix, i in enumerate(geno):
                        for j in geno[ix:]:
                                if i != j and i != 'N' and j != 'N':
					num_diff += 1

                pi_sum += (2 * num_diff) / float(len(geno) * len(geno))

        if  len(keep_sites) > 0:
                pi = pi_sum / float(len(keep_sites))	

        return pi


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--sp1", help="species1 for which to run analysis")
	parser.add_argument("--sp2", help="species2 for which to run analysis")
	parser.add_argument("--chr", help="chromosome for which to run anlaysis")
	args = parser.parse_args()
	sp1 = args.sp1
	sp2 = args.sp2
	chr = args.chr

	window = 10000
	
	chr_lengths = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                        'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                        'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                        'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                        'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                        'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                        'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                        'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                        'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351}

	vcfs = 	{	'ZF': '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.shared.noswitch.phased.vqsr2.vcf.gz' % chr,
			'LTF': '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.coverage.repeatmasked.filtered.phased.vqsr2.vcf.gz' % chr,
			'DBF': '/mnt/gluster/home/sonal.singhal1/DBF/after_vqsr/by_chr/gatk.ug.dbf.%s.filtered.coverage.vqsr.vcf.gz' % chr,
			'LTFh': '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.coverage.repeatmasked.filtered.phased.vqsr2.vcf.gz' % chr,
			'LTFa': '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.coverage.repeatmasked.filtered.phased.vqsr2.vcf.gz' % chr
		}
	
	if chr == 'chrZ':
		vcfs =  {       'ZF': '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.chrZ.coverage.repeatmasked.filtered.nomendel.shared.recodedsex.phased.vqsr2.vcf.gz',
                        'LTF': '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chrZ.coverage.repeatmasked.filtered.recodedsex.phased.vqsr2.vcf.gz',
                        'DBF': '/mnt/gluster/home/sonal.singhal1/DBF/after_vqsr/by_chr/gatk.ug.dbf.chrZ.filtered.coverage.vqsr.vcf.gz',
                        'LTFh': '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chrZ.coverage.repeatmasked.filtered.recodedsex.phased.vqsr2.vcf.gz',
                        'LTFa': '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chrZ.coverage.repeatmasked.filtered.recodedsex.phased.vqsr2.vcf.gz'
                }

	genomes = {	'ZF': '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.repeat_masked.switch_masked.fa', 
			'LTF': '/mnt/gluster/home/sonal.singhal1/LTF/masked_genome/LTF.masked_genome.repeat_masked.fa',
			'LTFh': '/mnt/gluster/home/sonal.singhal1/LTF/masked_genome/LTF.masked_genome.repeat_masked.fa',
			'LTFa': '/mnt/gluster/home/sonal.singhal1/LTF/masked_genome/LTF.masked_genome.repeat_masked.fa',
			'DBF': '/mnt/gluster/home/sonal.singhal1/DBF/masked_genome/DBF.masked_genome.fa'
		}
	
	ref_genome = '/mnt/gluster/home/sonal.singhal1/reference/Zfinch.fa'
	out_file = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/pop_gen/%s.%s_%s.pairwise_snps_divergence.csv' % (chr, sp1, sp2)

	(hap1, nind1) = get_variants(vcfs[sp1], sp1, chr)
	(hap2, nind2) = get_variants(vcfs[sp2], sp2, chr)

	o = open(out_file, 'w')
	o.write('chr,start,end,num_sites,dxy,dx,dy,da\n')
	for start in range(1, chr_lengths[chr] + 1, window):
		end = start + window - 1
		if end > chr_lengths[chr]:
			end = chr_lengths[chr]

		keep_sites = get_sequence(ref_genome, genomes[sp1], genomes[sp2], chr, start, end)
		dxy = get_divergence(start, end, hap1, nind1, hap2, nind2, keep_sites)
		dx = get_divergence_within(start, end, hap1, nind1, keep_sites)
		dy = get_divergence_within(start, end, hap2, nind2, keep_sites)
		da = dxy - (dx + dy) / 2.0
		
		o.write('%s,%s,%s,%s,%s,%s,%s,%s\n' % (chr, start, end, len(keep_sites), dxy, dx, dy, da))
	o.close()


if __name__ == "__main__":
    main()
