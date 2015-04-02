import re
import gzip
import glob
import subprocess
import argparse
import random

def get_variants(vcf):
	var = {}
	f = gzip.open(vcf, 'r')
	
	for l in f:
		if not re.search('#', l):

			d = re.split('\t', l.rstrip())

			# check not indel
			indel = False
			alleles = [d[3]] + re.split(',', d[4])
			for allele in alleles:
				if len(allele) > 1:
					indel = True

			if not indel:
				# check number of alleles
				genos = []
				for geno in d[9:]:
					geno = re.search('^([^:]+)', geno).group(1)
					genos += re.split('[|/]', geno)
				genos = [x for x in genos if not re.search('\.', x)]

				var[int(d[1])] = {}
				for ix, allele in enumerate(alleles):
					if genos.count(str(ix)) > 0:
						var[int(d[1])][allele] = genos.count(str(ix)) / float(len(genos))

	f.close()

	return var


def get_sequence(genome, chr, start, end):
	seq_call = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome, chr, start, end), shell=True, stdout=subprocess.PIPE)
	seq = ''
	for l in seq_call.stdout:
		if not re.match('>', l):
			seq += l.rstrip().upper()

	return list(seq)


def get_counts(seq):
	ancAT = seq.count('A') + seq.count('T')
	ancGC = seq.count('G') + seq.count('C')

	return (ancAT, ancGC)


def get_transitions(start, end, var1, var2, var3, aa_seq, ref_seq):
        counts = {'AT_AT': 0, 'AT_GC': 0, 'GC_GC': 0, 'GC_AT': 0}

        types = {       'A': {'A': 'AT_AT', 'G': 'AT_GC', 'T': 'AT_AT', 'C': 'AT_GC'},
                        'C': {'A': 'GC_AT', 'G': 'GC_GC', 'T': 'GC_AT', 'C': 'GC_GC'},
                        'T': {'A': 'AT_AT', 'G': 'AT_GC', 'T': 'AT_AT', 'C': 'AT_GC'},
                        'G': {'A': 'GC_AT', 'G': 'GC_GC', 'T': 'GC_AT', 'C': 'GC_GC'},
                }

        for ix, (a_bp, r_bp) in enumerate(zip(aa_seq, ref_seq)):
                pos = start + ix

		# no point in going through this if no ancestral identification
		if a_bp != 'N':
			allele = 'N'
			if pos in var1:
				# fixed in focal lineage
				if len(var1[pos]) == 1:
					allele = var1[pos].keys()[0]
			else:
				allele = r_bp

			# only interested in fixed subsitutions
			if allele != 'N':
				others = [a_bp]
	
				if pos in var2:
					if len(var2[pos]) == 1:
						# okay, fixed!
						others.append(var2[pos].keys()[0])
				else:
					others.append(r_bp)

				if pos in var3:
					if len(var3[pos]) == 1:
						# okay fixed
						others.append(var3[pos].keys()[0])
				else:
					others.append(r_bp)

				# if len == 3, then not polymorphic in any other species
				if len(others) == 3:
					# if lenght of set is 1, then all the other alleles are the same
					if len( set( others ) ) == 1:
						# check if allele in our base is diff from theirs
						if others[0] != allele:
							counts[types[a_bp][allele]] += 1
	return counts


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--sp", help="species for which to run analysis")
	parser.add_argument("--chr", help="chromosome for which to run anlaysis")
	args = parser.parse_args()
	sp = args.sp
	chr = args.chr

	window = 100000
	chr_lengths = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                        'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                        'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                        'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                        'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                        'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                        'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                        'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                        'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351}
	vcf_bases = {	'ZF': '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.shared.noswitch.phased.vqsr2.vcf.gz' % chr,
			'LTF': '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.coverage.repeatmasked.filtered.phased.vqsr2.vcf.gz' % chr,
			'DBF': '/mnt/gluster/home/sonal.singhal1/DBF/after_vqsr/by_chr/gatk.ug.dbf.%s.filtered.coverage.vqsr.vcf.gz' % chr
		}

	if chr == 'chrZ':
                vcf_bases =  { 'ZF': '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.chrZ.coverage.repeatmasked.filtered.nomendel.shared.recodedsex.phased.vqsr2.vcf.gz',
                        'LTF': '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chrZ.coverage.repeatmasked.filtered.recodedsex.phased.vqsr2.vcf.gz',
                        'DBF': '/mnt/gluster/home/sonal.singhal1/DBF/after_vqsr/by_chr/gatk.ug.dbf.chrZ.filtered.coverage.vqsr.vcf.gz'
                }

	ref_genome = '/mnt/gluster/home/sonal.singhal1/reference/Zfinch.fa'
	aa_genome = '/mnt/gluster/home/sonal.singhal1/reference/ancestral_genome.fa'
	out_file = '/mnt/gluster/home/sonal.singhal1/%s/analysis/mut_skew/%s.gcstar.csv' % (sp, chr)

	vcf1 = vcf_bases[sp]
	not_sp = [x for x in vcf_bases.keys() if x != sp]
	vcf2 = vcf_bases[not_sp[0]]
	vcf3 = vcf_bases[not_sp[1]]

	var1 = get_variants(vcf1)
	var2 = get_variants(vcf2)
	var3 = get_variants(vcf3)

	o = open(out_file, 'w')
	o.write('chr,start,end,gcstart,AT_GC,GC_AT,ancAT,ancGC\n')
	for start in range(1, chr_lengths[chr] + 1, window):
		end = start + window - 1
		if end > chr_lengths[chr]:
			end = chr_lengths[chr]
		aa_seq = get_sequence(aa_genome, chr, start, end)
		ref_seq = get_sequence(ref_genome, chr, start, end)

		(ancAT, ancGC) = get_counts(aa_seq)
		counts = get_transitions(start, end, var1, var2, var3, aa_seq, ref_seq)
		
		if ancAT > 0 and ancGC > 0 and (counts['AT_GC'] + counts['GC_AT']) > 0:
			gcstar = (counts['AT_GC'] / float(ancAT)) / ((counts['AT_GC'] / float(ancAT))+(counts['GC_AT'] / float(ancGC)))
		else:
			gcstar = 'NA'
		o.write('%s,%s,%s,%s,%s,%s,%s,%s\n' % (chr, start, end, gcstar, counts['AT_GC'], counts['GC_AT'], ancAT, ancGC))
	o.close()

if __name__ == "__main__":
    main()
