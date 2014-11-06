import glob
import re
import os
import subprocess
import gzip

# chrs = [ 'chr1', 'chr1A', 'chr1B', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
#        'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
#         'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
#         'chrLG2', 'chrLG5', 'chrLGE22']

chrs = ['chr1', 'chr1A', 'chr1B', 'chr2']

ref = '/mnt/gluster/home/sonal.singhal1/reference/taeGut1_60.bamorder.fasta'
mask_ZF = '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.fa'
mask_LTF = '/mnt/gluster/home/sonal.singhal1/LTF/masked_genome/LTF.masked_genome.fa'

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

def get_ancestral_allele(file, var):
	f = open(file, 'r')
	header = f.next()
	for l in f:
		d = re.split(',', l)
		if int(d[1]) in var:
			if not var[int(d[1])]['ref']:
				if d[3] == d[4]:
					if d[3] in ['A', 'T', 'C', 'G']:
						if re.search(d[3], d[2]):
							var[int(d[1])]['ref'] = var[int(d[1])]['triplet'][0] + d[3] + var[int(d[1])]['triplet'][2]
							bases = re.findall('(\S)\:', d[2])
							bases.remove(d[3])
							var[int(d[1])]['out'] = var[int(d[1])]['triplet'][0] + bases[0] + var[int(d[1])]['triplet'][2]
	f.close()
	return var


def get_zf_var(zf_vcf, var):
	vcf = gzip.open(zf_vcf, 'r')
	for l in vcf:
		if not re.search('^#', l):
			d = re.split('\t', l)
			pos = int(d[1])
			if pos in var:
				alt = len(re.findall('1\/', l)) + len(re.findall('\/1', l))
				ref = len(re.findall('0\/', l)) + len(re.findall('\/0', l))
				
				var[pos]['allele1']['bp'] = d[3]
				var[pos]['allele2']['bp'] = d[4]

				var[pos]['allele1']['ZF'] = ref
				var[pos]['allele2']['ZF'] = alt
	return var


def get_ltf_var(ltf_vcf, var):
        vcf = gzip.open(ltf_vcf, 'r')
        for l in vcf:
                if not re.search('^#', l):
                        d = re.split('\t', l)
                        pos = int(d[1])
                        if pos in var:
				if not var[pos]['allele1']['bp']:
					var[pos]['allele1']['bp'] = d[3]
	                                var[pos]['allele2']['bp'] = d[4]

                                	alt1 = len(re.findall('1\/', '\t'.join(d[9:19]))) + len(re.findall('\/1', '\t'.join(d[9:19])))
                                	alt2 = len(re.findall('1\/', '\t'.join(d[19:29]))) + len(re.findall('\/1', '\t'.join(d[19:29])))

					ref1 = len(re.findall('0\/', '\t'.join(d[9:19]))) + len(re.findall('\/0', '\t'.join(d[9:19])))
                                	ref2 = len(re.findall('0\/', '\t'.join(d[19:29]))) + len(re.findall('\/0', '\t'.join(d[19:29])))
                                
                                	var[pos]['allele1']['LTFh'] = ref1
                                	var[pos]['allele2']['LTFh'] = alt1

					var[pos]['allele1']['LTFa'] = ref2
                                	var[pos]['allele2']['LTFa'] = alt2
				elif var[pos]['allele1']['bp'] == d[3] and var[pos]['allele2']['bp'] == d[4]:
					alt1 = len(re.findall('1\/', '\t'.join(d[9:19]))) + len(re.findall('\/1', '\t'.join(d[9:19])))
                                        alt2 = len(re.findall('1\/', '\t'.join(d[19:29]))) + len(re.findall('\/1', '\t'.join(d[19:29])))

                                        ref1 = len(re.findall('0\/', '\t'.join(d[9:19]))) + len(re.findall('\/0', '\t'.join(d[9:19])))
                                        ref2 = len(re.findall('0\/', '\t'.join(d[19:29]))) + len(re.findall('\/0', '\t'.join(d[19:29])))

                                        var[pos]['allele1']['LTFh'] = ref1
                                        var[pos]['allele2']['LTFh'] = alt1

                                        var[pos]['allele1']['LTFa'] = ref2
                                        var[pos]['allele2']['LTFa'] = alt2
	vcf.close()
        return var


out_file = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/dadi/all_species.dadi_leftover.txt'
out = open(out_file, 'w')
out.write('In\tOut\tAllele1\tZF\tLTFh\tLTFa\tAllele2\tZF\tLTFh\tLTFa\tGene\tPosition\n')

for chr in chrs:
	var = {}

	chr_ref = get_chromosome(ref, chr)
	chr_ZF = get_chromosome(mask_ZF, chr)
	chr_LTF = get_chromosome(mask_LTF, chr)

	for pos, (bp1, bp2) in enumerate(zip(chr_ZF, chr_LTF)):
		bp1 = int(bp1)
		bp2 = int(bp2)
	
		if bp1 == 1 or bp2 == 1:
			if bp1 < 4 and bp2 < 4:
				triplet = chr_ref[(pos - 1):(pos+2)]
				var[pos+1] = {	'ref': '', 'out': '', 'triplet': triplet, 
						'allele1': {'bp': '', 'ZF': 38, 'LTFa': 20, 'LTFh': 20}, 
						'allele2': {'bp': '', 'ZF': 0, 'LTFa': 0, 'LTFh': 0} }

	del chr_ref
	del chr_ZF
	del chr_LTF

	zf_aa = '/mnt/gluster/home/sonal.singhal1/ZF/ancestral_allele/ancestral_allele.%s.csv' % chr
	ltf_aa = '/mnt/gluster/home/sonal.singhal1/LTF/ancestral_allele/ancestral_allele.%s.csv' % chr

	var = get_ancestral_allele(zf_aa, var)
	var = get_ancestral_allele(ltf_aa, var)	

	zf_vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.filtered.recoded_biallelicSNPs.nomendel.vcf.gz' % chr
	ltf_vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.filtered.coverage.recoded_biallelicSNPs.vqsr.vcf.gz' % chr
	var = get_zf_var(zf_vcf, var)
	var = get_ltf_var(ltf_vcf, var)
	
	for pos in sorted(var.keys()):
		if var[pos]['ref']:	
			out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (var[pos]['ref'], var[pos]['out'], var[pos]['allele1']['bp'], var[pos]['allele1']['ZF'], var[pos]['allele1']['LTFh'], var[pos]['allele1']['LTFa'], var[pos]['allele2']['bp'], var[pos]['allele2']['ZF'], var[pos]['allele2']['LTFh'], var[pos]['allele2']['LTFa'], chr, pos))

out.close()
