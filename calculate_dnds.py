import re
import gzip
import pandas as pd
import argparse
from itertools import izip
import subprocess
import os

def get_vcf(vcf, start, stop):
        var = {}
        f = gzip.open(vcf, 'r')
        for l in f:
                if not re.search('^#', l):
                        d = re.split('\t', l.rstrip())
                        
                        alleles = [d[3]] + re.split(',', d[4])
                        pos = int(d[1])
                        indel = False
                        
                        for allele in alleles:
                                if len(allele) > 1:
                                        indel = True

                        if not indel:
                                genos = []
                                for geno in d[9:]:
                                        geno = re.search('^([^:]+)', geno).group(1)
					genos.append(geno)
				genos = genos[start:stop]
				genos2 = []
				for geno in genos:
					genos2 += re.split('[|/]', geno)
				
				for ix, allele in enumerate(alleles):
					af = genos2.count(str(ix)) / float(len(genos2))
					if af >= 0.5 and allele != d[3]:
						var[int(d[1])] = allele
	f.close()
	return var


def get_genes(gff, chr):
	gff = pd.read_csv(gff, sep='\t', header=None, names=['chr', 'type', 'cds_mrna', 
                                                'start', 'stop', 'score', 
                                                'orientation', 'codon_pos', 'id'])
	gff['chr'] = 'chr' + gff.chr
	gff = gff[gff.chr == chr]
	gff['id'] = [x.replace('ID=', '') for x in gff.id]
	gff['id'] = [x.replace('Parent=', '') for x in gff.id]
	gff['id'] = [x.replace(';', '') for x in gff.id]

	genes = {}
	mrna = gff[gff.cds_mrna == 'mRNA']
        for id, orient in izip(mrna.id, mrna.orientation):
                genes[id] = {'orient': orient, 'exons': []}

	for id, group in gff[gff.cds_mrna == 'CDS'].groupby(['id']):
        	for start, stop in izip(group.start, group.stop):
                	genes[id]['exons'].append([start,stop])

	return genes


def get_seq(genome, chr):
	seqcall = subprocess.Popen('samtools faidx %s %s' % (genome, chr), shell=True, stdout=subprocess.PIPE)
	seq = ''
	for l in seqcall.stdout:
		if not re.match('>', l):
			seq += l.rstrip().upper()
	return list(seq)

def mutate_seq(seq, var):
	new_seq = list(seq)
	for pos in var:
		ix = pos - 1
		new_seq[ix] = var[pos]
	var = {}
	return new_seq

		
def rev_comp(seq):
	complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
	reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
	return reverse_complement


def translate(seq):
        map = { 'TTT':'F', 'TTC':'F', 'TTA':'L', 'TTG':'L',
                'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S',
                'TAT':'Y', 'TAC':'Y', 'TAA':'*', 'TAG':'*',
                'TGT':'C', 'TGC':'C', 'TGA':'*', 'TGG':'W',
                'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L',
                'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P',
                'CAT':'H', 'CAC':'H', 'CAA':'Q', 'CAG':'Q',
                'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R',
                'ATT':'I', 'ATC':'I', 'ATA':'I', 'ATG':'M',
                'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T',
                'AAT':'N', 'AAC':'N', 'AAA':'K', 'AAG':'K',
                'AGT':'S', 'AGC':'S', 'AGA':'R', 'AGG':'R',
                'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V',
                'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A',
                'GAT':'D', 'GAC':'D', 'GAA':'E', 'GAG':'E',
                'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G'}

        triplets = [seq[i:i+3] for i in range(0, len(seq), 3)]
        aa = ''
        for triplet in triplets:
                if triplet in map:
                        aa += map[triplet]
                else:
                        aa += 'X'
        return aa


def get_cds(cds, seq, gene, name):
	new_cds = ''
	for start, end in gene['exons']:
		new_cds += ''.join(seq[start - 1:end])
	if gene['orient'] == '-':
		new_cds = rev_comp(new_cds)
	aa = translate(new_cds)
	if not re.search('\S\*\S', aa):
		cds[name] = new_cds
	else:
		print '%s %s' % (gene, name)

	return cds


def print_genes(out_dir, genes, seq, zf_var, ltfa_var, ltfh_var, dbf_var):
	zf_seq = mutate_seq(seq, zf_var)
	ltfa_seq = mutate_seq(seq, ltfa_var)
	ltfh_seq = mutate_seq(seq, ltfh_var)
	dbf_seq = mutate_seq(seq, dbf_var)

	for gene in genes:
		cds = {}
		cds = get_cds(cds, zf_seq, genes[gene], 'zf')
		cds = get_cds(cds, ltfa_seq, genes[gene], 'ltfa')
		cds = get_cds(cds, ltfh_seq, genes[gene], 'ltfh')
		cds = get_cds(cds, dbf_seq, genes[gene], 'dbf')

		if len(cds) > 1:
			out = '%s/alignments/%s.phy' % (out_dir, gene)
			o = open(out, 'w')
			o.write(' %s %s\n' % (len(cds), len(cds.values()[0])))
			for id, seq in cds.items():
				id = id + ' ' * (10 - len(id))
				o.write('%s%s\n' % (id, seq))
			o.close()
		       		

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--chr", help="chromosome for which to run analysis")
	args = parser.parse_args()
	chr = args.chr

	# result file
	out_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/dnds/'

	# genome file
	genome = '/mnt/gluster/home/sonal.singhal1/reference/Zfinch.fa'

	# gff file
	gff = '/mnt/gluster/home/sonal.singhal1/reference/Taeniopygia_guttata.gff'

	# vcf files
	zf_vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.shared.noswitch.phased.vqsr2.vcf.gz' % chr
	if chr == 'chrZ':
        	zf_vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.chrZ.coverage.repeatmasked.filtered.nomendel.shared.recodedsex.phased.vqsr2.vcf.gz'
	ltf_vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.coverage.repeatmasked.filtered.phased.vqsr2.vcf.gz'  % chr
	if chr == 'chrZ':
        	ltf_vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chrZ.coverage.repeatmasked.filtered.recodedsex.phased.vqsr2.vcf.gz'
	dbf_vcf = '/mnt/gluster/home/sonal.singhal1/DBF/after_vqsr/by_chr/gatk.ug.dbf.%s.filtered.coverage.vqsr.vcf.gz' % chr
	
	zf_var = get_vcf(zf_vcf, 0, 19)
	ltfa_var = get_vcf(ltf_vcf, 10, 20)
	ltfh_var = get_vcf(ltf_vcf, 0, 10)
	dbf_var = get_vcf(dbf_vcf, 0, 1)

	genes = get_genes(gff, chr)
	seq = get_seq(genome, chr)
	print_genes(out_dir, genes, seq, zf_var, ltfa_var, ltfh_var, dbf_var)

if __name__ == "__main__":
    main()
