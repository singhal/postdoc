import argparse
import subprocess
import re

def get_chromosome(genome, chr):
	seq = ''
	out = subprocess.Popen('samtools faidx %s %s' % (genome, chr), shell=True, stdout=subprocess.PIPE)
	for l in out.stdout:
		if not re.match('>', l):
			seq += l.rstrip().upper()
	return list(seq)


def get_divergence(chr1, chr2):
	hets = { 'M': ['A', 'C'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['C', 'G'], 'Y': ['C', 'T'], 'K': ['G', 'T'],
		 'A': ['A', 'A'], 'G': ['G', 'G'], 'C': ['C', 'C'], 'T': ['T', 'T'] }
	legal = ['A', 'T', 'C', 'G']

	denom = diff = 0

	for bp1, bp2 in zip(chr1, chr2):
		if bp1 != 'N' and bp2 != 'N':
			denom += 1
			if bp1 in legal and bp2 in legal:
				if bp1 != bp2:
					diff += 1
			else:
				bases1 = hets[bp1]
				bases2 = hets[bp2]

				tmp_diff = 0

				for a in bases1:
					for b in bases2:
						if a != b:
							tmp_diff += 1
				diff += tmp_diff / 4.

	return (denom, diff)		
				

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--sp1", help="species 1 for which to run analysis")
	parser.add_argument("--sp2", help="species 2 for which to run analysis")
	args = parser.parse_args()
	sp1 = args.sp1
	sp2 = args.sp2

	genomes = {     'LTF': '/mnt/gluster/home/sonal.singhal1/reference/LTF_reference.fa',
                        'ZF': '/mnt/gluster/home/sonal.singhal1/reference/ZF_reference.fa',
                        'DBF': '/mnt/gluster/home/sonal.singhal1/reference/DBF_reference.fa',
                        'Gfor': '/mnt/gluster/home/sonal.singhal1/Darwin/g_fortis/Geospiza_fortis.fasta',
                        'Fic': '/mnt/gluster/home/sonal.singhal1/ficedula/Ficedula_albicollis.fasta'
                }

	genome1 = genomes[sp1]
	genome2 = genomes[sp2]

	chrs = [ 'chr1', 'chr1A', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         	'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
         	'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
         	'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ']

	out = '%s_%s.divergence.csv' % (sp1, sp2)
	o = open(out, 'w')
	o.write('chr,denominator,num_different,percent_different\n')
	for chr in chrs:
		chr1 = get_chromosome(genome1, chr)
		chr2 = get_chromosome(genome2, chr)

		denom, diff = get_divergence(chr1, chr2)
		o.write('%s,%s,%s,%s\n' % (chr, denom, diff, (diff / float(denom))))
	o.close()


if __name__ == "__main__":
    main()
