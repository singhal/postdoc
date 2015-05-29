import re
import subprocess
from itertools import izip
import argparse
import subprocess

def get_chr(genomes, chr):
	chrs = {}
	# get chromosome sequence
	for sp, genome in genomes.items():
		seq = ''
		seqcall = subprocess.Popen('samtools faidx %s %s' % (genome, chr), shell=True, stdout=subprocess.PIPE)
		for l in seqcall.stdout:
                        if not re.match('>', l):
				seq += l.upper().rstrip()
		chrs[sp] = seq
	
	# turn into a list
	for sp, seq in chrs.items():
		chrs[sp] = list(seq)

	return chrs


def most_common(lst):
    return max(set(lst), key=lst.count)


def get_ancestral_bp(bases):
        anc_base = 'N'

        bases = [x for x in bases if x in ['A', 'T', 'C', 'G']]

        # can't call an ancestral base
        # unless there is more than one base
        if len(bases) > 1:
                # all bases are the same
                if len(set(bases)) == 1:
                        anc_base = bases[0]
                else:
                        putative_anc = most_common(bases)
                        if bases.count(putative_anc) / float(len(bases)) > 0.5:
                                anc_base = putative_anc

        return anc_base


def get_anc(seq):
	anc = ''

	# get ancestral base
	for bases in izip(seq['ZF'], seq['LTF'], seq['DBF'], seq['Gfor'], seq['Gmag'], seq['Fic']):
                bases = list(bases)
                anc_base = get_ancestral_bp(bases)
                anc += anc_base
	# clear out the hash 
	seq = {}

	return list(anc)


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

        genomes = {     'LTF': '/mnt/gluster/home/sonal.singhal1/reference/LTF_reference.fa',
                        'ZF': '/mnt/gluster/home/sonal.singhal1/reference/ZF_reference.fa',
                        'DBF': '/mnt/gluster/home/sonal.singhal1/reference/DBF_reference.fa',
                        'Gfor': '/mnt/gluster/home/sonal.singhal1/Darwin/g_fortis/Geospiza_fortis.fasta',
                        'Gmag': '/mnt/gluster/home/sonal.singhal1/Darwin/g_magnirostris/Geospiza_magnirostris.fasta',
			'Fic': '/mnt/gluster/home/sonal.singhal1/ficedula/Ficedula_albicollis.fasta'
                }
	
	chrs = get_chr(genomes, chr)
	anc = get_anc(chrs)
	outfile = '/mnt/gluster/home/sonal.singhal1/reference/all_birds.ancestral_%s.fa' % chr
	print_chr(chr, anc, outfile)	
	print '%s is done!' % chr

if __name__ == "__main__":
    main()
