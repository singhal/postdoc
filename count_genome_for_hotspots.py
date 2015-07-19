import re
import subprocess
import gzip
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--chr", help="chromosome for which to run analysis")
parser.add_argument("--sp", help="species for which to run analysis")
args = parser.parse_args()
chr = args.chr
sp = args.sp

files = {'chimp': '/data/downloaded/PanMap/SNP-calls/chimp/2010_12_08/callability_mask_2011_01_26/chr%s.mask.fa.gz' % (chr), 
	'human': '/data/downloaded/PanMap/SNP-calls/chimp/SNP-calls-human-mapping/mask/mask_chr%s.hg18_mapping.fasta.gz' % (chr),
	'LTF': '/mnt/gluster/home/sonal.singhal1/LTF/masked_genome/LTF.masked_genome.repeat_masked.fa',
	'ZF': '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.repeat_masked.switch_masked.fa' }

genome = files[sp]

def get_seq_primate(file):
	f = gzip.open(file, 'r')
	header = f.next()

	seq = ''

	for l in f:
		seq += l.rstrip()

	return seq

def get_seq_birds(genome, chr):
	seq = ''
        out = subprocess.Popen('samtools faidx %s chr%s' % (genome, chr), shell=True, stdout=subprocess.PIPE)
        for l in out.stdout:
                if not re.search('>', l):
                        seq += l.rstrip()
	return seq	


if sp in ['human', 'chimp']:
	seq = get_seq_primate(genome)
else:
	seq = get_seq_birds(genome, chr)

good_bases = 0
seq = list(seq)

for i in range(1500, len(seq), 1000):
	subseq = [int(x) for x in seq[i-1000:i+1000]]
	bad = len([x for x in subseq if x > 2])
		
	if bad < 1600:
		good_bases += 1000

print '%s,chr%s,%s' % (sp, chr, good_bases)
	
