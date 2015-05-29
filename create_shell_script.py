import os

chrs = [ 'chr1', 'chr1A', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
         'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
         'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ']

long_chrs = [ 'chr1', 'chr1A', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15']

long_chrs = [ 'chr1', 'chr1A', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chrZ']

for chr in long_chrs:
	print "echo \'python get_qual.py --chr %s\' | qsub -l h_vmem=5g -cwd -V -j y -N %s" % (chr, chr)

