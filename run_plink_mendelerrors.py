import subprocess
import re

#chrs = [ 'chr1', 'chr1A', 'chr1B', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
#         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
#         'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
#         'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ' ]

chrs = ['chr1', 'chr2', 'chr3']

dir  = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/mendelian_errors/'

for chr in chrs:
	data_file = '%sall_zf.%s' % (dir, chr)
	out = '%sall_zf.me.%s' % (dir, chr)
	
	subprocess.call("echo '~/bin/plink-1.07-x86_64/plink --noweb --file %s --mendel --out %s' | qsub -l h_vmem=5g -cwd -V -j y -N \'me_%s\'" % (data_file, out, chr), shell = True)
