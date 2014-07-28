import subprocess
import glob
import re

dir = '/mnt/lustre/home/sonal.singhal1/ZF/phasing/family_approach/'

# chrs = [ 'chr16', 'chr1A', 'chr1B', 'chr2', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
#         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr1', 'chr17', 'chr18',
#         'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
#         'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ' ]

chrs_long = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr1A', 'chr1B']

chrs = ['chr26']

for chr in chrs:
	ped_file = dir + 'ped_map_files/all_zf.%s.ped' % chr
	map_file = dir + 'ped_map_files/all_zf.%s.map' % chr
	out_file = dir + '%s_haplotypes' % chr
	log_file = dir + '%s_haplotypes.log' % chr
	call = '/mnt/lustre/home/sonal.singhal1/bin/shapeit_v2r790/shapeit --input-ped %s %s --duohmm -W 3 --output-max %s -T 8 --output-log %s' % (ped_file, map_file, out_file, log_file)
	if chr == 'chr16':
		subprocess.call('echo \"%s\" | qsub -l h_vmem=2g -cwd -V -j y -N \"%s\"' % (call, chr + '_hap'), shell=True)
	elif chr in chrs_long:
		subprocess.call('echo \"%s\" | qsub -l h_vmem=40g -cwd -V -j y -N \"%s\"' % (call, chr + '_hap'), shell=True)
	else:
		subprocess.call('echo \"%s\" | qsub -l h_vmem=20g -cwd -V -j y -N \"%s\"' % (call, chr + '_hap'), shell=True)
