import os

chrs = [ 'chr1', 'chr1A', 'chr1B', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
         'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
         'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ']

long_chrs = [ 'chr1', 'chr1A', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15']


for chr in long_chrs:
	print "echo \'python plot_recombination_af1.py --chr %s\' | qsub -l h_vmem=5g -cwd -V -j y -N af%s" % (chr, chr)

# ix = 0
# for chr in chrs:
#	for window in [50000]:
#		for sp in ['LTF', 'ZF']:
#			window = int(window)
#			out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/maps/%s.window%s.bpen100.txt' % (sp, chr, window)
#			if not os.path.isfile(out):
#				print "echo \'python smooth_maps.py --chr %s --sp %s --window %s\' | qsub -l h_vmem=5g -cwd -V -j y -N smooth%s" % (chr, sp, window, ix)
#				ix += 1

