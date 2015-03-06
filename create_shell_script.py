import os

chrs = [ 'chr1', 'chr1A', 'chr1B', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
         'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
         'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ']

long_chrs = [ 'chr1', 'chr1A', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chrZ']

ix = 0
for chr in ['chr5', 'chr6']:
	for sp in ['LTF', 'ZF']:
		for block in [500, 1000, 2000]:
			for flank in [20000, 40000]:
				print 'echo \'python find_hotspots.py --sp %s --chr %s --block %s --flank %s\' | qsub -l h_vmem=3g -cwd -V -j y -N find%s' % (sp, chr, block, flank, ix)
				ix += 1

# for chr in long_chrs:
#	for sp in ['LTF', 'ZF']:
#		print 'echo \'python plot_recombination_hotspot.py --sp %s --chr %s\' | qsub -l h_vmem=2g -cwd -V -j y -N \'%s_%s_hot\'' % (sp, chr, sp, chr)

#for chr in reversed(chrs):
#	for sp in ['LTF', 'ZF']:
#		print 'echo \'python phase_vcfs_using_haps.py --sp %s --chr %s\' | qsub -l h_vmem=2g -cwd -V -j y -N %s_%s_ph' % (sp, chr, sp, chr)

#for chr in reversed(chrs):
#        print 'echo \'python 1D_sfs.py --sp LTF --chr %s\' | qsub -l h_vmem=2g -cwd -V -j y -N %s_sfs' % (chr, chr)
