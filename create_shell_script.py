import os

chrs = [ 'chr1', 'chr1A', 'chr1B', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
         'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
         'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ']

long_chrs = [ 'chr1', 'chr1A', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15']

# for chr in chrs:
#	print "echo \'python calc_dnds.py --chr %s\' | qsub -l h_vmem=5g -cwd -V -j y -N '%s'" % (chr, chr)

for chr in long_chrs:
	for sp in ['ZF', 'LTF']:
		print "echo \'python get_rho_across_genes.py --sp %s --chr %s\' | qsub -l h_vmem=5g -cwd -V -j y -N %s_%s" % (sp, chr, sp, chr)


#i = 0
#species = [('ZF', 'LTF'), ('LTFa', 'LTFh')]
#for chr in chrs:
#	for (sp1, sp2) in species:
#		print "echo \'python calculate_divergence_bysnp.py --sp1 %s --sp2 %s --chr %s\' | qsub -l h_vmem=2g -cwd -V -j y -N snp%s" % (sp1, sp2, chr, i)
#		i += 1

# ix = 0
# for chr in chrs:
#	for window in [50000]:
#		for sp in ['LTF', 'ZF']:
#			window = int(window)
#			out = '/mnt/gluster/home/sonal.singhal1/%s/analysis/LDhelmet/maps/%s.window%s.bpen100.txt' % (sp, chr, window)
#			if not os.path.isfile(out):
#				print "echo \'python smooth_maps.py --chr %s --sp %s --window %s\' | qsub -l h_vmem=5g -cwd -V -j y -N smooth%s" % (chr, sp, window, ix)
#				ix += 1

