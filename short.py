import glob

#files = glob.glob('x*')
#files = sorted(files)

#for ix, file in enumerate(files):
#	print 'echo \"sh %s\" | qsub -l h_vmem=1g -cwd -V -j y -N \'seqldhot%s\'' % (file, ix)

chrs = ['chr1', 'chr2', 'chr3', 'chr1A', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7',
	'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15']

ix = 0
for chr in chrs:
	for sp in ['ZF', 'LTF']:
		print "echo \"python hotspots_DAF.py --chr %s --sp %s\" | qsub -l h_vmem=4g -cwd -V -j y -N 'spots%s'" % (chr, sp, ix)
		ix += 1
