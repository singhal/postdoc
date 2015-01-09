import glob

files = glob.glob('x*')

for ix, file in enumerate(files):
	print 'echo \"sh %s\" | qsub -l h_vmem=10G -cwd -V -j y -N \'find%s\'' % (file, ix)
