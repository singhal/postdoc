import glob

files = glob.glob('la*')

for ix, file in enumerate(files):
	print 'echo \"sh %s\" | qsub -l h_vmem=5g -cwd -V -j y -N \'lfind%s\'' % (file, ix)
