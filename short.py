import glob

files = glob.glob('x**')
files = sorted(files)

for ix, file in enumerate(files):
	print 'echo \"sh %s\" | qsub -l h_vmem=500M -cwd -V -j y -N \'ph%s\'' % (file, ix)
