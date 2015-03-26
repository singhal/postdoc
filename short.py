import glob

files = glob.glob('b**')
files = sorted(files)

for ix, file in enumerate(files):
	print 'echo \"sh %s\" | qsub -l h_vmem=12g -cwd -V -j y -N \'bsim%s\'' % (file, ix)
