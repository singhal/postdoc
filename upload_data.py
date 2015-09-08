import glob

files = glob.glob('*bam')

for ix, file in enumerate(files[1:]):
	print "/data/tools/Aspera/bin/ascp -QT -l300M -L- %s Webin-41456@webin.ebi.ac.uk:." % file
