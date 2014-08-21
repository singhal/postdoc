import glob
import gzip
import re

# files = ['SRR448675', 'SRR448683', 'SRR448686', 'SRR448687']
files = ['SRR448686']
dir = '/mnt/lustre/home/sonal.singhal1/Darwin/g_fortis/'

for name in files:
	file = dir + name + '.fastq.gz'
	outfile1 = file.replace('.fastq', '_1.fastq')
	outfile2 = file.replace('.fastq', '_2.fastq')
	
	in_file = gzip.open(file, 'r')
	out1 = gzip.open(outfile1, 'w')
	out2 = gzip.open(outfile2, 'w')
	
	for l in in_file:
		if re.search('^[@|+]SR', l):
			read_id = re.search('^([@|+]SR\S+)', l).group(1)
			out1.write('%s_1\n' % read_id)
			out2.write('%s_2\n' % read_id)
		else:
			read = l.rstrip()
			out1.write(read[0:100] + '\n')
			out2.write(read[100:200] + '\n')
	
	out1.close()
	out2.close()
	file.close()
