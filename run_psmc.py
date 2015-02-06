import glob
import subprocess

files = glob.glob('/mnt/gluster/home/sonal.singhal1/*/analysis/PSMC/*chrZ*psmcfa')
for ix, file in enumerate(files):
	# out = file.replace('.fasta', '.30.psmcfa')
	# call = '~/bin/psmc/utils/fq2psmcfa -s 30 %s > %s' % (file, out)
	# subprocess.call('echo \'%s\' | qsub -l h_vmem=1g -cwd -V -j y -N \'psmc%s\'' % (call, ix), shell=True)

	out = file.replace('.psmcfa', '.psmc')
	call = '~/bin/psmc/psmc -N25 -t15 -r5 -p \"4+25*2+4+6\" -o %s %s' % (out, file)
	print call
	#subprocess.call('echo \'%s\' | qsub -l h_vmem=20g -cwd -V -j y -N \'psmc%s\'' % (call, ix), shell=True)

