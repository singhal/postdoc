import re
import glob
import subprocess
import os

dir = '/mnt/gluster/home/sonal.singhal1/simulations/shared/ZF/'
fasta_files = glob.glob('%shaplo/*fa' % dir)

for bpen in [5]:
	for ix, fasta_file in enumerate(fasta_files):
		ancestral = fasta_file.replace('haplo/', 'ancallele/')
		ancestral = ancestral.replace('haplo', 'ancallele')
		ancestral = ancestral.replace('.fa', '.txt')
		out = fasta_file.replace('haplo/', 'maps/')
		out = out.replace('haplo', 'recombination')
		out = out.replace('.fa', '')
		out = '%s_%s' % (out, bpen)		

		# theta = re.search('theta([0-9|\.]+)', fasta_file).group(1)
		# lk = '%ssimulation_theta%s.lk' % (dir, theta)
		# pade = '%ssimulation_theta%s.pade' % (dir, theta)
		
		lk = '%ssimulation.lk' % dir
		pade = '%ssimulation.pade' % dir

		if not os.path.isfile(out + '.txt'):
			call = '/mnt/lustre/home/sonal.singhal1/bin/LDhelmet_v1.6/ldhelmet rjmcmc --num_threads 12 -o %s -n 1000000 --burn_in 100000 -b %s -s %s -l %s -p %s -a %s -m /mnt/gluster/home/sonal.singhal1/ZF/analysis/mutation_matrix.txt -w 50' % (out, bpen, fasta_file, lk, pade, ancestral)
			#print call
			subprocess.call('echo \"%s\" | qsub -l h_vmem=15g -cwd -V -j y -N zsim_%s' % (call, ix), shell=True)
