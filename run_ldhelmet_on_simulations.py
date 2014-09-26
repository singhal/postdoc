import re
import glob
import subprocess

dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspot_power/'
fasta_files = glob.glob('%s*fa' % dir)
lk = '%ssimulation.lk' % dir
pade = '%ssimulation.pade' % dir

for bpen in [5, 10, 50, 100, 500]:
	for ix, fasta_file in enumerate(fasta_files):
		ancestral = fasta_file.replace('haplo', 'ancallele')
		ancestral = ancestral.replace('.fa', '.txt')
		out = fasta_file.replace('haplo', 'recombination')
		out = out.replace('.fa', '')
		
		call = '/mnt/lustre/home/sonal.singhal1/bin/LDhelmet_v1.6/ldhelmet rjmcmc --num_threads 12 -o %s_%s -n 1000000 --burn_in 100000 -b %s -s %s -l %s -p %s -a %s -m /mnt/gluster/home/sonal.singhal1/ZF/analysis/mutation_matrix.txt -w 50' % (out, bpen, bpen, fasta_file, lk, pade, ancestral)
		#print call

		subprocess.call('echo \"%s\" | qsub -l h_vmem=15g -cwd -V -j y -N s%s_%s' % (call, bpen, ix), shell=True)
