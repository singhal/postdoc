import re
import glob
import subprocess
import os

dir = '/mnt/gluster/home/sonal.singhal1/simulations/hotspot_simulations/sim_10_20_40_60_80_100/'
out_dir = '/mnt/gluster/home/sonal.singhal1/simulations/hotspot_simulations/sim_bpen/'
fasta_files = glob.glob('%shaplo/*fa' % dir)

for bpen in [5, 10, 50, 100, 500]:
	for ix, fasta_file in enumerate(fasta_files):
		ancestral = fasta_file.replace('haplo/', 'ancallele/')
		ancestral = ancestral.replace('haplo', 'ancallele')
		ancestral = ancestral.replace('.fa', '.txt')

		sim_num = int(re.search('(\d+).fa', fasta_file).group(1))

		out = fasta_file.replace('haplo/', 'maps/')
		out = out.replace('haplo', 'recombination')
		out = out.replace('.fa', '')
		out = out.replace(dir, out_dir)
		out = '%s_bpen%s' % (out, bpen)		

		# theta = re.search('theta([0-9|\.]+)', fasta_file).group(1)

		# lk = '%ssimulation_theta%s.lk' % (dir, theta)
		# pade = '%ssimulation_theta%s.pade' % (dir, theta)
		
		lk = '%ssimulation.lk' % dir
		pade = '%ssimulation.pade' % dir

		if not os.path.isfile(out + '.txt'):
			if sim_num < 12:
				call = '/mnt/lustre/home/sonal.singhal1/bin/LDhelmet_v1.6/ldhelmet rjmcmc --num_threads 12 -o %s -n 1000000 --burn_in 100000 -b %s -s %s -l %s -p %s -a %s -m /mnt/gluster/home/sonal.singhal1/ZF/analysis/mutation_matrix.txt -w 50' % (out, bpen, fasta_file, lk, pade, ancestral)
				print call
			#subprocess.call('echo \"%s\" | qsub -l h_vmem=15g -cwd -V -j y -N s%s_%s' % (call, bpen, ix), shell=True)
