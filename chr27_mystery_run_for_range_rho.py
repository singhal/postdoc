import subprocess

rho = [0.000005, 0.00001, 0.00005, 0.0001, 0.0005,  0.001, 0.005, 0.01, 0.05, 0.1, 0.5]

for ix, rho_val in enumerate(rho):
	out1 = 'chr27_haplotypes%s_subset2_now'  % ix
	out2 = 'chr27_haplotypes%s_subset2_w' % ix

	call1 = '/mnt/lustre/home/sonal.singhal1/bin/shapeit_v2r790/shapeit -assemble --input-vcf /mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.chr27.coverage.filtered.recoded_biallelicSNPs.nomendel.vcf_subset2.gz --input-pir /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/chr27_PIRlist -O /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s -L /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s.log --thread 8 --output-graph /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s.graph --rho %s' % (out1, out1, out1, rho_val)

	call2 = '/mnt/lustre/home/sonal.singhal1/bin/shapeit_v2r790/shapeit -assemble --input-vcf /mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.chr27.coverage.filtered.recoded_biallelicSNPs.nomendel.vcf.gz --input-pir /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/chr27_PIRlist -O /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s -L /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s.log --thread 8 --output-graph /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s.graph --rho %s --window 0.5' % (out2, out2, out2, rho_val)

	subprocess.call('echo \"%s\" | qsub -l h_vmem=2g -cwd -V -j y -N %s' % (call1, out1), shell=True)
	subprocess.call('echo \"%s\" | qsub -l h_vmem=2g -cwd -V -j y -N %s' % (call2, out2), shell=True)

