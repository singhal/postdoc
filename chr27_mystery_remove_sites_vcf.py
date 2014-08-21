import gzip
import re
import subprocess

start = 138840
end = 139830
vcf  = 'gatk.ug.unrel_zf.chr27.coverage.filtered.recoded_biallelicSNPs.nomendel.vcf.gz'

header = []
snps = []
v_open = gzip.open(vcf, 'r')
for l in v_open:
	l = l.rstrip()
	if re.match('^#', l):
		header.append(l)
	else:
		snps.append(l)
v_open.close() 
	
out_file1 = vcf.replace('vcf', 'vcf_subset1')
out_file2 = vcf.replace('vcf', 'vcf_subset2')
        
out_f1 = gzip.open(out_file1, 'w')
out_f2 = gzip.open(out_file2, 'w')
for l in header:
	out_f1.write('%s\n' % l)
	out_f2.write('%s\n' % l)

for ix, l in enumerate(snps):
	if ix < start:
		out_f1.write('%s\n' % l)
	else:
		if ix > end:
			out_f1.write('%s\n' % l)
		else:
			out_f2.write('%s\n' % l)
out_f1.close()
out_f2.close()

#name = 'chr27_haplotypes_subset1'
#call = '/mnt/lustre/home/sonal.singhal1/bin/shapeit_v2r790/shapeit -assemble --input-vcf %s --input-pir /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/chr27_PIRlist -O /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s -L /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s.log --window 0.5 --thread 8 --rho 0.0008 --output-graph /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s.graph' % (out_file1, name, name, name)
#subprocess.call('echo \"%s\" | qsub -l h_vmem=2g -cwd -V -j y -N %s' % (call, name), shell=True)
