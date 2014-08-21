import gzip
import re
import subprocess

chop = 3
start = 138845
end = 140810
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

chop_n = int((end - start) / float(chop))
	
for ix, file_start in enumerate(range(start, end, chop_n)):
	file_end = file_start + chop_n
	if file_end > end:
		file_end = end
	out_file = vcf.replace('vcf', 'vcf%s' % ix)
        
        out_f = gzip.open(out_file, 'w')
        for l in header:
                out_f.write('%s\n' % l)
	for l in snps[file_start:file_end]:
		out_f.write('%s\n' % l)
	out_f.close()

	name = 'chr27_haplotypes_s%s_e%s' % (file_start, file_end)
	call = '/mnt/lustre/home/sonal.singhal1/bin/shapeit_v2r790/shapeit -assemble --input-vcf %s --input-pir /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/chr27_PIRlist -O /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s -L /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s.log --window 0.5 --thread 8 --rho 0.0008 --output-graph /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s.graph' % (out_file, name, name, name)

	subprocess.call('echo \"%s\" | qsub -l h_vmem=2g -cwd -V -j y -N %s' % (call, name), shell=True)
