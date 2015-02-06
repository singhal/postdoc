import gzip
import re
import subprocess

chr = 'chr2'
chop = 10
start = 1
end = 15641252
vcf  = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/for_shapeit/gatk.ug.finch21.%s.allfilters.recoded_biallelicSNPs.vcf.gz' % chr

chop_n = int((end - start) / float(chop))
	
for ix, file_start in enumerate(range(start, end, chop_n)):
	file_end = file_start + chop_n
	if file_end > end:
		file_end = end
	out_file = vcf.replace('vcf.gz', 'vcf%s.gz' % ix)
	out_f = gzip.open(out_file, 'w')        

	v_open = gzip.open(vcf, 'r')
	for l in v_open:
        	l = l.rstrip()
        	if re.match('^#', l):
                	out_f.write('%s\n' % l)
		else:
			d = re.split('\t', l)
			pos = int(d[1])
			if pos >= file_start and pos <= file_end:	
				out_f.write('%s\n' % l)
	out_f.close()
	v_open.close()

	name = '%s_haplotypes_s%s_e%s' % (chr, file_start, file_end)
	call = '/mnt/lustre/home/sonal.singhal1/bin/shapeit_v2r790/shapeit -assemble --input-vcf %s --input-pir /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/finch21/%s_PIRlist -O /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/finch21/%s -L /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/finch21/%s.log --window 0.5 --thread 8 --rho 0.0008' % (out_file, chr, name, name)

	print 'echo \"%s\" | qsub -l h_vmem=2g -cwd -V -j y -N %s' % (call, name)
