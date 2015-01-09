import gzip
import re
import subprocess

chr = 'chr1'
chop = 1
start = 54374254
end = 54690383
vcf  = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/gatk.ug.all_zf.%s.coverage.filtered.repeatmasked.recoded_biallelicSNPs.nomendel.vcf.gz' % chr

chop_n = int((end - start) / float(chop))
	
for ix, file_start in enumerate(range(start, end, chop_n)):
	file_end = file_start + chop_n
	if file_end > end:
		file_end = end
	out_file = vcf.replace('vcf', 'vcf%s' % ix)
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
	call = '/mnt/lustre/home/sonal.singhal1/bin/shapeit_v2r790/shapeit -assemble --input-vcf %s --input-pir /mnt/gluster/home/sonal.singhal1/ZF/phasing/phasing_uncertainty/switch_error/%s_PIRlist -O /mnt/gluster/home/sonal.singhal1/ZF/phasing/phasing_uncertainty/switch_error/%s -L /mnt/gluster/home/sonal.singhal1/ZF/phasing/phasing_uncertainty/switch_error/%s.log --window 0.5 --thread 8 --rho 0.0008' % (out_file, chr, name, name)

	print 'echo \"%s\" | qsub -l h_vmem=5g -cwd -V -j y -N %s' % (call, name)
