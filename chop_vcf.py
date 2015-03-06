import re
import gzip
import glob

files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/gatk.ug.all_zf.*.coverage.repeatmasked.filtered.nomendel.shared.noswitch.vqsr2.vcf.gz')
for file in ['gatk.ug.all_zf.chrZ.coverage.repeatmasked.filtered.nomendel.shared.recodedsex.vqsr2.vcf.gz']:
	out = file.replace('all_zf', 'unrel_zf')
	f = gzip.open(file, 'r')
	o = gzip.open(out, 'w')
	for l in f:
		if re.search('^#CHROM', l) or re.search('^chr', l):
			d = re.split('\t', l.rstrip())
			o.write('\t'.join(d[:-5]) + '\n')
		else:
			o.write(l)
	o.close()
	f.close()
