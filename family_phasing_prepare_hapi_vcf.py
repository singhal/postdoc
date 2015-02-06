import re
import argparse
import gzip

parser = argparse.ArgumentParser(description='Make trimmed VCF and reference files for family phasing in SHAPEit from HAPI output.')
parser.add_argument('--chr', help='chromosome')
args = parser.parse_args()
chr = args.chr

unrel_vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/for_shapeit/gatk.ug.finch19.%s.allfilters.recoded_biallelicSNPs.vcf.gz'  % chr
locs_file = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/input_files/%s.hapi.noswitch.locs' % chr
# note will only use parents 4 haplotypes because the kids haplotypes are related
hapi_file = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/output_files/%s.noswitch.csv' % chr

# reference files for shapeit
hap_ref = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/for_shapeit/%s.hap.gz' % chr
leg_ref = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/for_shapeit/%s.legend.gz' % chr
samp_ref = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/all_vcf/for_shapeit/%s.sample' % chr

fam_sites = {}
f = open(locs_file, 'r')
for l in f:
	d = re.split('\s+', l)
	fam_sites[ int(d[2]) ] = 1
f.close()

snps = {}
v = gzip.open(unrel_vcf, 'r')
out = unrel_vcf.replace('.vcf.gz', '.nomissing.vcf.gz')
o = gzip.open(out, 'w')
for l in v:
	if not re.search('^#', l):
                d = re.split('\t', l.rstrip())
		if int(d[1]) in fam_sites:
			num_missing = len(re.findall('\.\/\.', l))
			if num_missing < len(d[9:]):
				snps[ int(d[1]) ] = {'ref': d[3], 'alt': d[4]}
				o.write(l)
	else:
		o.write(l)
v.close()
o.close()

ordered_fam_sites = sorted(fam_sites.keys())
haps = {}
h_f = open(hapi_file, 'r')
for l in h_f:
	if re.search('^1', l):
		d = re.split(',', l.rstrip())
		d = d[1:]
		d = d[:-1]		

		for pos, site in zip(ordered_fam_sites, d):
			if pos in snps:
				if pos not in haps:
					haps[pos] = []
				site = int(site) - 1
				if site == -1:
					site = '?'
				haps[pos].append(site)
h_f.close()

samp_r = open(samp_ref, 'w')
samp_r.write('sample population group sex\n')
for ix in range(2):
	samp_r.write('hap%s hap%s hap%s 1\n' % (ix, ix, ix))
samp_r.close()

leg_r = gzip.open(leg_ref, 'w')
leg_r.write('id position a0 a1\n')
for ix, pos in enumerate(sorted(snps.keys())):
	if pos in fam_sites:
		leg_r.write('SNP%s %s %s %s\n' % (ix, pos, snps[pos]['ref'], snps[pos]['alt']))
leg_r.close()

hap_r = gzip.open(hap_ref, 'w')
for pos in sorted(haps.keys()):
	hap_r.write(' '.join([str(x) for x in haps[pos]]) + '\n')
hap_r.close()
