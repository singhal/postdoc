import re
import argparse
import gzip

parser = argparse.ArgumentParser(description='Make trimmed VCF and reference files for family phasing in SHAPEit from HAPI output.')
parser.add_argument('--chr', help='chromosome')
args = parser.parse_args()
chr = args.chr

fam_vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/gatk.ug.all_zf.%s.coverage.filtered.repeatmasked.vqsr.vcf.gz' % chr
unrel_vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/for_shapeit/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.recoded_biallelicSNPs.nomendel.vcf.gz'  % chr
error_file = '/mnt/gluster/home/sonal.singhal1/ZF/mendelian_errors/plink_results/all_zf.me.%s.mendel' % chr
sites_file = '/mnt/gluster/home/sonal.singhal1/ZF/mendelian_errors/all_zf.%s.map' % chr
# note will only use parents 4 haplotypes because the kids haplotypes are related
hapi_file = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/hapi.%s.csv' % chr

# reference files for shapeit
hap_ref = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/for_shapeit/%s.hap.gz' % chr
leg_ref = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/for_shapeit/%s.legend.gz' % chr
samp_ref = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/for_shapeit/%s.sample' % chr

sites_f = open(sites_file, 'r')
sites = {}
for l in sites_f:
        d = re.split('\s+', l)
        sites[d[1]] = int(d[3])
sites_f.close()
        
errors = {}
error_f = open(error_file, 'r')
junk = error_f.next()
for l in error_f:
        d = re.split('\s+', l)
        errors[sites[d[4]]] = 1
error_f.close()
sites = {}

fam_sites = {}
v = gzip.open(fam_vcf, 'r')
for l in v:
        if re.search('^#CHROM', l):
                d = re.split('\t', l.rstrip())
                inds = d[-5:]
        if not re.search('^#', l):
                d = re.split('\t', l.rstrip())
                alleles = [d[3]] + re.split(',', d[4])
                indel = False
                for allele in alleles:
                        if len(allele) > 1:
                                indel = True
                if not indel and len(alleles) == 2:
                        if int(d[1]) not in errors:
                                all_missing = True
                                for geno in d[-5:]:
                                        geno = re.search('^(\S/\S)', geno).group(1)
                                        if geno != './.':
                                                all_missing = False
                                if not all_missing:
					fam_sites[ int(d[1]) ] = 1
v.close()

snps = {}
v = gzip.open(unrel_vcf, 'r')
out = unrel_vcf.replace('.vcf.gz', '.trimmed.vcf.gz')
o = gzip.open(out, 'w')
for l in v:
	if not re.search('^#', l):
                d = re.split('\t', l.rstrip())
		if int(d[1]) in fam_sites:
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
		
