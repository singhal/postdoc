import glob
import sys
import re
import os
import gzip
import random
import subprocess
import numpy as np
import pandas as pd
import argparse

def compare_ab(a,b):
        if a != '-' and b != '-':
                if a != b:
                        return 'F'
                elif a == b:
                        return 'T'
        else:
                return '0'


def compare_haplo(hap1, hap2):
        compare = map(compare_ab, hap1, hap2)
        compare = filter(lambda x: x != '0', compare)
                  
        bk_pts = 0
        for x, y in zip(compare, compare[1:]):
                if x != y:
                        bk_pts += 1
        return bk_pts


def get_vcf_var(vcf_file):
	site_freq = {}
        vcf_f = gzip.open(vcf_file, 'r')
        
	for l in vcf_f:
                if not re.search('^#', l):
                        d = re.split('\t', l.rstrip())

                        allele1 = len(re.findall('0\/', l)) + len(re.findall('\/0', l))
                        allele2 = len(re.findall('1\/', l)) + len(re.findall('\/1', l))
                        min_allele = min(allele1, allele2)
                        site_freq[int(d[1])] = min_allele
        vcf_f.close()

        return site_freq

def get_uncertain(phasing_uncertain, max_uncertain):
	data = pd.read_csv(phasing_uncertain, sep=',', header=0, names=['ind','site1','site2','uncertainty','hets1','hets2','pir1','pir2'])

	binned_site = data.groupby('site1')
	uncertainty_vals = binned_site.aggregate(np.mean)
	snps = uncertainty_vals.index
	uncertainty = uncertainty_vals.uncertainty
	uncertain1 = snps[ uncertainty < max_uncertain]

	binned_site = data.groupby('site2')
	uncertainty_vals = binned_site.aggregate(np.mean)
	snps = uncertainty_vals.index
	uncertainty = uncertainty_vals.uncertainty
	uncertain2 = snps[ uncertainty < max_uncertain]

	uncertain = set(uncertain1.tolist() + uncertain2.tolist())
	uncertain = [int(x) for x in uncertain]	

	return uncertain

def get_sites(file):
	sites = []
	f1 = open(file, 'r')
	for l in f1:
	        d = re.split('\s+', l)
	        sites.append(d[2])
	f1.close()
	return sites


def get_haps(files, inds):
	haps = {}
	
	for ind in inds:
	        haps[ind] = list()

	for file in files:
	        f1 = open(file, 'r')
	
	        tmp_haps = {}
	        for ind in inds:
	                tmp_haps[ind] = list()
        
	        for l in f1:
	                d = re.split('\s+', l)
	                for ind in inds:
	                        tmp_haps[ind].append(d[ind])
	        f1.close()
	
	        for ind in tmp_haps:
	                haps[ind].append(tmp_haps[ind])
	return haps


def get_switch_error(out_file, snps, block_size, haps, sites, min_freq, uncertain, site_freq, inds):
	out_f = open(out_file, 'w')
	out_f.write('snp_start,switch_error\n')

	for snp_start in range(0, snps, block_size):
	
		snp_end = snp_start + block_size
		if snp_end > snps:
			snp_end = snps

		tot_bk_pts = 0
		tot_inf_sites = 0

		for ind in inds[0::2]:
			num_compare = 0

			while num_compare < 20:
				rand1 = random.randint(0,len(haps[ind]) - 1)
				rand2 = random.randint(0,len(haps[ind]) - 1)

				if rand1 != rand2:
					hap1 = ''
					hap2A = ''
					hap2B = ''
			
					for ix, (a, b, c) in enumerate(zip(haps[ind][rand1][snp_start:snp_end], \
							haps[ind][rand2][snp_start:snp_end], haps[ind+1][rand2][snp_start:snp_end])):
						pos = int(sites[snp_start + ix])
						if b != c:
							if site_freq[pos] > min_freq:
								if pos not in uncertain:
									hap1 += a
									hap2A += b
									hap2B += c
					
					if len(hap1) > 0:
						bk_pts1 = compare_haplo(hap1, hap2A)
						bk_pts2 = compare_haplo(hap1, hap2B)

						tot_bk_pts += min(bk_pts1, bk_pts2)
						tot_inf_sites += len(hap1)

					num_compare += 1
		diff = 0
		if tot_inf_sites > 1:
			diff = tot_bk_pts / float(tot_inf_sites)
		out_f.write('%s,%s,%s\n' % (sites[snp_start], sites[snp_end - 1], diff))
	out_f.close()
	return


def get_snp_count(file):
	snps = subprocess.Popen("wc %s" % file, shell=True, stdout=subprocess.PIPE)
	snps = int(re.match('^\s+(\d+)', snps.stdout.readline()).group(1))
	return snps


def sample_haplotypes(chr, graph_file, n_sam):
	out_files = []
	for i in range(n_sam):
		out_file = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/phasing_uncertainty/samples/%s.sample%s' % (chr, i)
		out_files.append(out_file + '.haps')
		if not os.path.isfile(out_file + '.haps'):
			subprocess.call("~/bin/shapeit_v2r790/shapeit -convert -T 8 --input-graph %s --output-sample %s --seed %s" % (graph_file, out_file, random.randint(0,1000)), shell=True)
	return out_files


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--chr", help="chromosome for which to run analysis")
	parser.add_argument("--samples", help="number of samples for which to run analysis")
	args = parser.parse_args()
	chr = args.chr
        n_sam = int(args.samples)

	block_size = 50
	min_freq = 1
	max_uncertainty = 0

	out_file = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/phasing_uncertainty/switch_error_rate_%s_mf%s_uncertain%s.csv' % (chr, min_freq, max_uncertainty)
	vcf_file = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.filtered.recoded_biallelicSNPs.nomendel.vcf.gz' % chr
	graph_file = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/results/%s_haplotypes.graph' % chr

	files = sample_haplotypes(chr, graph_file, n_sam)
	snps = get_snp_count(files[0])
	site_freq = get_vcf_var(vcf_file)
	# uncertain = get_uncertain(phasing_file, max_uncertainty)
	uncertain = []
	sites = get_sites(files[0])
	inds = range(5,43)
	haps = get_haps(files, inds)
	get_switch_error(out_file, snps, block_size, haps, sites, min_freq, uncertain, site_freq, inds)

if __name__ == "__main__":
    main()
