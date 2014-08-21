import re
import sys
import glob
import subprocess
import random
import gzip
import os


def get_vcf_var(vcf_file):
	hets = {}
	ids = []
	site_freq = {}
	vcf_f = gzip.open(vcf_file, 'r')
	for l in vcf_f:
		if (re.search('^#', l)):
			if (re.search('^#CHROM', l)):
				ids = re.split('\t', l.rstrip())
				ids = ids[9:len(ids)]
				for id in ids:
					hets[id] = list()
		else:
			d = re.split('\t', l.rstrip())
			for ix, id in enumerate(ids):
				if re.search('0/1', d[ix + 9]):
					hets[id].append(int(d[1]))
			allele1 = len(re.findall('0\/', l)) + len(re.findall('\/0', l))
			allele2 = len(re.findall('1\/', l)) + len(re.findall('\/1', l))
			min_allele = min(allele1, allele2)
			site_freq[int(d[1])] = min_allele
	vcf_f.close()
	return hets, ids, site_freq


def get_pir_info(pir_file, ids):
	pir = {}
	for id in ids:
		pir[id] = dict()
	pir_f = open(pir_file, 'r')
	num_sites = int(re.search('MAP\s+(\d+)', pir_f.next()).group(1))
	site_index = {}
	for i in range(num_sites):
		site = re.search('^(\d+)', pir_f.next()).group(1)
                site_index[i] = int(site)
	ind = 'NA'
	for l in pir_f:  
		if re.search('^([A-Z|0-9]+)\s+\d+\s*$', l):
			ind = re.search('^([A-Z|0-9]+)\s+\d+\s*$', l).group(1)	
		else:
			for match in re.findall('(\d+)\s+[A|T|C|G]', l):
				pir_pos = site_index[int(match)]
				if pir_pos not in pir[ind]:
					pir[ind][pir_pos] = 0
        	        	pir[ind][pir_pos] += 1
	pir_f.close()
	return pir


def summarize_uncertainty(hets, ids, site_freq, pir, chr, out_dir):
	summary = out_dir + 'phasing_uncertainty_%s.csv' % chr
	s = open(summary, 'w')
	for id in ids:
		out_file = out_dir + 'sets_%s' % (id)
		out_f = open(out_file, 'w')
		for ix, (i, j) in enumerate(zip(hets[id], hets[id][1:])):
		 	line = '%sset_%s_%s %s %s\n' % (out_dir, id, ix, i, j)
			out_f.write(line)
		out_f.close()
	
		subprocess.call("~/bin/shapeit_v2r790/shapeit -convert --input-graph %s --input-sets %s" % (graph, out_file), shell=True)
		for ix, (i, j) in enumerate(zip(hets[id], hets[id][1:])):
        		file = '%sset_%s_%s.pairs.gz'  % (out_dir, id, ix)
			o = gzip.open(file, 'r')
			values = []
			for l in o:
				if re.search('%s\s+\d+\s+\d+\s+\d+\s+([\.|0-9]+)' % ids[0], l):
					values.append(re.search('%s\s+\d+\s+\d+\s+\d+\s+([\.|0-9]+)' % ids[0], l).group(1))
			
			if int(i) in pir[id]:
				pir_i = pir[id][int(i)]
			else:
				pir_i = 0

			if int(j) in pir[id]:                                
				pir_j = pir[id][int(j)]
                        else:
				pir_j = 0


			s.write('%s,%s,%s,%s,%s,%s,%s,%s\n' % (id, i, j, max(values), site_freq[int(i)], site_freq[int(j)], pir_i, pir_j))
			o.close()
			os.remove(file)
			os.remove('%sset_%s_%s.freqs.gz'  % (out_dir, id, ix))
	s.close()


chrs = ['chrLGE22', 'chr24', 'chr10', 'chr1']
out_dir = '/mnt/gluster/home/sonal.singhal1/LTF/phasing/phasing_uncertainty/'
for chr in chrs:
	graph = '/mnt/gluster/home/sonal.singhal1/LTF/phasing/PIR_approach/%s_haplotypes.graph' % chr
	pir_file = '/mnt/gluster/home/sonal.singhal1/LTF/phasing/PIR_approach/%s_PIRlist' % chr
	vcf_file = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.filtered.coverage.recoded_biallelicSNPs.vqsr.vcf.gz' % chr
	
	hets, ids, site_freq = get_vcf_var(vcf_file)
	pir = get_pir_info(pir_file, ids)
	summarize_uncertainty(hets, ids, site_freq, pir, chr, out_dir)
