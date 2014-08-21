import glob
import sys
import re
import os

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

files = glob.glob('/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/results/*LGE22*_haplotypes.haps')

for file1 in files:
	chr = re.search('chr([a-z|A-Z|0-9]+)', file1).group(1)
        file2 = file1.replace('PIR', 'family')

        pir = {}
	for ind in range(5, 43, 2):
		pir[ind] = dict()
        site_index = {}
        pir_file =  '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/results/chr%s_PIRlist' % chr
        pirf = open(pir_file, 'r')
        num_sites = int(re.search('MAP\s+(\d+)', pirf.next()).group(1))
        for i in range(num_sites):
                site = re.search('^(\d+)', pirf.next()).group(1)
                site_index[i] = site
	ind = 3
	for l in pirf:	
		if re.search('^\d+\s+\d+\s*$', l):
			ind += 2
		phase = []
                for match in re.findall('(\d+)\s+[A|T|C|G]', l):
                        phase.append(site_index[int(match)])
		pir[ind]['_'.join(phase)] = phase
			
        pirf.close()

	if os.path.isfile(file2):       
        	# these are the columns in the haplotype file that represent the
                #       individuals we want to sample
                # Do not want to sample family zf
                # Only sampling one haplotype from each individual
                inds = range(5,43)

                haps1 = {}
                # initialize the dictionary
                for ind in inds:
                        haps1[ind] = {}
                     
                f1 = open(file1, 'r')
                for l in f1:
                        d = re.split('\s+', l)
                        
                        # create the haplotype
                        for ind in inds:
				haps1[ind][int(d[2])] = d[ind]

                haps2 = {}
                for ind in inds:
                        haps2[ind] = dict()
 
                f2 = open(file2, 'r')
                for l in f2:
                        d = re.split('\s+', l)
                        
                        # create the haplotype
                        for ind in inds:
                                haps2[ind][ int(d[2]) ] = d[ind]

		for ind in inds[0::2]:
			for haplo in pir[ind]:
				hap1A = ''
				hap1B = ''
				hap2 = ''
			
				for site in pir[ind][haplo]:
					site = int(site)
					if haps1[ind][site] != haps1[ind + 1][site]:
						hap1A += haps1[ind][site]
						hap1B += haps1[ind + 1][site]
						if site in haps2[ind]:
							hap2 += haps2[ind][site]
						else:
							hap2 += '-'

				if len(hap1A) > 3:
					bk_pts1 = compare_haplo(hap1A, hap2)
					bk_pts2 = compare_haplo(hap1B, hap2)

					bk_pts = min(bk_pts1, bk_pts2)
					diff = bk_pts / float( len(hap1A) )
				
					print '%s\t%s\t%s\t%s\t%s' % (chr, haplo, ind, len(hap1A), diff)
			

