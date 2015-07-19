import re
import subprocess
import random
import numpy as np

def get_alignments(n_align, chr_lengths, genome):
	genome_length = sum(chr_lengths.values())
	chr_prop = {}
	loci = {}
	for chr, length in chr_lengths.items():
		chr_prop[ chr ] = round( n_align * length / float(genome_length))
	for chr, n_loci in chr_prop.items():
		print chr
		tot_n = 0
		while tot_n < n_loci:
			start = random.randint(0, chr_lengths[chr])
			end = start + 999
			seq = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome, chr, start, end), shell=True, stdout=subprocess.PIPE)			
			n_count = 0
			for l in seq.stdout:
				if not re.match('>', l):
					n_count += list(l.rstrip()).count('0')

			if n_count > 750:
				tot_n += 1
				if chr not in loci:
					loci[chr] = dict()
				loci[chr][start] = 1	
	return loci

def make_fasta_files(loci, names, dir):
	for chr in loci:
		for start in loci[chr]:
			end = start + 999
			out_file = '%s%s_%s_%s.fasta' % (dir, chr, start, end)
			out_f = open(out_file, 'w')
			for species in names:
				genome = '%schromosomes/%s_%s_haplotypes.fasta' % (dir, species, chr)
				for haplo in names[species]:
					seq = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome, haplo, start, end), shell=True, stdout=subprocess.PIPE) 
					for l in seq.stdout:
                                		if re.match('>', l):
                                        		out_f.write('>%s\n' % names[species][haplo])
						else:
							out_f.write(l)
			out_f.close()
	return

def main():
	chr_lengths = {	'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
			'chr14': 16419078, 'chr15': 14428146, 'chr1A': 73657157,
			'chr1': 118548696, 'chr2': 156412533, 'chr3': 112617285,
			'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
			'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186}
	'''
	names = {'LTF': {'haplo0': 'LTFh_73783a', 'haplo1': 'LTFh_73783b', 'haplo2': 'LTFh_73788a', 'haplo3': 'LTFh_73788b', 'haplo4': 'LTFh_73790a',
			'haplo5': 'LTFh_73790b', 'haplo6': 'LTFh_73900a', 'haplo7': 'LTFh_73900b', 'haplo8': 'LTFh_73903a', 'haplo9': 'LTFh_73903b',
			'haplo10': 'LTFh_73907a', 'haplo11': 'LTFh_73907b', 'haplo12': 'LTFh_73933a', 'haplo13': 'LTFh_73933b', 'haplo14': 'LTFh_73942a',
			'haplo15': 'LTFh_73942b', 'haplo16': 'LTFh_73948a', 'haplo17': 'LTFh_73948b', 'haplo18': 'LTFh_73958a', 'haplo19': 'LTFh_73958b',
			'haplo20': 'LTFa_G111a', 'haplo21': 'LTFa_G111b', 'haplo22': 'LTFa_G118a', 'haplo23': 'LTFa_G118b', 'haplo24': 'LTFa_G163a',
			'haplo25': 'LTFa_G163b', 'haplo26': 'LTFa_G169a', 'haplo27': 'LTFa_G169b', 'haplo28': 'LTFa_G183a', 'haplo29': 'LTFa_G183b',
			'haplo30': 'LTFa_G250a', 'haplo31': 'LTFa_G250b', 'haplo32': 'LTFa_G276a', 'haplo33': 'LTFa_G276b', 'haplo34': 'LTFa_G294a',
			'haplo35': 'LTFa_G294b', 'haplo36': 'LTFa_W2703a', 'haplo37': 'LTFa_W2703b', 'haplo38': 'LTFa_W2994a', 'haplo39': 'LTFa_W2994b'},
		'ZF': {'haplo0': 'ZF_26462a', 'haplo1': 'ZF_26462b', 'haplo2': 'ZF_26516a', 'haplo3': 'ZF_26516b', 'haplo4': 'ZF_26721a',
                        'haplo5': 'ZF_26721b', 'haplo6': 'ZF_26733a', 'haplo7': 'ZF_26733b', 'haplo8': 'ZF_26781a', 'haplo9': 'ZF_26781b',
                        'haplo10': 'ZF_26792a', 'haplo11': 'ZF_26792b', 'haplo12': 'ZF_26795a', 'haplo13': 'ZF_26795b', 'haplo14': 'ZF_26820a',
                        'haplo15': 'ZF_26820b', 'haplo16': 'ZF_26881a', 'haplo17': 'ZF_26881b', 'haplo18': 'ZF_26896a', 'haplo19': 'ZF_26896b',
                        'haplo20': 'ZF_28016a', 'haplo21': 'ZF_28016b', 'haplo22': 'ZF_28078a', 'haplo23': 'ZF_28078b', 'haplo24': 'ZF_28313a',
                        'haplo25': 'ZF_28313b', 'haplo26': 'ZF_28339a', 'haplo27': 'ZF_28339b', 'haplo28': 'ZF_28353a', 'haplo29': 'ZF_28353b',
                        'haplo30': 'ZF_28402a', 'haplo31': 'ZF_28402b', 'haplo32': 'ZF_28404a', 'haplo33': 'ZF_28404b', 'haplo34': 'ZF_28456a',
                        'haplo35': 'ZF_28456b', 'haplo36': 'ZF_28481a', 'haplo37': 'ZF_28481b'},
		'DBF': {'haplo0': 'DBFa', 'haplo1': 'DBFb'},
		'Geospiza': {'haplo0': 'Gfortisa', 'haplo1': 'Gfortisb'} }
	'''

	names = {'LTF': {'haplo0': 'LTFh_73783a', 'haplo15': 'LTFh_73942b', 'haplo25': 'LTFa_G163b', 'haplo39': 'LTFa_W2994b'},
                'ZF': {'haplo0': 'ZF_26462a', 'haplo10': 'ZF_26792b', 'haplo25': 'ZF_28313b', 'haplo30': 'ZF_28402a'},
                'DBF': {'haplo0': 'DBFa', 'haplo1': 'DBFb'},
                'Geospiza': {'haplo0': 'Gfortisa', 'haplo1': 'Gfortisb'},
		'Ficedula': {'haplo0': 'Ficedulaa', 'haplo1': 'Ficedulab'}}	
		
	dir = '/mnt/gluster/home/sonal.singhal1/gene_trees/'
	genome = '/mnt/gluster/home/sonal.singhal1/reference/all.masked_genome.fa'
	loci = get_alignments(1000, chr_lengths, genome) 
	make_fasta_files(loci, names, dir)

if __name__ == "__main__":
    main()

