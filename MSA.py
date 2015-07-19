import re
import glob
import subprocess
import os

chrs = ['chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6',
	'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15']
chr_lengths = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrZ': 72861351}


dir = '/mnt/gluster/home/sonal.singhal1/MSA/'
zf_genome = '/mnt/gluster/home/sonal.singhal1/reference/ZF_reference.fa'
fic_genome = '/mnt/gluster/home/sonal.singhal1/ficedula/Ficedula_albicollis.FicAlb15.fa'
geo_genome = '/mnt/gluster/home/sonal.singhal1/Darwin/g_fortis/geoFor1.masked.fa'
blast = '/mnt/gluster/home/sonal.singhal1/Darwin/g_fortis/gfortis.zebrafinch.blastn.out'

def make_dirs():
	for chr in chrs:
		tmp = dir + chr
		if not os.path.exists(tmp):
			os.makedirs(tmp)

def put_simple(genome, name):
	for chr in chrs:
		out = '%s%s/%s' % (dir, chr, name)
		o = open(out, 'w')
		seqcall = subprocess.Popen('samtools faidx %s %s' % (genome, chr), shell=True, stdout=subprocess.PIPE)
		for l in seqcall.stdout:
			if re.search('>', l):
				o.write('>%s\n' %  name)
			else:
				o.write('%s\n' % l.rstrip())
		o.close()

def put_geospiza(genome, blast):
	b = open(blast, 'r')
	matches = {}
	for l in b:
		d = re.split('\t', l)
		if d[0] not in matches:
			matches[d[0]] = d[1]
	b.close()

	match2 = {}
	for c, chr in matches.items():
		if chr not in match2:
			match2[chr] = []
		match2[chr].append(c)

	for chr in chrs:
		out = '%s%s/geospiza' % (dir, chr)
                o = open(out, 'w')
		for c_num, c in enumerate(match2[chr]):
			seqcall = subprocess.Popen('samtools faidx %s %s' % (genome, c), shell=True, stdout=subprocess.PIPE)
                	for l in seqcall.stdout:
                	        if re.search('>', l):
                        	        o.write('>geospiza:contig%s:1:+:%s\n' % (c_num, chr_lengths[chr]))
                        	else:
                                	o.write('%s\n' % l.rstrip())
		o.close()

def generate_shell_script():
	for chr in chrs:
		out = '%s%s/%s_tba.sh' % (dir, chr, chr)
		subprocess.call('tba "((zebrafinch geospiza) ficedula)" %s%s/*.*.maf %s%s/%s.maf >& %s' % (dir, chr, dir, chr, chr, out), shell=True)
		#subprocess.call('all_bz - "((zebrafinch geospiza) ficedula)" >& %s' % out, shell=True)
		#print 'echo \"sh %s\" | qsub -l h_vmem=5g -cwd -V -j y -N %s' % (out, chr)


# make_dirs()
# put_simple(zf_genome, 'zebrafinch')
# put_simple(fic_genome, 'ficedula')
# put_geospiza(geo_genome, blast)
generate_shell_script()
