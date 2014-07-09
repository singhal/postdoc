import re
import sys
import gzip


def average_depth(cov_summary):
	'''
	Get average depth for the entire genome.
	'''

	avg_depth = 0
	r = open(cov_summary, 'r')
	for line in r:
		if re.search("average_depth:\s+([0-9|\.]+)", line): 
			avg_depth = re.search("average_depth:\s+([0-9|\.]+)", line).group(1)
			avg_depth = float(avg_depth)
	r.close()
	
	return avg_depth


def get_mendel(genome, mendel_file):
	'''
	Get positions that are in Mendelian error.
	'''

	me = dict()
	for chr in genome:
		me[chr] = dict() 
	if mendel_file:
		r = open(mendel_file, 'r')
		for line in r:
			if re.search('^(\S+)\s+(\S+)', line):
				match = re.search('^(\S+)\s+(\S+)', line)
				me[ match.group(1) ][int( match.group(2) )] = 1
		r.close()
	return me


def get_variant(genome, vcf_file):
	'''
	Get the variant info from VCF files.	
	'''

	var = dict()
	for chr in genome:
		var[chr] = dict()

	r = gzip.open(vcf_file, 'r')
	cur_chr = 'chr10'

	for line in r:
		# ditch the awkward header rows
		if not re.search('^#', line):
			if re.search('AF=', line):
				d = re.split("\t",line)			
				if d[0] in genome:
					if d[0] != cur_chr:
						print "var %s" % cur_chr
					cur_chr = d[0]

					if d[6] == 'PASS':
						var[d[0]][int(d[1])] = 'PASS'
					else:
						var[d[0]][int(d[1])] = 'LOWQUAL'
	r.close()
	return var


def ditched_sites(genome, cov_data, avg_depth):
	'''
	Get depth file and figure out ditched
	Uh oh! this assumes that the file starts with 'chr1'
		need to make it smarter
	'''

	ditch = dict()
	for chr in genome:
		ditch[chr] = dict()
	r = open(cov_data, 'r')
	cur_pos = 0
	cur_chr = 'chr10'
	for line in r:
		d = re.split("\t", line)
	
		if d[0] in genome:
			d[1] = int(d[1])
			d[2] = float(d[2])

			# this super kludgy section accounts for the small portion 
			#	of sites that never got called for depth ... maybe 
			#	because never a single read got mapped?
			if d[0] == cur_chr:
				if d[1] - cur_pos > 1:
					for i in range(cur_pos + 1, d[1]):
						ditch[d[0]][i] = 1
			else:
				print "ditch %s" % cur_chr

			cur_pos = d[1] 
			cur_chr = d[0]
	
			if (d[2] < 0.5 * avg_depth) or (d[2] > 2 *avg_depth):
				ditch[d[0]][d[1]] = 1
	r.close()
	return ditch


def write_to_file(genome, chrs, out_file, me, var, ditch):
	'''
	Write genome and genome summary data to file.
	Genome data are written as 60chr to line.
	'''

	out = open(out_file, 'w')
	# will want to report on the types of sites later!
	site_types = dict()
	for i in range(8):
		site_types[i] = 0

	# iterate through each chromosome one-by-one
	for chr in chrs:
		# get the chr size from the array
		chr_size = genome[chr]

		out.write('>%s\n' % chr)

		counter = 0
		for i in xrange(1, chr_size+1):
			id = 'N'
			# check if outside coverage limits
			if i in ditch[chr]:
				# a variable site
                        	if i in var[chr]:
                                	# passing site
                                	if var[chr][i] == 'PASS':               
                                        	# check for Mendelian error
                                        	if i in me[chr]:
                                      	        	# this is a bad cov, passing site with me error
                                                	id = 6                           
                                        	else:
                                                	# this is a bad cov, passing site w/o me error
                                                	id = 5                           
                                	else:
                                        	# this is a bad cov, lowqual site
                                        	id = 7
                        	# not variable site
                        	else:
                                	# this is a bad coverage, not variable site
                                	id = 4
			# coverage is okay!
			else:
				# a variable site
				if i in var[chr]:
					# passing site
					if var[chr][i] == 'PASS':
						# check for Mendelian error
						if i in me[chr]:
							# this is a good cov, passing site with me error
							id = 2
						else:
							# this is a good cov, passing site w/o me error
							id = 1
					else:
						# this is a good cov, lowqual site
						id = 3
				# not variable site
				else:
					# this is a good coverage, not variable site
					id = 0
		
			out.write('%s' % id)
			site_types[id] += 1
			counter += 1
			# have written out 60chr so time to start a new line
			if counter > 59:
				out.write('\n')
				counter = 0
		out.write('\n')
	out.close()
	
	# output a summary file with info on how many of each site I have
	summary_file = out_file + '.summary.out'
	s_out = open(summary_file, 'w')
	for i in site_types.keys():
		s_out.write('%s: %s\n' % (i, site_types[i]))
	s_out.close()
	
	return


def main():	
	'''
	Running all the subroutines.
	Currently requires editing of script for different runs ... will port to ArgParser if necessary.
	Note that this script is slow and requires a lot of memory ~32 g.
	'''
	
	# ideally, vcf file with variants only for time efficiency
	# can use vcf file with non-variants, but will be less efficient
	# this file has been through vqsr but has not yet been filtered
	# assumes this is gzipped
	vcf_file = '/mnt/gluster/home/emleffler/genotype_callsets/zebrafinch/zf_unrels/unified_genotyper/after_vqsr/gatk.ug.unrelzf.allchrs.snps.indels.vqsr2.vcf.gz' 
 
	# file with coverage summary, created by 'get_average_depth.pl'
	cov_summary = '/mnt/lustre/home/sonal.singhal1/ZF/masked_genome/zebrafinch_depth_summary.txt'
	# file with average coverage data, created by 'get_average_depth.pl'
	cov_data = '/mnt/lustre/home/sonal.singhal1/ZF/masked_genome/zebrafinch_avg_depth.txt' 

	# mendelian error file
	# if you don't have any info on mendelian errors, assign to ''
	# assumed to be in PLINK output format of TDT of chromosome and site w. header row
	mendel_file = '/mnt/gluster/home/emleffler/genotype_callsets/zebrafinch/zf_family/unified_genotyper/after_vqsr/mendelian_errors/gatk.ug.MP1-5.allchrs.vqsr.filtered.allsites.with.errors.positions'

	# out genome file -- beware that this will be big
	# we will automatically create a summary file for this file in the same directory
	out_file = '/mnt/lustre/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.fa'

	# sizes of chromosomes as used for all experiments 
	# equivalent to the main mapping genome used in /KG/sannareddy/reference/
	# note that all unknown and random chromosomes dropped -- no variants were called 
	#       on these chromosomes
	chrs = [ 'chr1', 'chr1A', 'chr1B', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
       	         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
                 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
                 'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ', 'chrM']
	genome = {'chr1': 118548696, 'chr1A': 73657157, 'chr1B': 1083483, 'chr2': 156412533,
                  'chr3': 112617285, 'chr4': 69780378, 'chr4A': 20704505, 'chr5': 62374962,
                  'chr6': 36305782, 'chr7': 39844632, 'chr8': 27993427,'chr9': 27241186,
                  'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                  'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                  'chr18': 11201131, 'chr19': 11587733, 'chr20': 15652063, 'chr21': 5979137,
                  'chr22': 3370227, 'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379,
                  'chr26': 4907541, 'chr27': 4618897, 'chr28': 4963201, 'chrLG2': 109741,
                  'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351, 'chrM': 16853}
		
	avg_depth = average_depth(cov_summary)
	me = get_mendel(genome, mendel_file)
	var = get_variant(genome, vcf_file)
	ditch = ditched_sites(genome, cov_data, avg_depth)
	write_to_file(genome, chrs, out_file, me, var, ditch)


if __name__ == "__main__":
	main()


