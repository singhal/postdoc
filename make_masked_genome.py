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


def ditched_sites(genome, cov_data, avg_depth):
	'''
	Get depth file and figure out ditched
	Uh oh! this assumes that the file starts with 'chr1'
		need to make it smarter
	'''

	ditch = dict()
	for chr in genome:
		ditch[chr] = dict()
	r = gzip.open(cov_data, 'r')
	cur_pos = 0
	cur_chr = 'chr1'
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


def repeat_sites(repeat, ditch):
	'''
	Remove repeat sites.
	'''
        f = open(repeat, 'r')
        for l in f:
                if re.search('^\s+\d+', l):
                        d = re.split('\s+', l.rstrip())
                        if d[5] in ditch:
                                start = int(d[6])
                                end = int(d[7]) + 1

                                for pos in range(start, end):
                                        ditch[d[5]][pos] = 1
        f.close()
        return ditch


def write_to_file(genome, out_file, ditch):
	'''
	Write genome and genome summary data to file.
	Genome data are written as 60chr to line.
	'''

	out = open(out_file, 'w')
	
	# iterate through each chromosome one-by-one
	for chr in genome:
		# get the chr size from the array
		chr_size = genome[chr]

		out.write('>%s\n' % chr)

		counter = 0
		for i in xrange(1, chr_size+1):
			id = 'N'
			# check if outside coverage limits or repeat
			if i in ditch[chr]:
				id = 8
			# coverage is okay and not repeat
			else:
				id = 0
		
			out.write('%s' % id)
			counter += 1
			# have written out 60chr so time to start a new line
			if counter > 59:
				out.write('\n')
				counter = 0
		out.write('\n')
	out.close()
	
	return


def main():	
	'''
	Running all the subroutines.
	Currently requires editing of script for different runs ... will port to ArgParser if necessary.
	Note that this script is slow and requires a lot of memory ~20 g.
	'''
	
	# file with coverage summary, created by 'get_average_depth.pl'
	cov_summary = '/mnt/gluster/home/sonal.singhal1/LTF/masked_genome/longtailed_depth_summary.txt'
	# file with average coverage data, created by 'get_average_depth.pl'
	cov_data = '/mnt/gluster/home/sonal.singhal1/LTF/masked_genome/longtailed_avg_depth.txt.gz' 

	# out genome file -- beware that this will be big
	# we will automatically create a summary file for this file in the same directory
	out_file = '/mnt/gluster/home/sonal.singhal1/LTF/masked_genome/LTF.masked_genome.repeat_masked.fa'

	# file with repeat data
        repeat_file = '/mnt/gluster/home/sonal.singhal1/reference/taeGut1.repeatMaskerBlast.repeatLibrary20140131.out'

	# sizes of chromosomes as used for all experiments 
	# note that all unknown and random chromosomes dropped -- no variants were called 
	#       on these chromosomes
	genome = {'chr1': 118548696, 'chr1A': 73657157, 'chr1B': 1083483, 'chr2': 156412533,
                  'chr3': 112617285, 'chr4': 69780378, 'chr4A': 20704505, 'chr5': 62374962,
                  'chr6': 36305782, 'chr7': 39844632, 'chr8': 27993427,'chr9': 27241186,
                  'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                  'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                  'chr18': 11201131, 'chr19': 11587733, 'chr20': 15652063, 'chr21': 5979137,
                  'chr22': 3370227, 'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379,
                  'chr26': 4907541, 'chr27': 4618897, 'chr28': 4963201, 'chrLG2': 109741,
                  'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351}
		
	avg_depth = average_depth(cov_summary)
	ditch = ditched_sites(genome, cov_data, avg_depth)
	ditch = repeat_sites(repeat_file, ditch)
	write_to_file(genome, out_file, ditch)


if __name__ == "__main__":
	main()
