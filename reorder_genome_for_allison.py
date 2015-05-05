import subprocess

# you need to have samtools installed and in your path

# original genome file
# you need to change this to your genome file, including path
original_file = '/mnt/gluster/home/sonal.singhal1/reference/Zfinch.fa'

# new genome file
# you need to change this to your desired genome file name, including path
new_file = '//mnt/gluster/home/sonal.singhal1/reference/test.fa'

ordered_contigs = ['chr1', 'chr1A', 'chr1B', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28', 'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ', 'chrM', 'chrUn', 'chr1A_random', 'chr1B_random', 'chr1_random', 'chr2_random', 'chr3_random', 'chr4A_random', 'chr4_random', 'chr5_random', 'chr6_random', 'chr7_random', 'chr8_random', 'chr9_random', 'chr10_random', 'chr11_random', 'chr12_random', 'chr13_random', 'chr14_random', 'chr15_random', 'chr16_random', 'chr17_random', 'chr18_random', 'chr19_random', 'chr20_random', 'chr21_random', 'chr22_random', 'chr23_random', 'chr24_random', 'chr25_random', 'chr26_random', 'chr27_random', 'chr28_random', 'chrLGE22_random', 'chrZ_random']

for contig in ordered_contigs:
	subprocess.call("samtools faidx %s %s >> %s" % (original_file, contig, new_file), shell=True)
