import subprocess

# you need to have samtools installed and in your path

# original genome file
# you need to change this to your genome file, including path
original_file = '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.repeat_masked.switch_masked.fa'

# new genome file
# you need to change this to your desired genome file name, including path
new_file = '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.repeat_masked.switch_masked.fa2'

ordered_contigs = ['chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28', 'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ']

for contig in ordered_contigs:
	subprocess.call("samtools faidx %s %s >> %s" % (original_file, contig, new_file), shell=True)
