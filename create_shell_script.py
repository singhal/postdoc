
chrs = [ 'chr1', 'chr1A', 'chr1B', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
         'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
         'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ' ]

for chr in reversed(chrs):
	print 'echo "python ~/scripts/make_haplotypes.py --chr \'%s\'" | qsub -l h_vmem=20g -cwd -V -j y -N \'ht_%s\'' % (chr, chr)
	#call = '/mnt/lustre/home/sonal.singhal1/bin/LDhelmet_v1.6/ldhelmet rjmcmc --num_threads 12 -o /mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/%s.bpen10 -n 1000000 --burn_in 100000 -b 10 -s /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s_haplotypes.fasta -l /mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/ZF.lk -p /mnt/gluster/home/sonal.singhal1/ZF/analysis/LDhelmet/ZF.pade -a /mnt/gluster/home/sonal.singhal1/ZF/ancestral_allele/ancestral_allele.%s.ldhelmet.txt -m /mnt/gluster/home/sonal.singhal1/ZF/analysis/mutation_matrix.txt -w 50' % (chr, chr, chr)
	#print 'echo \"%s\" | qsub -l h_vmem=20g -cwd -V -j y -N \'ld_%s\'' % (call, chr)
