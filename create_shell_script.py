import os

chrs = [ 'chr1', 'chr1A', 'chr1B', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
         'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
         'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ']

for chr in reversed(chrs):
	print 'echo \'python simple_ancestral_chromosomes.py --chr %s\' | qsub -l h_vmem=3g -cwd -V -j y -N %s' % (chr, chr)
	#print 'echo \'python /mnt/lustre/home/sonal.singhal1/scripts/make_ldhelmet_ancestral_allele.py --chr %s --sp ZF\' | qsub -l h_vmem=3g -cwd -V -j y -N \"aa_%s\"' % (chr, chr)
	#call = '/mnt/lustre/home/sonal.singhal1/bin/shapeit_v2r790/shapeit -assemble --input-vcf /mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/for_shapeit/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.recoded_biallelicSNPs.nomendel.trimmed.vcf.gz --input-pir /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s_PIRlist -O /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s_haplotypes -L /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s_haplotypes --window 0.5 --thread 8 --rho 0.0005 --output-graph /mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/%s_haplotypes.graph -R /mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/for_shapeit/%s.hap.gz /mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/for_shapeit/%s.legend.gz /mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/for_shapeit/%s.sample --no-mcmc' % (chr, chr, chr, chr, chr, chr, chr, chr)
	#print "echo \"%s\" | qsub -l h_vmem=10g -cwd -V -j y -N \"PIR2_%s\"" % (call, chr)
	#for sp in ['ZF', 'LTF']:
	#	out = '/mnt/gluster/home/sonal.singhal1/gene_trees/%s_%s_haplotypes.fasta' % (sp, chr)
	#	if not os.path.isfile(out):
	#		print "echo \"python ~/scripts/make_haplotypes.py --sp %s --chr %s\" | qsub -l h_vmem=20g -cwd -V -j y -N \"%s_%s\"" % (sp,chr,sp,chr)
