import subprocess
import re

chrs = [ 'chr1', 'chr1A', 'chr1B', 'chr2', 'chr3',  'chr4', 'chr4A', 'chr5', 'chr6', 'chr7', 'chr8',
         'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18',
         'chr19', 'chr20', 'chr21', 'chr22', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28',
         'chrLG2', 'chrLG5', 'chrLGE22', 'chrZ' ]

chrs_long = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr1A']

bam_files = {   '26462': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/101.mateFixed.realigned.bam',
		'28339': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/105.mateFixed.realigned.bam',
		'28353': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/109.mateFixed.realigned.bam',
		'26721': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/113.mateFixed.realigned.bam',
		'28456': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/117.mateFixed.realigned.bam',
		'28402': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/121.mateFixed.realigned.bam',
		'26516': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/129.mateFixed.realigned.bam',
		'28404': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/133.mateFixed.realigned.bam',
		'26820': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/137.mateFixed.realigned.bam',
		'26733': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/141.mateFixed.realigned.bam',
		'28481': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/145.mateFixed.realigned.bam',
		'26881': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/149.mateFixed.realigned.bam',
		'26781': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/153.mateFixed.realigned.bam',
		'26896': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/161.mateFixed.realigned.bam',
		'26792': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/165.mateFixed.realigned.bam',
		'28016': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/173.mateFixed.realigned.bam',
		'26795': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/177.mateFixed.realigned.bam',
		'28078': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/185.mateFixed.realigned.bam',
		'28313': '/scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/189.mateFixed.realigned.bam' }

'''
bam_files = {   '73783': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/73783.realigned.mateFixed.bam',
		'73788': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/73788.recal.bam',
		'73790': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/73790.recal.bam',
		'73900': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/73900.recal.bam',
		'73903': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/73903.recal.bam',
		'73907': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/73907.recal.bam',
		'73933': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/73933.recal.bam',
		'73942': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/73942.recal.bam',
		'73948': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/73948.recal.bam',
		'73958': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/73958.recal.bam',
		'G111': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/G111.recal.bam',
		'G118': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/G118.realigned.mateFixed.bam',
		'G163': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/G163.recal.bam',
		'G169': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/G169.recal.bam',
		'G183': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/G183.recal.bam',
		'G250': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/G250.recal.bam',
		'G276': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/G276.recal.bam',
		'G294': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/G294.recal.bam',
		'W2703': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/W2703.recal.bam',
		'W2994': '/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/W2994.recal.bam'}
'''

dir = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/PIR_approach/'
vcf_dir = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/'

for chr in ['chr27']:
	bam_list = dir + chr + '_bamlist'
	o = open(bam_list, 'w')
	for bam in bam_files:
		o.write("%s\t%s\t%s\n" % (bam, bam_files[bam], chr))
	o.close()
	
	vcf_file = vcf_dir + 'gatk.ug.unrel_zf.%s.coverage.filtered.recoded_biallelicSNPs.nomendel.vcf.gz' % chr
	pir_out = dir + chr + '_PIRlist'	
	hap_out = dir + chr + '_haplotypes'

	job_file = dir + 'PIR_%s.sh' % chr
	o = open(job_file, 'w')
	o.write("~/bin/extractPIRs.v1.r68.x86_64/extractPIRs --bam %s --vcf %s --out %s" % (bam_list, vcf_file, pir_out))
	o.close()

	subprocess.call('chmod a+x %s' % job_file, shell=True)
	if chr in chrs_long:
		subprocess.call('echo "%s" | qsub -l h_vmem=10g -cwd -V -j y -N "PIR_L%s"' % (job_file, chr), shell=True)
	else:
		subprocess.call('echo "%s" | qsub -l h_vmem=3g -cwd -V -j y -N "PIR_L%s"' % (job_file, chr), shell=True)

	
	#call = '/mnt/lustre/home/sonal.singhal1/bin/shapeit_v2r790/shapeit -assemble --input-vcf %s --input-pir %s -O %s -L %s --window 0.5 --thread 8 --rho 0.0008 --output-graph %s' % (vcf_file, pir_out, hap_out, hap_out + '.log', hap_out + '.graph')
	#job_file = dir + 'PIR2_%s.sh' % chr
	#o = open(job_file, 'w')
	#o.write(call)
	#o.close()
	#subprocess.call('chmod a+x %s' % job_file, shell=True)
	#if chr in chrs_long:
	# 	subprocess.call('echo "%s" | qsub -l h_vmem=40g -cwd -V -j y -N "PIR2_%s"' % (job_file, chr), shell=True)
	#else:
	#	subprocess.call('echo "%s" | qsub -l h_vmem=20g -cwd -V -j y -N "PIR2_%s"' % (job_file, chr), shell=True)
