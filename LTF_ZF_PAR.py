import glob

ids = ['101','105','109','113','117','121','129','133','137','141','145','149','153','161','165','173','177','185','189']

ids = ['/mnt/gluster/home/sonal.singhal1/ZF/bam_files/' + id + '.recal.bam' for id in ids]
id_string = '-I ' + ' -I '.join(ids)

print '/home/shyamg/bin/java Xmx4g -jar /home/sonal.singhal1/bin/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /mnt/gluster/home/sonal.singhal1/reference/Zfinch.fa %s -L chrZ_random -glm SNP -mbq 20 -hets 0.006 -out_mode EMIT_ALL_SITES -o /mnt/gluster/home/sonal.singhal1/PAR/gatk.ug.unrel_zf.chrZ_random.raw.snps.vcf -nct 4' % id_string

ids = glob.glob('/mnt/gluster/data/internal_restricted_supp/finches_2014/longtailedfinch/realigned_mateFixed_recal_bams/*bam')
id_string = '-I ' + ' -I '.join(ids)
print '/home/shyamg/bin/java Xmx4g -jar /home/sonal.singhal1/bin/GenomeAnalysisTK.jar -T UnifiedGenotyper -R /mnt/gluster/home/sonal.singhal1/reference/Zfinch.fa %s -L chrZ_random -glm SNP -mbq 20 -hets 0.006 -out_mode EMIT_ALL_SITES -o /mnt/gluster/home/sonal.singhal1/PAR/gatk.ug.ltf.chrZ_random.raw.snps.vcf -nct 4' % id_string


# for id in ids:
#	print '/home/shyamg/bin/java -Xmx4g -jar ~/bin/GenomeAnalysisTK.jar -T PrintReads -R /mnt/gluster/home/sonal.singhal1/reference/Zfinch.fa -I /scratch/sannareddyk/FinchSeq/realigned_bams/realigned_mateFixed_bams/%s.mateFixed.realigned.bam -BQSR /scratch/sannareddyk/FinchSeq/realigned_bams/realign_intervals_grps/recal_grps/%s.recal.grp -o /mnt/gluster/home/sonal.singhal1/ZF/bam_files/%s.recal.bam -nct 4' % (id, id, id)
