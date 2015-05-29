import re
import gzip
import argparse
import subprocess
import datetime

def get_chromosome(ref_genome, chr):
	seqcall = subprocess.Popen('samtools faidx %s %s' % (ref_genome, chr), shell=True, stdout=subprocess.PIPE)
	seq = ''
	for l in seqcall.stdout:
		if not re.match('>', l):
			seq += l.upper().rstrip()
	return list(seq)


def get_masked(seq, genome, chr):
        seqcall = subprocess.Popen('samtools faidx %s %s' % (genome, chr), shell=True, stdout=subprocess.PIPE)
        mask = ''
        for l in seqcall.stdout:
                if not re.match('>', l):
                        mask += l.upper().rstrip()
	mask = list(mask)

	for pos, mask_pos in enumerate(mask):
		if int(mask_pos) > 3:
			seq[pos] = 'N'

        return list(seq)


def get_variants(vcf, seq):
	f = gzip.open(vcf, 'r')
	
	for l in f:
		if not re.search('#', l):
                        d = re.split('\t', l.rstrip())

                        # check not indel
                        indel = False
                        alleles = [d[3]] + re.split(',', d[4])
                        for allele in alleles:
                                if len(allele) > 1:
                                        indel = True

                        if not indel:
                                # check number of alleles
                                genos = []
                                for geno in d[9:]:
                                        geno = re.search('^([^:]+)', geno).group(1)
                                        genos += re.split('[|/]', geno)
                                genos = [x for x in genos if not re.search('\.', x)]
				if len(genos) > 0:
					for ix, allele in enumerate(alleles):
						af = genos.count(str(ix)) / len(genos)
						if af > 0.975:
							seq[int(d[1]) - 1] = allele
	f.close()
	return seq


def print_seq(sp, chr, seq):
	out = '/mnt/gluster/home/sonal.singhal1/reference/%s_%s.fa' % (sp, chr)
	o = open(out, 'w')
	o.write('>%s\n' % chr)
	for i in xrange(0, len(seq), 60):
                o.write(''.join(seq[i:i+60]) + '\n')
	o.close()


def main():
        parser = argparse.ArgumentParser()
        parser.add_argument("--chr", help="chromosome for which to run analysis")
        parser.add_argument("--sp", help="species [ZF|LTF|DBF] for which to run")

	args = parser.parse_args()
	chr = args.chr
	sp = args.sp

	ref_genome = '/mnt/gluster/home/sonal.singhal1/reference/Zfinch.fa'
	if sp == 'LTF':
        	genome = '/mnt/gluster/home/sonal.singhal1/LTF/masked_genome/LTF.masked_genome.repeat_masked.fa'
	if sp == 'ZF':
        	genome = '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.repeat_masked.switch_masked.fa'
	if sp == 'DBF':
		genome = '/mnt/gluster/home/sonal.singhal1/DBF/masked_genome/DBF.masked_genome.fa'

	# vcf files
	if sp == 'ZF':
		vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.%s.coverage.repeatmasked.filtered.nomendel.shared.noswitch.phased.vqsr2.vcf.gz' % chr
		if chr == 'chrZ':
        		vcf = '/mnt/gluster/home/sonal.singhal1/ZF/after_vqsr/by_chr/unrel_vcf/gatk.ug.unrel_zf.chrZ.coverage.repeatmasked.filtered.nomendel.shared.recodedsex.phased.vqsr2.vcf.gz'
	if sp == 'LTF':
		vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.%s.coverage.repeatmasked.filtered.phased.vqsr2.vcf.gz'  % chr
		if chr == 'chrZ':
        		vcf = '/mnt/gluster/home/sonal.singhal1/LTF/after_vqsr/by_chr/gatk.ug.ltf.chrZ.coverage.repeatmasked.filtered.recodedsex.phased.vqsr2.vcf.gz'
	if sp == 'DBF':
		vcf = '/mnt/gluster/home/sonal.singhal1/DBF/after_vqsr/by_chr/gatk.ug.dbf.%s.filtered.coverage.vqsr.vcf.gz' % chr

	seq = get_chromosome(ref_genome, chr)
	seq = get_variants(vcf, seq)
	mask = get_masked(seq, genome, chr)
	print_seq(sp, chr, seq)
	print 'all done!'

if __name__ == "__main__":
    main()


