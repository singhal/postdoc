dir = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/'

for ix in range(1,23):
	geno = '%sinput_files/chr%s.hapi.noswitch.geno' % (dir, ix)
	sites = geno.replace('geno', 'sites')
	list = geno.replace('geno', 'list')

	o = open('hapi%s.sh' % ix, 'w')

	o.write("~/bin/hapi-1.03-x86_64/hapi-mr -h -d /mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/ %s %s %s\n" % (list, sites, geno))
	name = ix
	if name < 10:
		name = '0' + str(name)
	o.write("mv /mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/hap-fam1-12-10.%s.csv /mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/output_files/chr%s.noswitch.csv\n" % (name, ix))
	o.close()

	print "echo \'sh hapi%s.sh\' | qsub -l h_vmem=5g -cwd -V -j y -N hapi%s" % (ix,ix)

chrs = ['chr1A', 'chr1B', 'chr4A', 'chr23', 'chr24', 'chr25', 'chr26', 'chr27', 'chr28', 'chrLG2', 'chrLG5', 'chrLGE22']

o = open('hapiX.sh', 'w')
for chr in reversed(chrs):
	o.write("~/bin/hapi-1.03-x86_64/hapi-mr -h -d /mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/output_files/ /mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/input_files/%s.hapi.noswitch.list /mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/input_files/%s.hapi.noswitch.sites /mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/input_files/%s.hapi.noswitch.geno\n" % (chr, chr, chr))
	o.write("mv /mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/output_files/hap-fam1-12-10.01.csv /mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/output_files/%s.noswitch.csv\n" % chr)
o.close()
print "echo \'sh hapiX.sh\' | qsub -l h_vmem=5g -cwd -V -j y -N hapiX"
