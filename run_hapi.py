dir = '/mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/'

for ix in range(1,23):
	geno = '%schr%s.hapi.geno' % (dir, ix)
	sites = geno.replace('geno', 'sites')
	list = geno.replace('geno', 'list')

	print "echo \"~/bin/hapi-1.03-x86_64/hapi-mr -h -d /mnt/gluster/home/sonal.singhal1/ZF/phasing/family_approach/hapi/ %s %s %s\" | qsub -l h_vmem=3g -cwd -V -j y -N \'hapi%s\'" % (list, sites, geno, ix)
