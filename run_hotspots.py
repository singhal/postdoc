chrs = [	'chr1', 'chr1A', 'chr2', 'chr3', 'chr4', 'chr4A', 'chr5', 'chr6', \
	'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', \
	'chr14', 'chr15', 'chrZ']

for block in [500,1000,2000]:
	for flank in [20000,40000]:
		for chr in chrs:
			print 'python ~/scripts/find_hotspots.py --sp LTF --chr %s --block %s --flank %s' % (chr, block, flank)
