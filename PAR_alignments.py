import subprocess
import re

fc_contigs = ['AGTO01003702.1', 'JH603380.1', 'JH603441.1']
zf_contigs = ['chr24', 'chrZ', 'chrZ_random']
geo_contigs = ['scaffold143']

for fc in fc_contigs:
	for zf in zf_contigs:
		seq1 = '/mnt/gluster/home/sonal.singhal1/PAR/flycatcher/%s.fa' % fc
		seq2 = '/mnt/gluster/home/sonal.singhal1/PAR/zebrafinch/%s.fa' % zf
		out = '/mnt/gluster/home/sonal.singhal1/PAR/alignments/%s_%s.dotplot.out' % (fc, zf)
		subprocess.call("~/bin/lastz_1.02/lastz %s %s --format=rdotplot > %s" % (seq1, seq2, out), shell=True)
 
for fc in fc_contigs:
        for geo in geo_contigs:
                seq1 = '/mnt/gluster/home/sonal.singhal1/PAR/flycatcher/%s.fa' % fc
                seq2 = '/mnt/gluster/home/sonal.singhal1/PAR/groundfinch/%s.fa' % geo
                out = '/mnt/gluster/home/sonal.singhal1/PAR/alignments/%s_%s.dotplot.out' % (fc, geo)
                subprocess.call("~/bin/lastz_1.02/lastz %s %s --format=rdotplot > %s" % (seq1, seq2, out), shell=True)

for zf in zf_contigs:
        for geo in geo_contigs:
                seq1 = '/mnt/gluster/home/sonal.singhal1/PAR/zebrafinch/%s.fa' % zf
                seq2 = '/mnt/gluster/home/sonal.singhal1/PAR/groundfinch/%s.fa' % geo
                out = '/mnt/gluster/home/sonal.singhal1/PAR/alignments/%s_%s.dotplot.out' % (zf, geo)
                subprocess.call("~/bin/lastz_1.02/lastz %s %s --format=rdotplot > %s" % (seq1, seq2, out), shell=True)
