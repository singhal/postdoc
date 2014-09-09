import glob
import subprocess
import re
import os.path

dir = '/mnt/gluster/home/sonal.singhal1/LTF/analysis/LDhat/'
num_jobs = 20
files = glob.glob(dir + 'sites_locs/*locs')

num_cmds = int(len(files)/num_jobs) + 1

for ix in range(num_jobs):
	out = dir + 'ldhat_runs%s.sh' % ix
	o = open(out, 'w')
	for loc_file in files[ix * num_cmds: ix*num_cmds + num_cmds]:
		site_file = loc_file.replace('.locs', '.sites')
		prefix = re.search('(chr.*_chunk\d+)', site_file).group(1)

		expected_out = dir + 'results/' + prefix + 'rates.txt'

		if not os.path.isfile(expected_out):
			call = '/mnt/lustre/home/sonal.singhal1/bin/LDhat_v2.2/interval -seq %s -loc %s -lk /mnt/gluster/home/sonal.singhal1/LTF/analysis/LDhat/LTFnew_lk.txt -its 5000000 -samp 5000 -bpen 5.0 -prefix %s' % (site_file, loc_file, dir + 'results/' + prefix)

			o.write(call + '\n')
	o.close()
	subprocess.call('chmod a+x %s' % out, shell=True)
	subprocess.call('echo \"%s\" | qsub -l h_vmem=500M -cwd -V -j y -N \'%s\'' % (out, 'ldhat' + str(ix)), shell=True)
