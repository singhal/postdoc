import glob
import re
import subprocess
import random
import os

out_dir = '/mnt/gluster/home/sonal.singhal1/ZF/analysis/dnds/'
aligned_files = glob.glob('%salignments/*phy' % out_dir)

def make_trees(out_dir, aligned_files):
        out_dir = out_dir + 'trees/'
        for file in aligned_files:
                name = re.search('([A-Z|0-9]+).phy', file).group(1)
                rand_num = random.randint(1,10000)
                rand_num2 = random.randint(1,10000)
                subprocess.call('raxmlHPC -f a -m GTRCAT -n %s -s %s -w %s -N 100 -x %s -p %s' % (name, file, out_dir, rand_num, rand_num2), shell=True)

def run_paml(out_dir, aligned_files):
        constant_seq = '%sseq.phy' % out_dir
        constant_tree = '%sseq.tre' % out_dir
        constant_out = '%sout' % out_dir

        for file in aligned_files:
                name = re.search('([A-Z|0-9]+).phy', file).group(1)
                tree = '%strees/RAxML_bestTree.%s' % (out_dir, name)
                out = '%sresults/%s.paml.out' % (out_dir, name)

                if os.path.isfile(tree):
                        subprocess.call('cp %s %s' % (file, constant_seq), shell=True)
                        subprocess.call('cp %s %s' % (tree, constant_tree), shell=True)
                        subprocess.call('~/bin/paml4.8/bin/codeml ~/bin/paml4.8/bin/codeml.ctl', shell=True)
                        subprocess.call('cp %s %s' % (constant_out, out), shell=True)

# make_trees(out_dir, aligned_files)
run_paml(out_dir, aligned_files)

