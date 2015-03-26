import subprocess
import pandas as pd
import re

genome = '/mnt/gluster/home/sonal.singhal1/ZF/masked_genome/ZF.masked_genome.repeat_masked.switch_masked.fa'
d = pd.read_csv('/mnt/gluster/home/sonal.singhal1/ZF/analysis/hotspots/spot2kb_flank40kb.seqldhot_validate_hotspots.csv') 
d = d[d.zlk >= 10]
d = d.groupby('chr')

chr_lengths = { 'chr10': 20806668, 'chr11': 21403021, 'chr12': 21576510, 'chr13': 16962381,
                'chr14': 16419078, 'chr15': 14428146, 'chr16': 9909, 'chr17': 11648728,
                'chr18': 11201131, 'chr19': 11587733, 'chr1A': 73657157, 'chr1B': 1083483,
                'chr1': 118548696, 'chr20': 15652063, 'chr21': 5979137, 'chr22': 3370227,
                'chr23': 6196912, 'chr24': 8021379, 'chr25': 1275379, 'chr26': 4907541,
                'chr27': 4618897, 'chr28': 4963201, 'chr2': 156412533, 'chr3': 112617285,
                'chr4A': 20704505, 'chr4': 69780378, 'chr5': 62374962, 'chr6': 36305782,
                'chr7': 39844632, 'chr8': 27993427, 'chr9': 27241186, 'chrLG2': 109741,
                'chrLG5': 16416, 'chrLGE22': 883365, 'chrZ': 72861351}

for chr, chrgroup in d:
    for hot_start in chrgroup.zstart:
        start = hot_start - 25000
        if start < 1:
            start = 1
        end = hot_start + 25000
        if end > chr_lengths[chr]:
            end = chr_lengths[chr]
            
        seq = subprocess.Popen('samtools faidx %s %s:%s-%s' % (genome, chr, start, end), shell=True, stdout=subprocess.PIPE)
        n_count = 0
        for l in seq.stdout:
            if not re.match('>', l):
                n_count += len(filter(lambda x: int(x) < 4, list(l.rstrip())))
        print n_count / float(50000)
