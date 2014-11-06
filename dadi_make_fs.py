import dadi

dd = dadi.Misc.make_data_dict('/mnt/gluster/home/sonal.singhal1/ZF/analysis/dadi/all_species.dadi_pruned.txt')

fs1 = dadi.Spectrum.from_data_dict(dd, ['ZF', 'LTFh'], [38, 20])
fs2 = dadi.Spectrum.from_data_dict(dd, ['ZF', 'LTFa'], [38, 20])
fs3 = dadi.Spectrum.from_data_dict(dd, ['LTFh', 'LTFa'], [20, 20])
fs4 = dadi.Spectrum.from_data_dict(dd, ['ZF', 'LTFa', 'LTFh'], [38, 20, 20])

fs1.to_file('/mnt/gluster/home/sonal.singhal1/ZF/analysis/dadi/ZF_LTFh.fs')
fs2.to_file('/mnt/gluster/home/sonal.singhal1/ZF/analysis/dadi/ZF_LTFa.fs')
fs3.to_file('/mnt/gluster/home/sonal.singhal1/ZF/analysis/dadi/LTFa_LTFh.fs')
fs4.to_file('/mnt/gluster/home/sonal.singhal1/ZF/analysis/dadi/all.fs')
