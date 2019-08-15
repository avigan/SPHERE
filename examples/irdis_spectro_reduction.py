import vltpf.IRDIS as IRDIS

reduction = IRDIS.SpectroReduction('/Users/avigan/data/VLTPF-test-target/IRD/LSS/')


reduction.sort_files()
reduction.sort_frames()
