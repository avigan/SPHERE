import vltpf.IRDIS as IRDIS

reduction = IRDIS.SpectroReduction('/Users/avigan/data/VLTPF-test-target/IRD/LSS/')

# reduction.sort_files()
# reduction.sort_frames()
reduction.check_files_association()

# reduction.sph_ird_cal_dark()
# reduction.sph_ird_cal_detector_flat()
# reduction.sph_ird_wave_calib()
