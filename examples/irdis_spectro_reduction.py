import vltpf.IRDIS as IRDIS

reduction = IRDIS.SpectroReduction('/Users/avigan/data/VLTPF-test-target/IRD/LSS/')

#
# full reduction
#
reduction.config['combine_science_dim'] = 300
reduction.config['clean'] = True

reduction.full_reduction()

#
# manual reduction
#
# reduction.sort_files()
# reduction.sort_frames()
# reduction.check_files_association()

# reduction.sph_ird_cal_dark()
# reduction.sph_ird_cal_detector_flat()
# reduction.sph_ird_wave_calib()

# reduction.sph_ird_preprocess_science(subtract_background=True, fix_badpix=True, 
#                                      collapse_science=True, collapse_psf=True,
#                                      collapse_center=True)

# reduction.sph_ird_star_center(high_pass=False, display=True, save=True)
# reduction.sph_ird_wavelength_recalibration(fit_scaling=True, display=True, save=True)
# reduction.sph_ird_combine_data(cpix=True, psf_dim=80, science_dim=300, 
#                                correct_mrs_chromatism=True, split_posang=True,
#                                shift_method='fft', manual_center=None, skip_center=False)

# reduction.sph_ird_clean(delete_raw=False, delete_products=False)
