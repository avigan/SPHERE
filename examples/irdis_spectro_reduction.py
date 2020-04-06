import pysphere.IRDIS as IRDIS

####################################################@
# full reduction
#

#%% init reduction
reduction = IRDIS.SpectroReduction('/Users/avigan/data/pysphere-test-target/IRD/LSS/', log_level='info')

#%% configuration
reduction.config['combine_science_dim'] = 300
reduction.config['clean'] = True
reduction.show_config()

#%% reduction
reduction.full_reduction()

####################################################@
# manual reduction
#

#%% init reduction
reduction = IRDIS.SpectroReduction('/Users/avigan/data/pysphere-test-target/IRD/LSS/', log_level='info')

#%% sorting
reduction.sort_files()
reduction.sort_frames()
reduction.check_files_association()

#%% static calibrations
reduction.sph_ird_cal_dark(silent=True)
reduction.sph_ird_cal_detector_flat(silent=True)
reduction.sph_ird_cal_wave(silent=True)

#%% science pre-processing
reduction.sph_ird_preprocess_science(subtract_background=True, fix_badpix=True,
                                     collapse_science=True, collapse_psf=True,
                                     collapse_center=True)

#%% high-level science processing
reduction.sph_ird_star_center(high_pass=False, plot=True)
reduction.sph_ird_wavelength_recalibration(fit_scaling=True, plot=True)
reduction.sph_ird_combine_data(cpix=True, psf_dim=80, science_dim=300,
                               correct_mrs_chromatism=True, split_posang=True,
                               shift_method='fft', manual_center=None, coarse_centering=False)

#%% cleaning
reduction.sph_ird_clean(delete_raw=False, delete_products=False)
