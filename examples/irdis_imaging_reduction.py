import vltpf.IRDIS as IRDIS

reduction = IRDIS.ImagingReduction('/Users/avigan/data/VLTPF-test-target/IRD/DBI/')

#
# full reduction
#

#%% configuration
# reduction.config['combine_psf_dim']          = 80
# reduction.config['combine_science_dim']      = 400
# reduction.config['combine_shift_method']     = 'fft'
# reduction.config['preproc_collapse_science'] = True
# reduction.config['preproc_collapse_type']    = 'mean'

#%% full reduction
# reduction.full_reduction()

#
# manual reduction
#

#%% sorting
reduction.sort_files()
reduction.sort_frames()
reduction.check_files_association()

#%% static calibrations
reduction.sph_ird_cal_dark(silent=True)
reduction.sph_ird_cal_detector_flat(silent=True)

#%% science pre-processing
reduction.sph_ird_preprocess_science(subtract_background=True, fix_badpix=True,
                                     collapse_science=True, collapse_type='mean', coadd_value=2,
                                     collapse_psf=True, collapse_center=True)

#%% high-level science processing
reduction.sph_ird_star_center(high_pass=False, offset=(0, 0), plot=True)
reduction.sph_ird_combine_data(cpix=True, psf_dim=80, science_dim=200, correct_anamorphism=True,
                               shift_method='interp', manual_center=None, skip_center=False,
                               save_scaled=False)

#%% cleaning
reduction.sph_ird_clean(delete_raw=False, delete_products=False)
