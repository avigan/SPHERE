import vltpf.IRDIS as IRDIS

reduction = IRDIS.ImagingReduction('/Users/avigan/data/VLTPF-test-target/IRD/DBI/')

# reduction.config['combine_psf_dim']          = 200
# reduction.config['combine_science_dim']      = 200
# reduction.config['combine_shift_method']     = 'interp'
# reduction.config['preproc_collapse_science'] = True
# reduction.config['preproc_collapse_type']    = 'mean'
# reduction.full_reduction()

reduction.sort_files()
reduction.sort_frames()
# reduction.check_files_association()
# reduction.sph_ird_cal_dark(silent=True)
# reduction.sph_ird_cal_detector_flat(silent=True)
# reduction.sph_ird_preprocess_science(subtract_background=True, fix_badpix=True,
#                                      collapse_science=True, collapse_type='mean', coadd_value=2,
#                                      collapse_psf=True, collapse_center=True)
# reduction.sph_ird_star_center(high_pass=False, offset=(0,0), display=False, save=True)
# reduction.sph_ird_combine_data(cpix=True, psf_dim=200, science_dim=200, correct_anamorphism=True,
#                                shift_method='interp', manual_center=None, skip_center=False,
#                                save_scaled=False)
# reduction.sph_ird_clean(delete_raw=False, delete_products=False)
