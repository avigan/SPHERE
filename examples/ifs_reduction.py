import sphere.IFS as IFS

#####################################################################
# Starting in the present version of the pipeline, the default      #
# -1.75° true North offset is automatically added to the derotation #
# angles. The offset value can be modified in the configuration of  #
# the reduction:                                                    #
#                                                                   #
#   >>> reduction.config[\'cal_true_north\'] = xxx                  #
#                                                                   #
# To avoid any issues, make sure to:                                #
#   * either reprocess data previously processed with version <1.4  #
#   * or take into account the offset in your astrometric analysis  #
#####################################################################

####################################################@
# full reduction
#

#%% init reduction
reduction = IFS.Reduction('/Users/avigan/data/sphere-test-target/IFS/', log_level='info')

#%% configuration
reduction.config['preproc_collapse_science'] = True
reduction.config['preproc_collapse_type']    = 'coadd'
reduction.config['preproc_coadd_value']      = 2
reduction.config['combine_center_selection'] = 'first'
reduction.config['center_high_pass_waffle']  = True
reduction.config['clean']                    = False
reduction.show_config()

#%% reduction
reduction.full_reduction()

####################################################
# manual reduction
#

#%% init reduction
reduction = IFS.Reduction('/Users/avigan/data/sphere-test-target/IFS/', log_level='info')

#%% sorting
reduction.sort_files()
reduction.sort_frames()
reduction.check_files_association()

#%% static calibrations
reduction.sph_ifs_cal_dark(silent=True)
reduction.sph_ifs_cal_detector_flat(silent=True)
reduction.sph_ifs_cal_specpos(silent=True)
reduction.sph_ifs_cal_wave(silent=True)
reduction.sph_ifs_cal_ifu_flat(silent=True)

#%% science pre-processing
reduction.sph_ifs_preprocess_science(subtract_background=True, fix_badpix=True, correct_xtalk=True,
                                     collapse_science=True, collapse_type='mean', coadd_value=2,
                                     collapse_psf=True, collapse_center=True)
reduction.sph_ifs_preprocess_wave()
reduction.sph_ifs_science_cubes(silent=True)

#%% high-level science processing
reduction.sph_ifs_wavelength_recalibration(high_pass=True, offset=(-3, 0), plot=True)
reduction.sph_ifs_star_center(high_pass_psf=False, high_pass_waffle=True, offset=(-3, 0), plot=True)
reduction.sph_ifs_combine_data(cpix=True, psf_dim=80, science_dim=200, correct_anamorphism=True,
                               shift_method='interp', manual_center=None, center_selection='time',
                               coarse_centering=False, save_scaled=False)

#%% cleaning
reduction.sph_ifs_clean(delete_raw=False, delete_products=False)
