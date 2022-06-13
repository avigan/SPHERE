import sphere.IRDIS as IRDIS

####################################################@
# full reduction
#

#%% init reduction
reduction = IRDIS.SpectroReduction('/Users/avigan/data/sphere-test-target/IRD/LSS/', log_level='info')

###############################################################################
# It is possible to provide a default JSON configuration file to set some (or #
# all) of the reduction parameters to a default value different from the ones #
# hard-coded in the sphere package. This is done with the keyword:            #
#   user_config='... path to the file ...'                                    #
# The increasing priority for setting reduction parameters is the following:  #
#   0- default values hard-coded in the sphere package                        #
#   1- values contained in the file pointed by the user_config keyword, if a  #
#      file path is provided and exists                                       #
#   2- values contained in a reduction_config.json file left in the reduction #
#      directory by a previous reduction                                      #
#   3- values manually set by the user (see examples below)                   #
###############################################################################

#%% configuration
reduction.config['combine_science_dim'] = 300
reduction.config['combine_center_selection'] = 'first'
reduction.config['clean'] = False
print(reduction.config)

#%% reduction
reduction.full_reduction()

####################################################@
# manual reduction
#

#%% init reduction
reduction = IRDIS.SpectroReduction('/Users/avigan/data/sphere-test-target/IRD/LSS/', log_level='info')

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
reduction.sph_ird_star_center(high_pass_psf=False, high_pass_waffle=True, plot=True)
reduction.sph_ird_wavelength_recalibration(fit_scaling=True, plot=True)
reduction.sph_ird_combine_data(cpix=True, psf_dim=80, science_dim=300,
                               correct_mrs_chromatism=True, split_posang=True,
                               shift_method='fft', manual_center=None, center_selection='time',
                               coarse_centering=False)

#%% cleaning
reduction.sph_ird_clean(delete_raw=False, delete_products=False, delete_config=False)
