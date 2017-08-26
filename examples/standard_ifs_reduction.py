import sys
import pysphere.IFS as IFS


root_path = '/Users/avigan/data/pySPHERE-test/IFS_sub/'

red = IFS.IFSReduction(root_path)

# one-line reduction
# red.full_reduction()

# standard manual reduction
# red.init_reduction()
# red.create_static_calibrations()
# red.preprocess_science()
# red.process_science()
# red.clean()

# completely manual reduction (full control)
# red.sort_files()
# red.sort_frames()
# red.check_files_association()
# red.sph_ifs_cal_dark(silent=True)
# red.sph_ifs_cal_detector_flat()
# red.sph_ifs_cal_specpos(silent=True)
# red.sph_ifs_cal_wave(silent=True)
# red.sph_ifs_cal_ifu_flat(silent=True)
# red.sph_ifs_preprocess_science(subtract_background=True, fix_badpix=True, correct_xtalk=True,
#                                collapse_science=True, collapse_type='mean', coadd_value=2,
#                                collapse_psf=True, collapse_center=True)
# red.sph_ifs_preprocess_wave()
# red.sph_ifs_science_cubes(postprocess=True, silent=True)
# red.sph_ifs_wavelength_recalibration(high_pass=False, display=False, save=True)
# red.sph_ifs_star_center(high_pass=False, display=False, save=True)
# red.sph_ifs_combine_data(cpix=True, psf_dim=80, science_dim=290, save_scaled=True)
# red.sph_ifs_clean(delete_raw=False, delete_products=False)
