import sys
import pysphere.IRDIS as IRDIS


root_path = '/Users/avigan/data/pySPHERE-test/IRD/'

red = IRDIS.ImagingReduction(root_path)

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
# red.sph_ird_cal_dark(silent=True)
# red.sph_ird_cal_detector_flat(silent=True)
# red.sph_ird_preprocess_science(subtract_background=True, fix_badpix=True,
#                                collapse_science=False, collapse_type='mean', coadd_value=2,
#                                collapse_psf=True, collapse_center=True)
# red.sph_ird_star_center(high_pass=False, display=False, save=True)
# red.sph_ird_combine_data(cpix=True, psf_dim=100, science_dim=400, save_scaled=False)
# red.sph_ird_clean(delete_raw=False, delete_products=False)
