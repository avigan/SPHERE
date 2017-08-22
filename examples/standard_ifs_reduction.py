import sys
sys.path.append('/Users/avigan/Work/GitHub/pySPHERE/')
import pysphere.IFS as IFS


root_path = '/Users/avigan/data/pySPHERE-test/IFS/'

# files_info = IFS.sort_files(root_path)
# frames_info = IFS.sort_frames(root_path, files_info)

# IFS.check_files_association(root_path, files_info)

# files_info, frames_info, frames_info_preproc = IFS.read_info(root_path)
# IFS.sph_ifs_cal_dark(root_path, files_info)
# IFS.sph_ifs_cal_detector_flat(root_path, files_info)
# IFS.sph_ifs_cal_specpos(root_path, files_info)
# IFS.sph_ifs_cal_wave(root_path, files_info)
# IFS.sph_ifs_cal_ifu_flat(root_path, files_info)

# files_info, frames_info, frames_info_preproc = IFS.read_info(root_path)
# IFS.sph_ifs_preprocess_science(root_path, files_info, frames_info,
#                            subtract_background=True, fix_badpix=True, correct_xtalk=True,
#                            collapse_science=False, collapse_type='mean', coadd_value=2,
#                            collapse_psf=True, collapse_center=True)
# IFS.sph_ifs_preprocess_wave(root_path, files_info)

# files_info, frames_info, frames_info_preproc = IFS.read_info(root_path)
# IFS.sph_ifs_science_cubes(root_path, files_info, frames_info_preproc)

# wave = IFS.sph_ifs_wavelength_recalibration(root_path)

# IFS.sph_ifs_star_center(root_path)

# IFS.sph_ifs_combine_data(root_path, save_scaled=True)

red = IFS.IFSReduction(root_path)

red.full_reduction()
