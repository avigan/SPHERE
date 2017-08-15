import os
import glob
import pandas as pd
import subprocess
import numpy as np
import astropy.coordinates as coord
import astropy.units as units
import scipy.ndimage as ndimage

import imutils

from astropy.io import fits
from astropy.time import Time

# keywords to be saved
keywords = [
    # standard
    'INSTRUME', 
    'OBJECT', 'DATE-OBS', 'DATE', 'HIERARCH ESO DET FRAM UTC',

    # DPR
    'HIERARCH ESO DPR CATG', 'HIERARCH ESO DPR TYPE', 'HIERARCH ESO DPR TECH',
        
    # coordinates
    'HIERARCH ESO TEL GEOLAT', 'HIERARCH ESO TEL GEOLON', 'HIERARCH ESO TEL GEOELEV',
    'HIERARCH ESO INS4 DROT2 RA', 'HIERARCH ESO INS4 DROT2 DEC',
    'HIERARCH ESO TEL ALT', 'HIERARCH ESO TEL AZ',

    # SAXO
    'HIERARCH ESO AOS TTLOOP STATE', 'HIERARCH ESO AOS HOLOOP STATE',
    'HIERARCH ESO AOS IRLOOP STATE', 'HIERARCH ESO AOS PUPLOOP STATE',
    'HIERARCH ESO AOS VISWFS MODE',

    # CPI
    'HIERARCH ESO SEQ ARM',
    'HIERARCH ESO INS COMB ICOR', 'HIERARCH ESO INS COMB IFLT', 'HIERARCH ESO INS COMB POLA',
    'HIERARCH ESO INS4 FILT2 NAME',
    'HIERARCH ESO INS4 DROT2 BEGIN', 'HIERARCH ESO INS4 DROT2 END', 'HIERARCH ESO INS4 DROT2 MODE', 
    
    # IFS
    'HIERARCH ESO INS2 MODE', 'HIERARCH ESO INS2 COMB IFS',

    # IRDIS
    'HIERARCH ESO INS1 MODE',
    'HIERARCH ESO INS1 FILT NAME', 'HIERARCH ESO INS1 OPTI2 NAME',
    'HIERARCH ESO INS1 PAC X', 'HIERARCH ESO INS1 PAC Y',
    
    # detector
    'HIERARCH ESO DET SEQ1 DIT', 'HIERARCH ESO DET NDIT',

    # observing conditions
    'HIERARCH ESO TEL AIRM START', 'HIERARCH ESO TEL AIRM END',
    'HIERARCH ESO TEL AMBI FWHM START', 'HIERARCH ESO TEL AMBI FWHM END', 'HIERARCH ESO TEL IA FWHM',
    'HIERARCH ESO TEL AMBI TAU0', 'HIERARCH ESO TEL AMBI TEMP',
    'HIERARCH ESO TEL AMBI WINDSP', 'HIERARCH ESO TEL AMBI WINDDIR'
]

# short keywords
keywords_short = keywords.copy()
for idx in range(len(keywords_short)):
    key = keywords_short[idx]
    if key.find('HIERARCH ESO ') != -1:
        keywords_short[idx] = key[13:]
        

def read_info(root_path):
    '''
    Read the files, calibs and frames information

    Parameters
    ----------
    root_path : str
        Path to the dataset

    Returns
    -------
    files_info : dataframe
        The data frame with all the information on raw science files

    calibs_info : dataframe
        The data frame with all the information on calibration files

    frames_info : dataframe
        The data frame with all the information on science frames

    frames_info_preproc : dataframe
        The data frame with all the information on science frames after pre-processing
    '''

    # files info
    fname = os.path.join(root_path, 'files.csv')
    if os.path.exists(fname):
        files_info = pd.read_csv(fname, index_col=0)

        # convert times
        files_info['DATE-OBS'] = pd.to_datetime(files_info['DATE-OBS'], utc=True)
        files_info['DATE'] = pd.to_datetime(files_info['DATE'], utc=True)
        files_info['DET FRAM UTC'] = pd.to_datetime(files_info['DET FRAM UTC'], utc=True)
    else:
        files_info = None

    fname = os.path.join(root_path, 'calibs.csv')
    if os.path.exists(fname):
        calibs_info = pd.read_csv(fname, index_col=0)
    else:
        calibs_info = None

    fname = os.path.join(root_path, 'frames.csv')
    if os.path.exists(fname):
        frames_info = pd.read_csv(fname, index_col=(0, 1))

        # convert times
        frames_info['DATE-OBS'] = pd.to_datetime(frames_info['DATE-OBS'], utc=True)
        frames_info['DATE'] = pd.to_datetime(frames_info['DATE'], utc=True)
        frames_info['DET FRAM UTC'] = pd.to_datetime(frames_info['DET FRAM UTC'], utc=True)
        frames_info['TIME START'] = pd.to_datetime(frames_info['TIME START'], utc=True)
        frames_info['TIME'] = pd.to_datetime(frames_info['TIME'], utc=True)
        frames_info['TIME END'] = pd.to_datetime(frames_info['TIME END'], utc=True)
    else:
        frames_info = None

    fname = os.path.join(root_path, 'frames_preproc.csv')
    if os.path.exists(fname):
        frames_info_preproc = pd.read_csv(fname, index_col=(0, 1))

        # convert times
        frames_info_preproc['DATE-OBS'] = pd.to_datetime(frames_info_preproc['DATE-OBS'], utc=True)
        frames_info_preproc['DATE'] = pd.to_datetime(frames_info_preproc['DATE'], utc=True)
        frames_info_preproc['DET FRAM UTC'] = pd.to_datetime(frames_info_preproc['DET FRAM UTC'], utc=True)
        frames_info_preproc['TIME START'] = pd.to_datetime(frames_info_preproc['TIME START'], utc=True)
        frames_info_preproc['TIME'] = pd.to_datetime(frames_info_preproc['TIME'], utc=True)
        frames_info_preproc['TIME END'] = pd.to_datetime(frames_info_preproc['TIME END'], utc=True)
    else:
        frames_info_preproc = None
           
    return files_info, calibs_info, frames_info, frames_info_preproc

    
def sort_files(root_path):
    '''
    Sort all raw files and save result in a data frame

    Parameters
    ----------
    root_path : str
        Path to the dataset

    Returns
    -------
    files_info : dataframe
        Data frame with the information on raw files
    '''

    print('Sorting raw files')
    
    # list files
    files = glob.glob(os.path.join(root_path, 'raw', '*.fits'))
    files = [os.path.splitext(os.path.basename(f))[0] for f in files]
    
    if len(files) == 0:
        raise ValueError('No raw FITS files in root_path directory')

    print(' * found {0} FITS files in {1}'.format(len(files), root_path))

    # files table
    files_info = pd.DataFrame(index=pd.Index(files, name='FILE'), columns=keywords_short, dtype='float')

    raw_path = os.path.join(root_path, 'raw/')
    for f in files:
        hdu = fits.open(os.path.join(raw_path, f+'.fits'))
        hdr = hdu[0].header

        for k, sk in zip(keywords, keywords_short):
            files_info.loc[f, sk] = hdr.get(k)

        hdu.close()

    # processed column
    files_info.insert(len(files_info.columns), 'PROCESSED', False)

    # convert times
    files_info['DATE-OBS'] = pd.to_datetime(files_info['DATE-OBS'], utc=True)
    files_info['DATE'] = pd.to_datetime(files_info['DATE'], utc=True)
    files_info['DET FRAM UTC'] = pd.to_datetime(files_info['DET FRAM UTC'], utc=True)
    
    # save files_info
    files_info.to_csv(os.path.join(root_path, 'files.csv'))    
    
    return files_info


def parallatic_angle(ha, dec, geolat):
    '''
    Parallactic angle of a source in degrees

    Parameters
    ----------
    ha : array_like
        Hour angle, in hours

    dec : float
        Declination, in degrees

    geolat : float
        Observatory declination, in degrees

    Returns
    -------
    pa : array_like
        Parallactic angle values
    '''
    pa = -np.arctan2(-np.sin(ha),
                     np.cos(dec) * np.tan(geolat) - np.sin(dec) * np.cos(ha))

    if (dec >= geolat):
        pa[ha < 0] += 360*units.degree
    
    return np.degrees(pa)


def compute_times(frames_info):
    '''
    Compute the various timestamps associated to frames

    Parameters
    ----------
    frames_info : dataframe
        The data frame with all the information on science frames
    '''

    # get necessary values
    time_start = frames_info['DATE-OBS'].values
    time_end   = frames_info['DET FRAM UTC'].values
    time_delta = (time_end - time_start) / frames_info['DET NDIT'].values.astype(np.int)
    DIT        = np.array(frames_info['DET SEQ1 DIT'].values.astype(np.float)*1000, dtype='timedelta64[ms]')

    # calculate time stamps
    idx = frames_info.index.get_level_values(1).values
    ts_start = time_start + time_delta * idx
    ts       = time_start + time_delta * idx + DIT/2
    ts_end   = time_start + time_delta * idx + DIT

    # update frames_info
    frames_info['TIME START'] = ts_start
    frames_info['TIME']       = ts
    frames_info['TIME END']   = ts_end


def compute_angles(frames_info):
    '''
    Compute the various angles associated to frames: RA, DEC, parang,
    pupil offset, final derotation angle

    Parameters
    ----------
    frames_info : dataframe
        The data frame with all the information on science frames
    '''

    # RA/DEC
    ra_drot = frames_info['INS4 DROT2 RA'].values.astype(np.float)
    ra_drot_h = np.floor(ra_drot/1e4)
    ra_drot_m = np.floor((ra_drot - ra_drot_h*1e4)/1e2)
    ra_drot_s = ra_drot - ra_drot_h*1e4 - ra_drot_m*1e2
    ra = coord.Angle((ra_drot_h, ra_drot_m, ra_drot_s), units.hour)
    frames_info['RA'] = ra

    dec_drot = frames_info['INS4 DROT2 DEC'].values.astype(np.float)
    sign = np.sign(dec_drot)
    udec_drot  = np.abs(dec_drot)
    dec_drot_d = np.floor(udec_drot/1e4)
    dec_drot_m = np.floor((udec_drot - dec_drot_d*1e4)/1e2)
    dec_drot_s = udec_drot - dec_drot_d*1e4 - dec_drot_m*1e2
    dec_drot_d *= sign
    dec = coord.Angle((dec_drot_d, dec_drot_m, dec_drot_s), units.degree)
    frames_info['DEC'] = dec
    
    # calculate parallactic angles
    geolon = coord.Angle(frames_info['TEL GEOLON'].values[0], units.degree)
    geolat = coord.Angle(frames_info['TEL GEOLAT'].values[0], units.degree)
    geoelev = frames_info['TEL GEOELEV'].values[0]

    utc = Time(frames_info['TIME START'].values.astype(str), scale='utc', location=(geolon, geolat, geoelev))
    lst = utc.sidereal_time('apparent')
    ha  = lst - ra
    pa  = parallatic_angle(ha, dec[0], geolat)    
    frames_info['PARANG START'] = pa

    utc = Time(frames_info['TIME'].values.astype(str), scale='utc', location=(geolon, geolat, geoelev))
    lst = utc.sidereal_time('apparent')
    ha  = lst - ra
    pa  = parallatic_angle(ha, dec[0], geolat)    
    frames_info['PARANG'] = pa

    utc = Time(frames_info['TIME END'].values.astype(str), scale='utc', location=(geolon, geolat, geoelev))
    lst = utc.sidereal_time('apparent')
    ha  = lst - ra
    pa  = parallatic_angle(ha, dec[0], geolat)    
    frames_info['PARANG END'] = pa

    # pupil offset
    # PA_on-sky = PA_detector + PARANGLE + True_North + PUPOFFSET + IFSOFFSET
    #   PUPOFFSET = 135.99±0.11
    #   IFSOFFSET = 100.48±0.0.10
    drot_mode = frames_info['INS4 DROT2 MODE'].unique()
    if len(drot_mode) != 1:
        raise ValueError('Derotator mode has several values in the sequence')
    if drot_mode == 'ELEV':
        pupoff = 135.99 - 100.48
    elif drot_mode == 'SKY':
        pupoff = -100.48
    else:
        raise ValueError('Unknown derotator mode {0}'.format(drot_mode))

    frames_info['PUPIL OFFSET'] = pupoff

    # final derotation value
    frames_info['DEROT ANGLES'] = frames_info['PARANG'] + pupoff
    

def sort_frames(root_path, files_info):
    '''
    Extract the frames information from the science files

    Parameters
    ----------
    root_path : str
        Path to the dataset
    
    files_info : dataframe
        The data frame with all the information on raw files

    Returns
    -------
    calibs : dataframe
        A data frame with the information on all frames
    '''

    print('Extracting frames information')

    # convert times
    files_info['DATE-OBS'] = pd.to_datetime(files_info['DATE-OBS'], utc=True)
    files_info['DATE'] = pd.to_datetime(files_info['DATE'], utc=True)
    files_info['DET FRAM UTC'] = pd.to_datetime(files_info['DET FRAM UTC'], utc=True)
        
    # science files
    sci_files = files_info[(files_info['DPR CATG'] == 'SCIENCE') & (files_info['DPR TYPE'] != 'SKY')]    

    # build indices
    files = []
    img   = []
    for file, finfo in sci_files.iterrows():
        NDIT = int(finfo['DET NDIT'])

        files.extend(np.repeat(file, NDIT))
        img.extend(list(np.arange(NDIT)))

    # create new dataframe
    frames_info = pd.DataFrame(columns=sci_files.columns, index=pd.MultiIndex.from_arrays([files, img], names=['FILE', 'IMG']))
    
    # expand files_info into frames_info
    frames_info = frames_info.align(files_info, level=0)[1]    

    # compute timestamps
    compute_times(frames_info)
    
    # compute angles (ra, dec, parang)
    compute_angles(frames_info)
        
    # save
    frames_info.to_csv(os.path.join(root_path, 'frames.csv'))
    
    return frames_info


def files_association(root_path, files_info):
    '''
    Performs the calibration files association and a sanity check

    Parameters
    ----------
    root_path : str
        Path to the dataset
    
    files_info : dataframe
        The data frame with all the information on raw files

    Returns
    -------
    calibs : dataframe
        A data frame with the information on all calibrations
    '''

    print('Performing file association for calibrations')
    
    # IFS obs mode
    modes = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS2 COMB IFS']
    if len(modes.unique()) != 1:
        raise ValueError('Sequence is mixing YJ and YJH observations.')

    mode = modes.unique()[0]    
    if mode == 'OBS_YJ':
        mode_short = 'YJ'
    elif mode == 'OBS_H':
        mode_short = 'YJH'
    else:
        raise ValueError('Unknown IFS mode {0}'.format(mode))

    # specific data frame for calibrations
    # keep static calibrations and sky backgrounds
    calibs = files_info[(files_info['DPR CATG'] == 'CALIB') |
                        ((files_info['DPR CATG'] == 'SCIENCE') & (files_info['DPR TYPE'] == 'SKY'))]
    
    ###############################################
    # static calibrations not dependent on science
    ###############################################
    error_flag = 0
    warning_flag = 0
    
    # white flat
    cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == 'CAL_BB_2_{0}'.format(mode_short))]
    if len(cfiles) != 2:
        error_flag += 1
        print(' * Error: there should be 2 flat files for white lamp!')
        
    # 1020 nm flat
    cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == 'CAL_NB1_1_{0}'.format(mode_short))]
    if len(cfiles) != 2:
        error_flag += 1
        print(' * Error: there should be 2 flat files for 1020 nm filter!')

    # 1230 nm flat
    cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == 'CAL_NB2_1_{0}'.format(mode_short))]
    if len(cfiles) != 2:
        error_flag += 1
        print(' * Error: there should be 2 flat files for 1230 nm filter!')

    # 1300 nm flat
    cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == 'CAL_NB3_1_{0}'.format(mode_short))]
    if len(cfiles) != 2:
        error_flag += 1
        print(' * Error: there should be 2 flat files for 1300 nm filter!')

    # 1550 nm flat (YJH mode only)
    if mode_short == 'YJH':
        cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == 'CAL_NB4_2_{0}'.format(mode_short))]
        if len(cfiles) != 2:
            error_flag += 1
            print(' * Error: there should be 2 flat files for 1550 nm filter!')

    # spectra position
    cfiles = calibs[(calibs['DPR TYPE'] == 'SPECPOS,LAMP') & (calibs['INS2 COMB IFS'] == mode)]
    if len(cfiles) != 1:
        error_flag += 1
        print(' * Error: there should be 1 spectra position file!')

    # wavelength
    cfiles = calibs[(calibs['DPR TYPE'] == 'WAVE,LAMP') & (calibs['INS2 COMB IFS'] == mode)]
    if len(cfiles) != 1:
        error_flag += 1
        print(' * Error: there should be 1 wavelength calibration file!')
    
    # IFU flat
    cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == mode)]
    if len(cfiles) != 1:
        error_flag += 1
        print(' * Error: there should be 1 IFU flat file!')
    
    # calibs dark file
    cfiles = calibs[((calibs['DPR TYPE'] == 'DARK') | (calibs['DPR TYPE'] == 'DARK,BACKGROUND')) &
                    (calibs['DET SEQ1 DIT'].round(2) == 1.65)]
    if len(cfiles) != 1:
        warning_flag += 1
        print(' * Warning: there is no dark/background for the basic calibrations (DIT=1.6 sec). ' +
              'It is *highly recommended* to include one to obtain the best data reduction. ' +
              'A single dark/background file is sufficient, and it can easily be downloaded ' +
              'from the ESO archive')
    
    ##################################################
    # static calibrations that depend on science (DIT)
    ##################################################

    obj = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'DPR TYPE'].apply(lambda s: s[0:6])
    DITs = files_info.loc[(files_info['DPR CATG'] == 'SCIENCE') & (obj == 'OBJECT'), 'DET SEQ1 DIT'].unique().round(2)

    # handle darks in a slightly different way because there might be several different DITs
    for DIT in DITs:
        # instrumental backgrounds
        cfiles = calibs[((calibs['DPR TYPE'] == 'DARK') | (calibs['DPR TYPE'] == 'DARK,BACKGROUND')) &
                        (calibs['DET SEQ1 DIT'].round(2) == DIT)]
        if len(cfiles) == 0:
            warning_flag += 1
            print(' * Warning: there is no dark/background for science files with DIT={0} sec. '.format(DIT) +
                  'it is *highly recommended* to include one to obtain the best data reduction. ' +
                  'A single dark/background file is sufficient, and it can easily be downloaded ' +
                  'from the ESO archive')
            
        # sky backgrounds
        cfiles = files_info[(files_info['DPR TYPE'] == 'SKY') & (files_info['DET SEQ1 DIT'].round(2) == DIT)]
        if len(cfiles) == 0:
            warning_flag += 1
            print(' * Warning: there is no sky background for science files with DIT={0} sec. '.format(DIT) +
                  'Using a sky background instead of an internal instrumental background can ' +
                  'usually provide a cleaner data reduction')

    # keep only some columns
    calibs = calibs[['DPR CATG', 'DPR TYPE', 'DPR TECH', 'INS2 COMB IFS', 'DET SEQ1 DIT', 'DET NDIT', 'PROCESSED']]
            
    # save
    calibs.to_csv(os.path.join(root_path, 'calibs.csv'))

    # error reporting
    print('There are {0} warnings and {1} error in the classification of files'.format(warning_flag, error_flag))
    if error_flag:
        raise ValueError('There is {0} errors that should be solved before proceeding'.format(error_flag))
    
    return calibs


def sph_ifs_cal_dark(root_path, calibs_info, silent=True):
    '''
    Create the dark and background calibrations

    Parameters
    ----------
    root_path : str
        Path to the dataset

    calibs_info : dataframe
        The data frame with all the information on calibration files

    silent : bool
        Suppress esorex output. Optional, default is True
    '''

    print('Creating darks and backgrounds')
    
    # check directories
    raw_path = os.path.join(root_path, 'raw/')

    calib_path = os.path.join(root_path, 'calib/')
    if not os.path.exists(calib_path):
        os.makedirs(calib_path)

    sof_path = os.path.join(root_path, 'sof/')
    if not os.path.exists(sof_path):
        os.makedirs(sof_path)
        
    tmp_path = os.path.join(root_path, 'tmp/')
    if not os.path.exists(tmp_path):
        os.makedirs(tmp_path)
        
    # get list of files
    raw = calibs_info[np.logical_not(calibs_info['PROCESSED'])]
    calibs = raw[(calibs_info['DPR TYPE'] == 'DARK') | (calibs_info['DPR TYPE'] == 'DARK,BACKGROUND') |
                 (calibs_info['DPR TYPE'] == 'SKY')]

    # loops on type and DIT value
    types = ['DARK', 'DARK,BACKGROUND', 'SKY']
    DITs  = calibs['DET SEQ1 DIT'].unique().round(2)

    for ctype in types:
        for DIT in DITs:            
            cfiles = calibs[(calibs['DPR TYPE'] == ctype) & (calibs['DET SEQ1 DIT'].round(2) == DIT)]
            files = cfiles.index

            # skip non-existing combinations
            if len(cfiles) == 0:
                continue
            
            print(' * {0} with DIT={1:.2f} sec ({2} files)'.format(ctype, DIT, len(cfiles)))

            # create sof
            sof = os.path.join(sof_path, 'dark_DIT={0:.2f}.sof'.format(DIT))
            file = open(sof, 'w')
            for f in files:
                file.write('{0}{1}.fits     {2}\n'.format(raw_path, f, 'IFS_DARK_RAW'))
            file.close()

            # products
            dark_file = 'dark_DIT={0:.2f}'.format(DIT)
            bpm_file  = 'dark_bpm_DIT={0:.2f}'.format(DIT)

            # execute esorex    
            args = ['esorex',
                    '--msg-level=debug',
                    'sph_ifs_master_dark',
                    '--ifs.master_dark.coll_alg=2',
                    '--ifs.master_dark.sigma_clip=5.0',
                    '--ifs.master_dark.smoothing=5',
                    '--ifs.master_dark.min_acceptable=0.0',
                    '--ifs.master_dark.max_acceptable=2000.0',
                    '--ifs.master_dark.outfilename={0}{1}.fits'.format(calib_path, dark_file),
                    '--ifs.master_dark.badpixfilename={0}{1}.fits'.format(calib_path, bpm_file),
                    sof]

            if silent:
                proc = subprocess.run(args, cwd=tmp_path, stdout=subprocess.DEVNULL)
            else:
                proc = subprocess.run(args, cwd=tmp_path)

            if proc.returncode != 0:
                raise ValueError('esorex process was not successful')

            # store products
            calibs_info.loc[dark_file, 'DPR CATG'] = cfiles['DPR CATG'][0]
            calibs_info.loc[dark_file, 'DPR TYPE'] = cfiles['DPR TYPE'][0]
            calibs_info.loc[dark_file, 'INS2 COMB IFS'] = cfiles['INS2 COMB IFS'][0]
            calibs_info.loc[dark_file, 'DET SEQ1 DIT'] = cfiles['DET SEQ1 DIT'][0]
            calibs_info.loc[dark_file, 'PROCESSED'] = True
            calibs_info.loc[dark_file, 'PRO CATG'] = 'IFS_MASTER_DARK'

            calibs_info.loc[bpm_file, 'DPR CATG'] = cfiles['DPR CATG'][0]
            calibs_info.loc[bpm_file, 'DPR TYPE'] = cfiles['DPR TYPE'][0]
            calibs_info.loc[bpm_file, 'INS2 COMB IFS'] = cfiles['INS2 COMB IFS'][0]
            calibs_info.loc[bpm_file, 'PROCESSED'] = True
            calibs_info.loc[bpm_file, 'PRO CATG']  = 'IFS_STATIC_BADPIXELMAP'

    # save
    calibs_info.to_csv(os.path.join(root_path, 'calibs.csv'))


def compute_bad_pixel_map(bpm_files, dtype=np.uint8):
    '''
    Compute a combined bad pixel map provided a list of files

    Parameters
    ----------
    bpm_files : list
        List of names for the bpm files

    dtype : data type
        Data type for the final bpm
    
    Returns
    bpm : array_like
        Combined bad pixel map
    '''

    # star with empty bpm
    bpm = np.zeros((2048, 2048), dtype=np.uint8)

    # fill if files are provided
    if len(bpm_files) != 0:
        for f in bpm_files:
            data = fits.getdata(f)
            bpm = np.logical_or(bpm, data)

    bpm = bpm.astype(dtype)

    return bpm


def compute_detector_flat(raw_flat_files, bpm_files=[]):
    '''
    Compute a master detector flat and associated bad pixel map

    Parameters
    ----------
    raw_flat_files : list
        List of 2 raw flat files

    bpm_files : list
        List of bad pixel map files

    Returns
    -------
    flat : array
        Master detector flat

    bpm : array
        Bad pixel map from flat
    '''

    # read bad pixel maps
    bpm_in = compute_bad_pixel_map(bpm_files, dtype=np.uint8)

    # read data
    ff0, hdr0 = fits.getdata(raw_flat_files[0], header=True)
    ff1, hdr1 = fits.getdata(raw_flat_files[1], header=True)

    # flatten if needed
    if ff0.ndim == 3:
        ff0 = np.median(ff0, axis=0)

    if ff1.ndim == 3:
        ff1 = np.median(ff1, axis=0)

    # create master flat
    DIT0 = hdr0['HIERARCH ESO DET SEQ1 DIT']
    DIT1 = hdr1['HIERARCH ESO DET SEQ1 DIT']
    
    if DIT0 > DIT1:
        flat = ff0 - ff1
    else:
        flat = ff1 - ff0
    
    # bad pixels correction    
    flat = imutils.fix_badpix(flat, bpm_in, box=5)    
    flat = imutils.sigma_filter(flat, box=5, nsigma=3, iterate=True)
    flat = imutils.sigma_filter(flat, box=7, nsigma=3, iterate=True)

    # normalized flat
    flat = flat / np.median(flat)

    # additional rounad of bad pixels correction
    bpm = (flat <= 0.9) | (flat >= 1.1)
    bpm = bpm.astype(np.uint8)
    flat = imutils.fix_badpix(flat, bpm, box=5)

    # final products
    flat = flat / np.median(flat)
    
    bpm = (flat <= 0.9) | (flat >= 1.1)
    bpm = bpm.astype(np.uint8)
    
    return flat, bpm

    
def sph_ifs_cal_detector_flat(root_path, calibs_info, silent=True):
    '''
    Create the detector flat calibrations

    Parameters
    ----------
    root_path : str
        Path to the dataset

    calibs_info : dataframe
        The data frame with all the information on calibration files

    silent : bool
        Suppress esorex output. Optional, default is True
    '''

    print('Creating flats')

    # check directories
    raw_path = os.path.join(root_path, 'raw/')
    
    calib_path = os.path.join(root_path, 'calib/')
    if not os.path.exists(calib_path):
        os.makedirs(calib_path)

    sof_path = os.path.join(root_path, 'sof/')
    if not os.path.exists(sof_path):
        os.makedirs(sof_path)

    tmp_path = os.path.join(root_path, 'tmp/')
    if not os.path.exists(tmp_path):
        os.makedirs(tmp_path)

    # get list of files
    raw = calibs_info[np.logical_not(calibs_info['PROCESSED'])]
    calibs = raw[(calibs_info['DPR TYPE'] == 'FLAT,LAMP') | (calibs_info['DPR TECH'] == 'IMAGE')]

    # IFS obs mode
    mode = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS2 COMB IFS'].unique()[0]    
    if mode == 'OBS_YJ':
        mode_short = 'YJ'
    elif mode == 'OBS_H':
        mode_short = 'YJH'
    else:
        raise ValueError('Unknown IFS mode {0}'.format(mode))

    # bpm files
    cfiles = calibs_info[calibs_info['PRO CATG'] == 'IFS_STATIC_BADPIXELMAP'].index
    bpm_files = [os.path.join(calib_path, f+'.fits') for f in cfiles]
    
    # loop on wavelengths
    waves = [         0,        1020,        1230,        1300,        1550]
    combs = ['CAL_BB_2', 'CAL_NB1_1', 'CAL_NB2_1', 'CAL_NB3_1', 'CAL_NB4_2']
    lamps = [         5,           1,           2,           3,           4]
    
    for wave, comb, lamp in zip(waves, combs, lamps):
        print(' * flat for wavelength {0} nm (filter {1}, lamp {2})'.format(wave, comb, lamp))
        
        cfiles = calibs[calibs['INS2 COMB IFS'] == '{0}_{1}'.format(comb, mode_short)]
        files = [os.path.join(raw_path, f+'.fits') for f in cfiles.index]

        if len(files) == 0:
            continue
        elif len(files) != 2:
            raise ValueError('There should be exactly 2 raw flat files. Found {0}.'.format(len(files)))

        # create the flat and bpm
        flat, bpm = compute_detector_flat(files, bpm_files=bpm_files)
    
        # products
        if wave == 0:
            wav = 'white'
        else:
            wav = str(int(wave))
        flat_file = 'master_detector_flat_{0}_l{1}'.format(wav, lamp)
        bpm_file  = 'dff_badpixels_{0}_l{1}'.format(wav, lamp)        

        hdu = fits.open(os.path.join(raw_path, files[0]))
        fits.writeto(os.path.join(calib_path, flat_file+'.fits'), flat, header=hdu[0].header, output_verify='silentfix', overwrite=True)
        fits.writeto(os.path.join(calib_path, bpm_file+'.fits'), bpm, header=hdu[0].header, output_verify='silentfix', overwrite=True)
        hdu.close()

        # store products
        calibs_info.loc[flat_file, 'DPR CATG'] = cfiles['DPR CATG'][0]
        calibs_info.loc[flat_file, 'DPR TYPE'] = cfiles['DPR TYPE'][0]
        calibs_info.loc[flat_file, 'INS2 COMB IFS'] = cfiles['INS2 COMB IFS'][0]
        calibs_info.loc[flat_file, 'DET SEQ1 DIT'] = cfiles['DET SEQ1 DIT'][0]
        calibs_info.loc[flat_file, 'PROCESSED'] = True
        calibs_info.loc[flat_file, 'PRO CATG'] = 'IFS_MASTER_DFF'

        calibs_info.loc[bpm_file, 'DPR CATG'] = cfiles['DPR CATG'][0]
        calibs_info.loc[bpm_file, 'DPR TYPE'] = cfiles['DPR TYPE'][0]
        calibs_info.loc[bpm_file, 'INS2 COMB IFS'] = cfiles['INS2 COMB IFS'][0]
        calibs_info.loc[bpm_file, 'PROCESSED'] = True
        calibs_info.loc[bpm_file, 'PRO CATG']  = 'IFS_STATIC_BADPIXELMAP'

    # save
    calibs_info.to_csv(os.path.join(root_path, 'calibs.csv'))


def sph_ifs_cal_specpos(root_path, calibs_info, silent=True):
    '''
    Create the specpos calibration

    Parameters
    ----------
    root_path : str
        Path to the dataset

    calibs_info : dataframe
        The data frame with all the information on calibration files

    silent : bool
        Suppress esorex output. Optional, default is True
    '''

    print('Creating specpos')

    # check directories
    raw_path = os.path.join(root_path, 'raw/')
    
    calib_path = os.path.join(root_path, 'calib/')
    if not os.path.exists(calib_path):
        os.makedirs(calib_path)

    sof_path = os.path.join(root_path, 'sof/')
    if not os.path.exists(sof_path):
        os.makedirs(sof_path)

    tmp_path = os.path.join(root_path, 'tmp/')
    if not os.path.exists(tmp_path):
        os.makedirs(tmp_path)
    
    # get list of files
    specpos_file = calibs_info[np.logical_not(calibs_info['PROCESSED']) & (calibs_info['DPR TYPE'] == 'SPECPOS,LAMP')]
    if len(specpos_file) != 1:
        raise ValueError('There should be exactly 1 raw specpos files. Found {0}.'.format(len(specpos_file)))
    
    dark_file = calibs_info[calibs_info['PROCESSED'] & (calibs_info['PRO CATG'] == 'IFS_MASTER_DARK') & 
                            (calibs_info['DPR CATG'] == 'CALIB') & (calibs_info['DET SEQ1 DIT'].round(2) == 1.65)]
    if len(dark_file) == 0:
        raise ValueError('There should at least 1 dark file for calibrations. Found none.')

    # IFS obs mode
    mode = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS2 COMB IFS'].unique()[0]    
    if mode == 'OBS_YJ':
        Hmode = 'FALSE'
    elif mode == 'OBS_H':
        Hmode = 'TRUE'
    else:
        raise ValueError('Unknown IFS mode {0}'.format(mode))

    # create sof
    sof = os.path.join(sof_path, 'specpos.sof')
    file = open(sof, 'w')
    file.write('{0}{1}.fits     {2}\n'.format(raw_path, specpos_file.index[0], 'IFS_SPECPOS_RAW'))
    file.write('{0}{1}.fits     {2}\n'.format(calib_path, dark_file.index[0], 'IFS_MASTER_DARK'))
    file.close()
    
    # products
    specp_file = 'spectra_positions'

    # execute esorex    
    args = ['esorex',
            'sph_ifs_spectra_positions',
            '--ifs.spectra_positions.hmode={0}'.format(Hmode),
            '--ifs.spectra_positions.outfilename={0}{1}.fits'.format(calib_path, specp_file),
            sof]

    if silent:
        proc = subprocess.run(args, cwd=tmp_path, stdout=subprocess.DEVNULL)
    else:
        proc = subprocess.run(args, cwd=tmp_path)

    if proc.returncode != 0:
        raise ValueError('esorex process was not successful')

    # store products
    calibs_info.loc[specp_file, 'DPR CATG'] = specpos_file['DPR CATG'][0]
    calibs_info.loc[specp_file, 'DPR TYPE'] = specpos_file['DPR TYPE'][0]
    calibs_info.loc[specp_file, 'INS2 COMB IFS'] = specpos_file['INS2 COMB IFS'][0]
    calibs_info.loc[specp_file, 'DET SEQ1 DIT'] = specpos_file['DET SEQ1 DIT'][0]
    calibs_info.loc[specp_file, 'PROCESSED'] = True
    calibs_info.loc[specp_file, 'PRO CATG'] = 'IFS_SPECPOS'

    # save
    calibs_info.to_csv(os.path.join(root_path, 'calibs.csv'))


def sph_ifs_cal_wave(root_path, calibs_info, silent=True):
    '''
    Create the wavelength calibration

    Parameters
    ----------
    root_path : str
        Path to the dataset

    calibs_info : dataframe
        The data frame with all the information on calibration files

    silent : bool
        Suppress esorex output. Optional, default is True
    '''

    print('Creating wavelength calibration')

    # check directories
    raw_path = os.path.join(root_path, 'raw/')
    
    calib_path = os.path.join(root_path, 'calib/')
    if not os.path.exists(calib_path):
        os.makedirs(calib_path)

    sof_path = os.path.join(root_path, 'sof/')
    if not os.path.exists(sof_path):
        os.makedirs(sof_path)
    
    tmp_path = os.path.join(root_path, 'tmp/')
    if not os.path.exists(tmp_path):
        os.makedirs(tmp_path)
        
    # get list of files
    wave_file = calibs_info[np.logical_not(calibs_info['PROCESSED']) & (calibs_info['DPR TYPE'] == 'WAVE,LAMP')]
    if len(wave_file) != 1:
        raise ValueError('There should be exactly 1 raw wavelength calibration file. Found {0}.'.format(len(wave_file)))
    
    specpos_file = calibs_info[calibs_info['PROCESSED'] & (calibs_info['PRO CATG'] == 'IFS_SPECPOS')]
    if len(specpos_file) != 1:
        raise ValueError('There should be exactly 1 specpos file. Found {0}.'.format(len(specpos_file)))
    
    dark_file = calibs_info[calibs_info['PROCESSED'] & (calibs_info['PRO CATG'] == 'IFS_MASTER_DARK') & 
                            (calibs_info['DPR CATG'] == 'CALIB') & (calibs_info['DET SEQ1 DIT'].round(2) == 1.65)]
    if len(dark_file) == 0:
        raise ValueError('There should at least 1 dark file for calibrations. Found none.')

    # IFS obs mode
    mode = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS2 COMB IFS'].unique()[0]            
        
    # create sof
    sof = os.path.join(sof_path, 'wave.sof')
    file = open(sof, 'w')
    file.write('{0}{1}.fits     {2}\n'.format(raw_path, wave_file.index[0], 'IFS_WAVECALIB_RAW'))
    file.write('{0}{1}.fits     {2}\n'.format(calib_path, specpos_file.index[0], 'IFS_SPECPOS'))
    file.write('{0}{1}.fits     {2}\n'.format(calib_path, dark_file.index[0], 'IFS_MASTER_DARK'))
    file.close()
    
    # products
    wav_file = 'wave_calib'

    # execute esorex
    if mode == 'OBS_YJ':
        args = ['esorex',
                'sph_ifs_wave_calib',
                '--ifs.wave_calib.number_lines=3',
                '--ifs.wave_calib.wavelength_line1=0.9877',
                '--ifs.wave_calib.wavelength_line2=1.1237',
                '--ifs.wave_calib.wavelength_line3=1.3094',
                '--ifs.wave_calib.outfilename={0}{1}.fits'.format(calib_path, wav_file),
                sof]
    elif mode == 'OBS_YJH':
        args = ['esorex',
                'sph_ifs_wave_calib',
                '--ifs.wave_calib.number_lines=3',
                '--ifs.wave_calib.wavelength_line1=0.9877',
                '--ifs.wave_calib.wavelength_line2=1.1237',
                '--ifs.wave_calib.wavelength_line3=1.3094',
                '--ifs.wave_calib.wavelength_line4=1.5451',
                '--ifs.wave_calib.outfilename={0}{1}.fits'.format(calib_path, wav_file),
                sof]

    if silent:
        proc = subprocess.run(args, cwd=tmp_path, stdout=subprocess.DEVNULL)
    else:
        proc = subprocess.run(args, cwd=tmp_path)

    if proc.returncode != 0:
        raise ValueError('esorex process was not successful')

    # store products
    calibs_info.loc[wav_file, 'DPR CATG'] = wave_file['DPR CATG'][0]
    calibs_info.loc[wav_file, 'DPR TYPE'] = wave_file['DPR TYPE'][0]
    calibs_info.loc[wav_file, 'INS2 COMB IFS'] = wave_file['INS2 COMB IFS'][0]
    calibs_info.loc[wav_file, 'DET SEQ1 DIT'] = wave_file['DET SEQ1 DIT'][0]
    calibs_info.loc[wav_file, 'PROCESSED'] = True
    calibs_info.loc[wav_file, 'PRO CATG'] = 'IFS_WAVECALIB'

    # save
    calibs_info.to_csv(os.path.join(root_path, 'calibs.csv'))

    
def sph_ifs_cal_ifu_flat(root_path, calibs_info, silent=True):
    '''
    Create the IFU flat calibration

    Parameters
    ----------
    root_path : str
        Path to the dataset

    calibs_info : dataframe
        The data frame with all the information on calibration files

    silent : bool
        Suppress esorex output. Optional, default is True
    '''

    print('Creating IFU flat')

    # check directories
    raw_path = os.path.join(root_path, 'raw/')
    
    calib_path = os.path.join(root_path, 'calib/')
    if not os.path.exists(calib_path):
        os.makedirs(calib_path)

    sof_path = os.path.join(root_path, 'sof/')
    if not os.path.exists(sof_path):
        os.makedirs(sof_path)

    tmp_path = os.path.join(root_path, 'tmp/')
    if not os.path.exists(tmp_path):
        os.makedirs(tmp_path)
        
    # IFS obs mode
    mode = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'INS2 COMB IFS'].unique()[0]            
    if mode == 'OBS_YJ':
        mode_short = 'YJ'
    elif mode == 'OBS_H':
        mode_short = 'YJH'
    else:
        raise ValueError('Unknown IFS mode {0}'.format(mode))

    # get list of files
    ifu_flat_file = calibs_info[np.logical_not(calibs_info['PROCESSED']) & (calibs_info['DPR TYPE'] == 'FLAT,LAMP') &
                                (calibs_info['DPR TECH'] == 'IFU')]
    if len(ifu_flat_file) != 1:
        raise ValueError('There should be exactly 1 raw IFU flat file. Found {0}.'.format(len(ifu_flat_file)))
    
    wave_file = calibs_info[calibs_info['PROCESSED'] & (calibs_info['PRO CATG'] == 'IFS_WAVECALIB')]
    if len(wave_file) != 1:
        raise ValueError('There should be exactly 1 wavelength calibration file. Found {0}.'.format(len(wave_file)))
    
    dark_file = calibs_info[calibs_info['PROCESSED'] & (calibs_info['PRO CATG'] == 'IFS_MASTER_DARK') &
                            (calibs_info['DPR CATG'] == 'CALIB') & (calibs_info['DET SEQ1 DIT'].round(2) == 1.65)]
    if len(dark_file) == 0:
        raise ValueError('There should at least 1 dark file for calibrations. Found none.')

    flat_white_file = calibs_info[calibs_info['PROCESSED'] & (calibs_info['PRO CATG'] == 'IFS_MASTER_DFF') &
                                  (calibs_info['INS2 COMB IFS'] == 'CAL_BB_2_{0}'.format(mode_short))]
    if len(flat_white_file) != 1:
        raise ValueError('There should be exactly 1 white flat file. Found {0}.'.format(len(flat_white_file)))
    
    flat_1020_file = calibs_info[calibs_info['PROCESSED'] & (calibs_info['PRO CATG'] == 'IFS_MASTER_DFF') &
                                 (calibs_info['INS2 COMB IFS'] == 'CAL_NB1_1_{0}'.format(mode_short))]
    if len(flat_1020_file) != 1:
        raise ValueError('There should be exactly 1 1020 nm flat file. Found {0}.'.format(len(flat_1020_file)))
    
    flat_1230_file = calibs_info[calibs_info['PROCESSED'] & (calibs_info['PRO CATG'] == 'IFS_MASTER_DFF') &
                                 (calibs_info['INS2 COMB IFS'] == 'CAL_NB2_1_{0}'.format(mode_short))]
    if len(flat_1230_file) != 1:
        raise ValueError('There should be exactly 1 1230 nm flat file. Found {0}.'.format(len(flat_1230_file)))

    flat_1300_file = calibs_info[calibs_info['PROCESSED'] & (calibs_info['PRO CATG'] == 'IFS_MASTER_DFF') &
                                 (calibs_info['INS2 COMB IFS'] == 'CAL_NB3_1_{0}'.format(mode_short))]
    if len(flat_1300_file) != 1:
        raise ValueError('There should be exactly 1 1300 nm flat file. Found {0}.'.format(len(flat_1300_file)))
    
    if mode == 'OBS_YJH':
        flat_1550_file = calibs_info[calibs_info['PROCESSED'] & (calibs_info['PRO CATG'] == 'IFS_MASTER_DFF') &
                                     (calibs_info['INS2 COMB IFS'] == 'CAL_NB4_2_{0}'.format(mode_short))]
        if len(flat_1550_file) != 1:
            raise ValueError('There should be exactly 1 1550 nm flat file. Found {0}.'.format(len(flat_1550_file)))
    
    # create sof
    sof = os.path.join(sof_path, 'ifu_flat.sof')
    file = open(sof, 'w')
    file.write('{0}{1}.fits     {2}\n'.format(raw_path, ifu_flat_file.index[0], 'IFS_FLAT_FIELD_RAW'))
    file.write('{0}{1}.fits     {2}\n'.format(calib_path, wave_file.index[0], 'IFS_WAVECALIB'))
    file.write('{0}{1}.fits     {2}\n'.format(calib_path, dark_file.index[0], 'IFS_MASTER_DARK'))
    file.write('{0}{1}.fits     {2}\n'.format(calib_path, flat_white_file.index[0], 'IFS_MASTER_DFF_SHORT'))
    file.write('{0}{1}.fits     {2}\n'.format(calib_path, flat_white_file.index[0], 'IFS_MASTER_DFF_LONGBB'))
    file.write('{0}{1}.fits     {2}\n'.format(calib_path, flat_1020_file.index[0], 'IFS_MASTER_DFF_LONG1'))
    file.write('{0}{1}.fits     {2}\n'.format(calib_path, flat_1230_file.index[0], 'IFS_MASTER_DFF_LONG2'))
    file.write('{0}{1}.fits     {2}\n'.format(calib_path, flat_1300_file.index[0], 'IFS_MASTER_DFF_LONG3'))
    if mode == 'OBS_YJH':
        file.write('{0}{1}     {2}\n'.format(calib_path, flat_1550_file.index[0], 'IFS_MASTER_DFF_LONG4'))
    file.close()

    # products
    ifu_file = 'ifu_flat'

    # execute esorex
    args = ['esorex',
            'sph_ifs_instrument_flat',
            '--ifs.instrument_flat.nofit=TRUE',
            '--ifs.instrument_flat.ifu_filename={0}{1}.fits'.format(calib_path, ifu_file),
            sof]

    if silent:
        proc = subprocess.run(args, cwd=tmp_path, stdout=subprocess.DEVNULL)
    else:
        proc = subprocess.run(args, cwd=tmp_path)

    if proc.returncode != 0:
        raise ValueError('esorex process was not successful')

    # store products
    calibs_info.loc[ifu_file, 'DPR CATG'] = ifu_flat_file['DPR CATG'][0]
    calibs_info.loc[ifu_file, 'DPR TYPE'] = ifu_flat_file['DPR TYPE'][0]
    calibs_info.loc[ifu_file, 'INS2 COMB IFS'] = ifu_flat_file['INS2 COMB IFS'][0]
    calibs_info.loc[ifu_file, 'DET SEQ1 DIT'] = ifu_flat_file['DET SEQ1 DIT'][0]
    calibs_info.loc[ifu_file, 'PROCESSED'] = True
    calibs_info.loc[ifu_file, 'PRO CATG'] = 'IFS_IFU_FLAT_FIELD'

    # save
    calibs_info.to_csv(os.path.join(root_path, 'calibs.csv'))


def sph_ifs_correct_spectral_xtalk(img):
    '''
    Corrects a IFS frame from the spectral crosstalk

    This routines corrects for the SPHERE/IFS spectral crosstalk at
    small scales and (optionally) at large scales. This correction is
    necessary to correct the signal that is "leaking" between
    lenslets. See Antichi et al. (2009ApJ...695.1042A) for a
    theoretical description of the IFS crosstalk. Some informations
    regarding its correction are provided in Vigan et al. (2015), but
    this procedure still lacks a rigorous description and performance
    analysis.

    Since the correction of the crosstalk involves a convolution by a
    kernel of size 41x41, the values at the edges of the frame depend
    on how you choose to apply the convolution. Current implementation
    is EDGE_TRUNCATE. In other parts of the image (i.e. far from the
    edges), the result is identical to original routine by Dino
    Mesa. Note that in the original routine, the convolution that was
    coded did not treat the edges in a clean way defined
    mathematically. The scipy.ndimage.convolve() function offers
    different possibilities for the edges that are all documented.
    
    Parameters
    ----------
    img : array_like
        Input IFS science frame

    Returns
    -------
    img_corr : array_like
        Science frame corrected from the spectral crosstalk

    '''
    
    # definition of the dimension of the matrix
    sepmax = 20
    dim    = sepmax*2+1
    bfac   = 0.727986/1.8

    # defines a matrix to be used around each pixel
    # (the value of the matrix is lower for greater
    # distances form the center.
    x, y = np.meshgrid(np.arange(dim)-sepmax, np.arange(dim)-sepmax)
    rdist  = np.sqrt(x**2 + y**2)
    kernel = 1 / (1+rdist**3 / bfac**3)
    kernel[(np.abs(x) <= 1) & (np.abs(y) <= 1)] = 0

    # convolution and subtraction
    conv = ndimage.convolve(img, kernel, mode='reflect')
    img_corr = img - conv

    return img_corr


def collapse_frames_info(finfo, fname, collapse_type, coadd_value=2):
    '''
    Collapse frame info to match the collapse operated on the data

    Parameters
    ----------
    finfo : dataframe
        The data frame with all the information on science frames

    fname : str
       The name of the current file
    
    collapse_type : str
        Type of collapse. Possible values are mean or coadd. Default
        is mean.

    coadd_value : int
        Number of consecutive frames to be coadded when collapse_type
        is coadd. Default is 2

    Returns
    -------
    nfinfo : dataframe
        Collapsed data frame
    '''
    
    print('   ==> collapse frames information')

    nfinfo = None
    if collapse_type == 'none':
        nfinfo = finfo
    elif collapse_type == 'mean':
        index = pd.MultiIndex.from_arrays([[fname], [0]], names=['FILE', 'IMG'])
        nfinfo = pd.DataFrame(columns=finfo.columns, index=index)

        # get min/max indices
        imin = finfo.index.get_level_values(1).min()
        imax = finfo.index.get_level_values(1).max()

        # copy data
        nfinfo.loc[(fname, 0)] = finfo.loc[(fname, imin)]
        
        # update time values
        nfinfo.loc[(fname, 0), 'DET NDIT'] = 1
        nfinfo.loc[(fname, 0), 'TIME START'] = finfo.loc[(fname, imin), 'TIME START']
        nfinfo.loc[(fname, 0), 'TIME END'] = finfo.loc[(fname, imax), 'TIME END']
        nfinfo.loc[(fname, 0), 'TIME'] = finfo.loc[(fname, imin), 'TIME START'] + \
                                         (finfo.loc[(fname, imax), 'TIME END'] - finfo.loc[(fname, imin), 'TIME START']) / 2
        
        # recompute angles
        compute_angles(nfinfo)            
    elif collapse_type == 'coadd':
        coadd_value = int(coadd_value)
        NDIT = len(finfo)
        NDIT_new = NDIT // coadd_value

        index = pd.MultiIndex.from_arrays([np.full(NDIT_new, fname), np.arange(NDIT_new)], names=['FILE', 'IMG'])
        nfinfo = pd.DataFrame(columns=finfo.columns, index=index)

        for f in range(NDIT_new):
            # get min/max indices
            imin = int(f*coadd_value)
            imax = int((f+1)*coadd_value-1)

            # copy data
            nfinfo.loc[(fname, f)] = finfo.loc[(fname, imin)]

            # update time values
            nfinfo.loc[(fname, f), 'DET NDIT'] = 1
            nfinfo.loc[(fname, f), 'TIME START'] = finfo.loc[(fname, imin), 'TIME START']
            nfinfo.loc[(fname, f), 'TIME END'] = finfo.loc[(fname, imax), 'TIME END']
            nfinfo.loc[(fname, f), 'TIME'] = finfo.loc[(fname, imin), 'TIME START'] + \
                                             (finfo.loc[(fname, imax), 'TIME END'] - finfo.loc[(fname, imin), 'TIME START']) / 2

        # recompute angles
        compute_angles(nfinfo)
    else:
        raise ValueError('Unknown collapse type {0}'.format(collapse_type))        

    return nfinfo
    

    
def sph_ifs_preprocess(root_path, files_info, calibs_info, frames_info,
                       subtract_background=True, fix_badpix=True, correct_xtalk=True,
                       collapse_science=False, collapse_type='mean', coadd_value=2,
                       collapse_psf=True, collapse_center=True):
    '''Pre-processes the science frames.

    This function can perform multiple steps:
      - collapse of the frames according to different schemes
      - subtract the background
      - correct bad pixels
      - correct the spectral crosstalk

    For the science, 2 collapse methods are available: mean or
    coadd. With mean, the full cubes are mean-combined into a single
    frame. With coadd, the frames are coadded following the
    coadd_value. This can result in lost frames if the number of NDIT
    is not a multiple of coadd_value.

    For the PSFs and star center frames, there is either no collapse
    or a mean collapse.
    
    Parameters
    ----------
    root_path : str
        Path to the dataset

    files_info : dataframe
        The data frame with all the information on raw science files

    calibs_info : dataframe
        The data frame with all the information on calibration files

    frames_info : dataframe
        The data frame with all the information on science frames

    subtract_background : bool
        Performs background subtraction. Default is True

    fix_badpix : bool
        Performs correction of bad pixels. Default is True

    correct_xtalk : bool
        Performs spectral crosstalk correction. Default is True

    collapse_science :  bool
        Collapse data for OBJECT cubes. Default is False

    collapse_type : str
        Type of collapse. Possible values are mean or coadd. Default
        is mean.

    coadd_value : int
        Number of consecutive frames to be coadded when collapse_type
        is coadd. Default is 2

    collapse_psf :  bool
        Collapse data for OBJECT,FLUX cubes. Default is True. Note
        that the collapse type is mean and cannot be changed.
    
    collapse_center :  bool
        Collapse data for OBJECT,CENTER cubes. Default is True. Note
        that the collapse type is mean and cannot be changed.    
    '''

    print('Pre-processing science files')

    # check directories
    raw_path = os.path.join(root_path, 'raw/')
    
    calib_path = os.path.join(root_path, 'calib/')
    if not os.path.exists(calib_path):
        os.makedirs(calib_path)

    preproc_path = os.path.join(root_path, 'preproc/')
    if not os.path.exists(preproc_path):
        os.makedirs(preproc_path)

    # bpm
    if fix_badpix:
        bpm_files = calibs_info[calibs_info['PRO CATG'] == 'IFS_STATIC_BADPIXELMAP'].index
        bpm_files = [os.path.join(calib_path, f+'.fits') for f in bpm_files]

        bpm = compute_bad_pixel_map(bpm_files)

    # final dataframe
    index = pd.MultiIndex(names=['FILE', 'IMG'], levels=[[], []], labels=[[], []])
    frames_info_collapse = pd.DataFrame(index=index, columns=frames_info.columns)
        
    # loop on the different type of science files
    sci_types = ['OBJECT,CENTER', 'OBJECT,FLUX', 'OBJECT']
    dark_types = ['SKY', 'DARK,BACKGROUND', 'DARK']
    for typ in sci_types:    
        # science files
        sci_files = files_info[(files_info['DPR CATG'] == 'SCIENCE') & (files_info['DPR TYPE'] == typ)]
        sci_DITs = list(sci_files['DET SEQ1 DIT'].round(2).unique())
        
        if len(sci_files) == 0:
            continue        

        for DIT in sci_DITs:
            sfiles = sci_files[sci_files['DET SEQ1 DIT'].round(2) == DIT]
            
            print('{0} files of type {1} with DIT={2} sec'.format(len(sfiles), typ, DIT))

            if subtract_background:
                # look for sky, then background, then darks
                # normally there should be only one with the proper DIT
                dfiles = []
                for d in dark_types:
                    dfiles = calibs_info[(calibs_info['PRO CATG'] == 'IFS_MASTER_DARK') &
                                         (calibs_info['DPR TYPE'] == d) & (calibs_info['DET SEQ1 DIT'].round(2) == DIT)]
                    if len(dfiles) != 0:
                        break
                print('   ==> found {0} corresponding {1} file'.format(len(dfiles), d))

                if len(dfiles) != 1:
                    raise ValueError('Unexpected number of background fliles ({0})'.format(len(dfiles)))

                bkg = fits.getdata(os.path.join(calib_path, dfiles.index[0]+'.fits'))
                
            # process files
            for fname, finfo in sci_files.iterrows():
                print(' * {0}'.format(fname))

                # frames_info extract
                finfo = frames_info.loc[(fname, slice(None)), :]
                
                # read data
                print('   ==> read data')
                img, hdr = fits.getdata(os.path.join(raw_path, fname+'.fits'), header=True)
                
                # add extra dimension to single images to make cubes
                if img.ndim == 2:
                    img = img[np.newaxis, ...]
                
                # collapse
                if (typ == 'OBJECT,CENTER'):
                    if collapse_center:
                        print('   ==> collapse: mean')
                        img = np.mean(img, axis=0, keepdims=True)
                        frames_info_new = collapse_frames_info(finfo, fname, 'mean')
                    else:
                        frames_info_new = collapse_frames_info(finfo, fname, 'none')
                elif (typ == 'OBJECT,FLUX'):
                    if collapse_psf:
                        print('   ==> collapse: mean')
                        img = np.mean(img, axis=0, keepdims=True)
                        frames_info_new = collapse_frames_info(finfo, fname, 'mean')
                    else:
                        frames_info_new = collapse_frames_info(finfo, fname, 'none')
                elif (typ == 'OBJECT'):
                    if collapse_science:                        
                        if collapse_type == 'mean':
                            print('   ==> collapse: mean ({0} -> 1 frame, 0 dropped)'.format(len(img)))
                            img = np.mean(img, axis=0, keepdims=True)

                            frames_info_new = collapse_frames_info(finfo, fname, 'mean')
                        elif collapse_type == 'coadd':
                            if (not isinstance(coadd_value, int)) or (coadd_value <= 1):
                                raise TypeError('coadd_value must be an integer >1')
                            
                            coadd_value = int(coadd_value)
                            NDIT = len(img)
                            NDIT_new = NDIT // coadd_value
                            dropped = NDIT % coadd_value

                            if coadd_value > NDIT:
                                raise ValueError('coadd_value ({0}) must be < NDIT ({1})'.format(coadd_value, NDIT))
                            
                            print('   ==> collapse: mean ({0} -> {1} frames, {2} dropped)'.format(NDIT, NDIT_new, dropped))

                            # coadd frames
                            nimg = np.empty((NDIT_new, 2048, 2048), dtype=img.dtype)
                            for f in range(NDIT_new):
                                nimg[f] = np.mean(img[f*coadd_value:(f+1)*coadd_value], axis=0)
                            img = nimg

                            frames_info_new = collapse_frames_info(finfo, fname, 'coadd', coadd_value=coadd_value)
                        else:
                            raise ValueError('Unknown collapse type {0}'.format(collapse_type))
                    else:
                        frames_info_new = collapse_frames_info(finfo, fname, 'none')

                # merge collapse collapsed frames_info
                frames_info_collapse = pd.concat((frames_info_collapse, frames_info_new))
                        
                # background subtraction
                if subtract_background:
                    print('   ==> subtract background')
                    for f in range(len(img)):
                        img[f] -= bkg

                # bad pixels correction
                if fix_badpix:
                    print('   ==> correct bad pixels')
                    for f in range(len(img)):
                        frame = img[f]
                        frame = imutils.fix_badpix(frame, bpm, box=5)
                        frame = imutils.sigma_filter(frame, box=5, nsigma=3, iterate=True)
                        frame = imutils.sigma_filter(frame, box=7, nsigma=3, iterate=True)
                        img[f] = frame

                # spectral crosstalk correction
                if correct_xtalk:
                    print('   ==> correct spectral crosstalk')
                    for f in range(len(img)):
                        frame = img[f]
                        frame = sph_ifs_correct_spectral_xtalk(frame)
                        img[f] = frame

                # save DITs individually
                for f in range(len(img)):
                    fits.writeto(os.path.join(preproc_path, fname+'_frame{0:03d}_preproc.fits'.format(f)), img, hdr, overwrite=True, output_verify='silentfix')
                                        
                print()
                
        print()

    # save final dataframe
    frames_info_collapse.to_csv(os.path.join(preproc_path, 'frames_preproc.csv'))
        


def clean(root_path):
    '''
    Clean everything exact raw data

    Parameters
    ----------
    root_path : str
        Path to the dataset
    '''
    pass
    
    
    
root_path = '/Users/avigan/data/pySPHERE-test/IFS/'

# files_info = sort_files(root_path)
# frames_info = sort_frames(root_path, files_info)
# calibs_info = files_association(root_path, files_info)

# files_info, calibs_info, frames_info, frames_info_collapse = read_info(root_path)
# sph_ifs_cal_dark(root_path, calibs_info)
# sph_ifs_cal_detector_flat(root_path, calibs_info)
# sph_ifs_cal_specpos(root_path, calibs_info)
# sph_ifs_cal_wave(root_path, calibs_info)
# sph_ifs_cal_ifu_flat(root_path, calibs_info)

files_info, calibs_info, frames_info, frames_info_collapse = read_info(root_path)
sph_ifs_preprocess(root_path, files_info, calibs_info, frames_info,
                   subtract_background=True, fix_badpix=False, correct_xtalk=False,
                   collapse_science=True, collapse_type='coadd', coadd_value=2,
                   collapse_psf=True, collapse_center=False)

# files_info, calibs_info, frames_info, frames_info_collapse = read_info(root_path)
