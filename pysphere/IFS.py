import os
import glob
import pandas as pd
import subprocess

from ast import literal_eval
from astropy.io import fits


# keywords to be saved
keywords = [
    # standard
    'INSTRUME', 
    'OBJECT', 'DATE-OBS', 'DATE', 'HIERARCH ESO DET FRAM UTC',

    # DPR
    'HIERARCH ESO DPR CATG', 'HIERARCH ESO DPR TYPE', 'HIERARCH ESO DPR TECH',
        
    # coordinates
    'HIERARCH ESO TEL GEOLAT', 'HIERARCH ESO TEL GEOLON',
    'HIERARCH ESO INS4 DROT1 RA', 'HIERARCH ESO INS4 DROT1 DEC',
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
    
    if len(files) == 0:
        raise ValueError('No raw FITS files in root_path directory')

    print(' * found {0} FITS files in {1}'.format(len(files), root_path))

    # files table
    files_info = pd.DataFrame(index=files, columns=keywords_short, dtype='float')
    
    for f in files:
        hdu = fits.open(f)
        hdr = hdu[0].header

        for k, sk in zip(keywords, keywords_short):
            files_info.loc[f, sk] = hdr.get(k)

        hdu.close()

    # save files_info
    files_info.to_csv(os.path.join(root_path, 'files.csv'))
    
    return files_info


def files_association(root_path, files_info):
    '''
    Performs the calibration files association and a sanity check

    Parameters
    ----------
    root_path : str
        Path to the dataset
    
    files_info : dataframe
        The data frame with all the information on raw files
    '''

    print('Performing file associtation for calibrations')
    
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

    # files associtation dictionary
    files_assoc = []
    
    ###############################################
    # static calibrations not dependent on science
    ###############################################

    calibs = files_info[files_info['DPR CATG'] == 'CALIB']
    
    # white flat
    cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == 'CAL_BB_2_{0}'.format(mode_short))]
    if len(cfiles) != 2:
        print(' * Error: there should be 2 flat files for white lamp!')
    files_assoc.append({'type': 'flat_white', 'DIT': 0, 'files': list(cfiles.index.values), 'products': []})
        
    # 1020 nm flat
    cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == 'CAL_NB1_1_{0}'.format(mode_short))]
    if len(cfiles) != 2:
        print(' * Error: there should be 2 flat files for 1020 nm filter!')
    files_assoc.append({'type': 'flat_1020', 'DIT': 0, 'files': list(cfiles.index.values), 'products': []})

    # 1230 nm flat
    cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == 'CAL_NB2_1_{0}'.format(mode_short))]
    if len(cfiles) != 2:
        print(' * Error: there should be 2 flat files for 1230 nm filter!')
    files_assoc.append({'type': 'flat_1230', 'DIT': 0, 'files': list(cfiles.index.values), 'products': []})

    # 1300 nm flat
    cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == 'CAL_NB3_1_{0}'.format(mode_short))]
    if len(cfiles) != 2:
        print(' * Error: there should be 2 flat files for 1300 nm filter!')
    files_assoc.append({'type': 'flat_1300', 'DIT': 0, 'files': list(cfiles.index.values), 'products': []})

    # 1550 nm flat (YJH mode only)
    if mode_short == 'YJH':
        cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == 'CAL_NB4_2_{0}'.format(mode_short))]
        if len(cfiles) != 2:
            print(' * Error: there should be 2 flat files for 1550 nm filter!')
        files_assoc.append({'type': 'flat_1550', 'DIT': 0, 'files': list(cfiles.index.values), 'products': []})

    # spectra position
    cfiles = calibs[(calibs['DPR TYPE'] == 'SPECPOS,LAMP') & (calibs['INS2 COMB IFS'] == mode)]
    if len(cfiles) != 1:
        print(' * Error: there should be 1 spectra position file!')
    files_assoc.append({'type': 'specpos', 'DIT': 0, 'files': list(cfiles.index.values), 'products': []})

    # wavelength
    cfiles = calibs[(calibs['DPR TYPE'] == 'WAVE,LAMP') & (calibs['INS2 COMB IFS'] == mode)]
    if len(cfiles) != 1:
        print(' * Error: there should be 1 wavelength calibration file!')
    files_assoc.append({'type': 'wave', 'DIT': 0, 'files': list(cfiles.index.values), 'products': []})
    
    # IFU flat
    cfiles = calibs[(calibs['DPR TYPE'] == 'FLAT,LAMP') & (calibs['INS2 COMB IFS'] == mode)]
    if len(cfiles) != 1:
        print(' * Error: there should be 1 IFU flat file!')
    files_assoc.append({'type': 'flat_ifu', 'DIT': 0, 'files': list(cfiles.index.values), 'products': []})
    
    # calibs dark file
    cfiles = calibs[((calibs['DPR TYPE'] == 'DARK') | (calibs['DPR TYPE'] == 'DARK,BACKGROUND')) &
                    (calibs['DET SEQ1 DIT'].round(2) == 1.65)]
    if len(cfiles) != 1:
        print(' * Warning: there is no dark/background for the basic calibrations (DIT=1.6 sec). ' +
              'It is *highly recommended* to include one to obtain the best data reduction. ' +
              'A single dark/background file is sufficient, and it can easily be downloaded ' +
              'from the ESO archive')
    files_assoc.append({'type': 'dark_calibs', 'DIT': 0, 'files': list(cfiles.index.values), 'products': []})
    
    ##################################################
    # static calibrations that depend on science (DIT)
    ##################################################
    
    DITs = files_info.loc[files_info['DPR CATG'] == 'SCIENCE', 'DET SEQ1 DIT'].unique().round(2)

    # handle darks in a slightly different way because there might be several different DITs
    for DIT in DITs:
        # instrumental backgrounds
        cfiles = calibs[((calibs['DPR TYPE'] == 'DARK') | (calibs['DPR TYPE'] == 'DARK,BACKGROUND')) &
                        (calibs['DET SEQ1 DIT'].round(2) == DIT)]
        if len(cfiles) == 0:
            print(' * Warning: there is no dark/background for science files with DIT={0} sec. '.format(DIT) +
                  'it is *highly recommended* to include one to obtain the best data reduction. ' +
                  'A single dark/background file is sufficient, and it can easily be downloaded ' +
                  'from the ESO archive')
        else:
            files_assoc.append({'type': 'background', 'DIT': DIT, 'files': list(cfiles.index.values), 'products': []})
            
        # sky backgrounds
        cfiles = files_info[(files_info['DPR TYPE'] == 'SKY') & (files_info['DET SEQ1 DIT'].round(2) == DIT)]
        if len(cfiles) == 0:
            print(' * Warning: there is no sky background for science files with DIT={0} sec. '.format(DIT) +
                  'Using a sky background instead of an internal instrumental background can ' +
                  'usually provide a cleaner data reduction')
        else:
            files_assoc.append({'type': 'sky', 'DIT': DIT, 'files': list(cfiles.index.values), 'products': []})

    # transform array into dataframe
    files_assoc = pd.DataFrame(files_assoc)
    
    # save
    files_assoc.to_json(os.path.join(root_path, 'calibs.json'))
    files_assoc.to_csv(os.path.join(root_path, 'calibs.csv'))
    
    return files_assoc


def sph_ifs_dark(root_path, calibs_info):
    '''
    Create the static backgrounds

    Parameters
    ----------
    root_path : str
        Path to the dataset

    calibs_info : dataframe
        The data frame with all the information on raw calibration files
    '''

    # check directories
    calib_path = os.path.join(root_path, 'calib/')
    if not os.path.exists(calib_path):
        os.makedirs(calib_path)

    sof_path = os.path.join(root_path, 'sof/')
    if not os.path.exists(sof_path):
        os.makedirs(sof_path)

    print('Creating dark and backgrounds')
        
    ###################################
    # calibs dark
    ###################################

    print(' * dark for calibrations')
    
    # get files
    ctype = 'dark_calibs'
    cfiles = calibs_info.loc[calibs_info['type'] == ctype, 'files']
    files = cfiles.values[0]
    
    # create sof
    sof = os.path.join(sof_path, '{0}.sof'.format(ctype))
    file = open(sof, 'w')
    for f in files:
        file.write('{0}     {1}\n'.format(f, 'IFS_DARK_RAW'))
    file.close()

    # products
    dark_file = '{0}dark_{1}.fits'.format(calib_path, ctype)
    bpm_file  = '{0}bpm_{1}.fits'.format(calib_path, ctype)
    
    # execute esorex    
    args = ['esorex',
            '--msg-level=debug',
            'sph_ifs_master_dark',
            '--ifs.master_dark.coll_alg=2',
            '--ifs.master_dark.sigma_clip=5.0',
            '--ifs.master_dark.smoothing=5',
            '--ifs.master_dark.min_acceptable=0.0',
            '--ifs.master_dark.max_acceptable=2000.0',
            '--ifs.master_dark.outfilename={0}'.format(dark_file),
            '--ifs.master_dark.badpixfilename={0}'.format(bpm_file),
            sof]
    proc = subprocess.run(args, stdout=subprocess.DEVNULL)

    if proc.returncode != 0:
        raise ValueError('esorex process was not successful')

    # store products
    calibs_info.set_value(cfiles.index[0], 'products', [dark_file, bpm_file])
    
    ###################################
    # science dark/backgrounds
    ###################################

    # get list of files
    calibs = calibs_info[(calibs_info['type'] == 'background') | (calibs_info['type'] == 'sky')]

    for idx, cal in calibs.iterrows():
        ctype = cal['type']
        DIT = cal['DIT']
        files = cal['files']

        print(' * {0} for science with DIT={1:.2f} sec'.format(ctype, DIT))
        
        # create sof
        sof = os.path.join(sof_path, '{0}_DIT={1:.2f}.sof'.format(ctype, DIT))
        file = open(sof, 'w')
        for f in files:
            file.write('{0}     {1}\n'.format(f, 'IFS_DARK_RAW'))
        file.close()

        # products
        dark_file = '{0}dark_{1}_DIT={2:.2f}.fits'.format(calib_path, ctype, DIT)
        bpm_file  = '{0}bpm_{1}_DIT={2:.2f}.fits'.format(calib_path, ctype, DIT)
    
        # execute esorex    
        args = ['esorex',
                '--msg-level=debug',
                'sph_ifs_master_dark',
                '--ifs.master_dark.coll_alg=2',
                '--ifs.master_dark.sigma_clip=5.0',
                '--ifs.master_dark.smoothing=5',
                '--ifs.master_dark.min_acceptable=0.0',
                '--ifs.master_dark.max_acceptable=2000.0',
                '--ifs.master_dark.outfilename={0}'.format(dark_file),
                '--ifs.master_dark.badpixfilename={0}'.format(bpm_file),
                sof]
        proc = subprocess.run(args, stdout=subprocess.DEVNULL)

        if proc.returncode != 0:
            raise ValueError('esorex process was not successful')

        # store products
        calibs_info.set_value(idx, 'products', [dark_file, bpm_file])
        
    
    
    

root_path = '/Users/avigan/data/pySPHERE-test/IFS/'

# files_info = sort_files(root_path)

files_info = pd.read_csv(root_path+'files.csv', index_col=0)
calibs_info = files_association(root_path, files_info)

calibs_info = pd.read_json(root_path+'calibs.json')
sph_ifs_dark(root_path, calibs_info)
