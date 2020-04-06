'''
VLT/SPHERE filters and neutral densities transmissions
'''

import os
import numpy as np

from scipy.interpolate import interp1d

# definition of IRDIS filter combinations
combinations = {
    # broad-band
    'BB_Y':     {'CFW': 'B_Y',     'DFW': None, 'Wavelength': (1043, 1043), 'Bandwidth': (140, 140)},
    'BB_J':     {'CFW': 'B_J',     'DFW': None, 'Wavelength': (1245, 1245), 'Bandwidth': (240, 240)},
    'BB_H':     {'CFW': 'B_H',     'DFW': None, 'Wavelength': (1625, 1625), 'Bandwidth': (290, 290)},
    'BB_Ks':    {'CFW': 'B_Ks',    'DFW': None, 'Wavelength': (2182, 2182), 'Bandwidth': (300, 300)},

    # dual-band
    'DB_Y23':   {'CFW': 'B_Y',     'DFW': 'D_Y23', 'Wavelength': (1022, 1076), 'Bandwidth': (49, 50)},
    'DB_J23':   {'CFW': 'B_J',     'DFW': 'D_J23', 'Wavelength': (1190, 1273), 'Bandwidth': (42, 46)},
    'DB_H23':   {'CFW': 'B_H',     'DFW': 'D_H23', 'Wavelength': (1593, 1667), 'Bandwidth': (52, 54)},
    'DB_NDH23': {'CFW': 'B_ND-H',  'DFW': 'D_H23', 'Wavelength': (1593, 1667), 'Bandwidth': (52, 54)},
    'DB_H34':   {'CFW': 'B_Hcmp2', 'DFW': 'D_H34', 'Wavelength': (1667, 1733), 'Bandwidth': (54, 57)},
    'DB_K12':   {'CFW': 'B_Ks',    'DFW': 'D_K12', 'Wavelength': (2110, 2251), 'Bandwidth': (102, 109)},

    # narrow-band
    'NB_BrG':   {'CFW': 'N_BrG',   'DFW': None, 'Wavelength': (2170, 2170), 'Bandwidth': (31, 31)},
    'NB_CO':    {'CFW': 'N_CO',    'DFW': None, 'Wavelength': (2290, 2290), 'Bandwidth': (33, 33)},
    'NB_CntH':  {'CFW': 'N_CntH',  'DFW': None, 'Wavelength': (1573, 1573), 'Bandwidth': (23, 23)},
    'NB_CntJ':  {'CFW': 'N_CntJ',  'DFW': None, 'Wavelength': (1213, 1213), 'Bandwidth': (17, 17)},
    'NB_CntK1': {'CFW': 'N_CntK1', 'DFW': None, 'Wavelength': (2091, 2091), 'Bandwidth': (34, 34)},
    'NB_CntK2': {'CFW': 'N_CntK2', 'DFW': None, 'Wavelength': (2266, 2266), 'Bandwidth': (32, 32)},
    'NB_FeII':  {'CFW': 'N_FeII',  'DFW': None, 'Wavelength': (1642, 1642), 'Bandwidth': (24, 24)},
    'NB_H2':    {'CFW': 'N_H2',    'DFW': None, 'Wavelength': (2124, 2124), 'Bandwidth': (31, 31)},
    'NB_HeI':   {'CFW': 'N_HeI',   'DFW': None, 'Wavelength': (1085, 1085), 'Bandwidth': (14, 14)},
    'NB_PaB':   {'CFW': 'N_PaB',   'DFW': None, 'Wavelength': (1283, 1283), 'Bandwidth': (18, 18)},

    # 0-90 deg polarisers
    'DP_0_BB_Y':      {'CFW': 'B_Y',      'DFW': 'P0-90', 'Wavelength': (1043, 1043), 'Bandwidth': (140, 140)},
    'DP_0_BB_J':      {'CFW': 'B_J',      'DFW': 'P0-90', 'Wavelength': (1245, 1245), 'Bandwidth': (240, 240)},
    'DP_0_BB_H':      {'CFW': 'B_H',      'DFW': 'P0-90', 'Wavelength': (1625, 1625), 'Bandwidth': (290, 290)},
    'DP_0_BB_Ks':     {'CFW': 'B_Ks',     'DFW': 'P0-90', 'Wavelength': (2182, 2182), 'Bandwidth': (300, 300)},
    'DP_0_NB_BrG':    {'CFW': 'N_BrG',    'DFW': 'P0-90', 'Wavelength': (2170, 2170), 'Bandwidth': (31, 31)},
    'DP_0_NB_CO':     {'CFW': 'N_CO',     'DFW': 'P0-90', 'Wavelength': (2290, 2290), 'Bandwidth': (33, 33)},
    'DP_0_NB_ContH':  {'CFW': 'N_ContH',  'DFW': 'P0-90', 'Wavelength': (1573, 1573), 'Bandwidth': (23, 23)},
    'DP_0_NB_ContJ':  {'CFW': 'N_ContJ',  'DFW': 'P0-90', 'Wavelength': (1213, 1213), 'Bandwidth': (17, 17)},
    'DP_0_NB_ContK1': {'CFW': 'N_ContK1', 'DFW': 'P0-90', 'Wavelength': (2091, 2091), 'Bandwidth': (34, 34)},
    'DP_0_NB_ContK2': {'CFW': 'N_ContK2', 'DFW': 'P0-90', 'Wavelength': (2266, 2266), 'Bandwidth': (32, 32)},
    'DP_0_NB_FeII':   {'CFW': 'N_FeII',   'DFW': 'P0-90', 'Wavelength': (1642, 1642), 'Bandwidth': (24, 24)},
    'DP_0_NB_H2':     {'CFW': 'N_H2',     'DFW': 'P0-90', 'Wavelength': (2124, 2124), 'Bandwidth': (31, 31)},
    'DP_0_NB_HeI':    {'CFW': 'N_HeI',    'DFW': 'P0-90', 'Wavelength': (1085, 1085), 'Bandwidth': (14, 14)},
    'DP_0_NB_PaB':    {'CFW': 'N_PaB',    'DFW': 'P0-90', 'Wavelength': (1283, 1283), 'Bandwidth': (18, 18)},

    # 45-135 deg polarisers
    'DP_45_BB_Y':      {'CFW': 'B_Y',      'DFW': 'P45-135', 'Wavelength': (1043, 1043), 'Bandwidth': (140, 140)},
    'DP_45_BB_J':      {'CFW': 'B_J',      'DFW': 'P45-135', 'Wavelength': (1245, 1245), 'Bandwidth': (240, 240)},
    'DP_45_BB_H':      {'CFW': 'B_H',      'DFW': 'P45-135', 'Wavelength': (1625, 1625), 'Bandwidth': (290, 290)},
    'DP_45_BB_Ks':     {'CFW': 'B_Ks',     'DFW': 'P45-135', 'Wavelength': (2182, 2182), 'Bandwidth': (300, 300)},
    'DP_45_NB_BrG':    {'CFW': 'N_BrG',    'DFW': 'P45-135', 'Wavelength': (2170, 2170), 'Bandwidth': (31, 31)},
    'DP_45_NB_CO':     {'CFW': 'N_CO',     'DFW': 'P45-135', 'Wavelength': (2290, 2290), 'Bandwidth': (33, 33)},
    'DP_45_NB_ContH':  {'CFW': 'N_ContH',  'DFW': 'P45-135', 'Wavelength': (1573, 1573), 'Bandwidth': (23, 23)},
    'DP_45_NB_ContJ':  {'CFW': 'N_ContJ',  'DFW': 'P45-135', 'Wavelength': (1213, 1213), 'Bandwidth': (17, 17)},
    'DP_45_NB_ContK1': {'CFW': 'N_ContK1', 'DFW': 'P45-135', 'Wavelength': (2091, 2091), 'Bandwidth': (34, 34)},
    'DP_45_NB_ContK2': {'CFW': 'N_ContK2', 'DFW': 'P45-135', 'Wavelength': (2266, 2266), 'Bandwidth': (32, 32)},
    'DP_45_NB_FeII':   {'CFW': 'N_FeII',   'DFW': 'P45-135', 'Wavelength': (1642, 1642), 'Bandwidth': (24, 24)},
    'DP_45_NB_H2':     {'CFW': 'N_H2',     'DFW': 'P45-135', 'Wavelength': (2124, 2124), 'Bandwidth': (31, 31)},
    'DP_45_NB_HeI':    {'CFW': 'N_HeI',    'DFW': 'P45-135', 'Wavelength': (1085, 1085), 'Bandwidth': (14, 14)},
    'DP_45_NB_PaB':    {'CFW': 'N_PaB',    'DFW': 'P45-135', 'Wavelength': (1283, 1283), 'Bandwidth': (18, 18)}
}

# dictionary for already read transmissions
wave_grid = np.arange(900, 2501)
transmissions = {}


def _reinterpolate(tr, wave, new_wave):
    '''
    Reinterpolate a transmission curve on a fixed, regular wavelength grid

    Parameters
    ----------
    tr : array_like
        Transmission, normalized to 1

    wave : array_like
        Wavelengths vector, in nanometers

    new_wave : array_like
        New wavelengths vector, in nanometers

    Returns
    -------
    tr_regular : array_like
        The reinterpolated transmission
    '''

    interp_fun = interp1d(wave, tr, bounds_error=False, fill_value=np.nan)

    return interp_fun(new_wave)


def _load(type, name):
    '''
    Reads the transmission curves from files on disk, reinterpolate
    and save for later use

    Parameters
    ----------
    type : str
        Type of transmission curve to read. Possible values are 'ndf',
        'cfw', 'dfw' and 'ird_ndf'.

    name : str
        Name of the filter.

    Returns
    -------
    value : array_like
        Transmission curve of the filter, reinterpolated on the wave_grid
        reference grid
    '''

    if type == 'ndf':
        # transmission for NDF

        # find file
        package_directory = os.path.dirname(os.path.abspath(__file__))
        filter_file = os.path.join(package_directory, 'data', 'SPHERE_CPI_ND.txt')

        # load data
        ndf_all_tr = np.loadtxt(filter_file, unpack=False).T

        # save for later calls
        transmissions['OPEN']   = _reinterpolate(ndf_all_tr[1], ndf_all_tr[0], wave_grid)
        transmissions['ND_1.0'] = _reinterpolate(ndf_all_tr[2], ndf_all_tr[0], wave_grid)
        transmissions['ND_2.0'] = _reinterpolate(ndf_all_tr[3], ndf_all_tr[0], wave_grid)
        transmissions['ND_3.5'] = _reinterpolate(ndf_all_tr[4], ndf_all_tr[0], wave_grid)

        return transmissions[name]
    elif type == 'cfw':
        # transmission for CFW

        # find file
        package_directory = os.path.dirname(os.path.abspath(__file__))
        filter_file = os.path.join(package_directory, 'data', 'SPHERE_IRDIS_{0}.txt'.format(name))

        # load data
        cfw_tr = np.loadtxt(filter_file, unpack=False).T

        # save for later calls
        transmissions[name] = _reinterpolate(cfw_tr[1], cfw_tr[0], wave_grid)

        return transmissions[name]
    elif type == 'ird_ndf':
        # transmission for IRDIS ND

        # find file
        package_directory = os.path.dirname(os.path.abspath(__file__))
        filter_file = os.path.join(package_directory, 'data', 'SPHERE_IRDIS_ND.txt')

        # load data
        ird_ndf_tr = np.loadtxt(filter_file, unpack=False).T

        # save for later calls
        transmissions['IRD-ND'] = _reinterpolate(ird_ndf_tr[1], ird_ndf_tr[0], wave_grid)

        return transmissions['IRD-ND']
    elif type == 'dfw':
        # transmission for DFW

        # find file
        package_directory = os.path.dirname(os.path.abspath(__file__))
        filter_file = os.path.join(package_directory, 'data', 'SPHERE_IRDIS_{0}.txt'.format(name))

        # load data
        dfw_tr_tmp = np.loadtxt(filter_file, unpack=False).T

        dfw_tr = np.zeros((2, wave_grid.size), dtype=np.float)
        dfw_tr[0] = _reinterpolate(dfw_tr_tmp[1], dfw_tr_tmp[0], wave_grid)
        dfw_tr[1] = _reinterpolate(dfw_tr_tmp[2], dfw_tr_tmp[0], wave_grid)

        # save for later calls
        transmissions[name] = dfw_tr

        return transmissions[name]
    else:
        raise ValueError('Unknown type {0}'.format(type))


def irdis_nd(combination, nd_filter):
    '''
    Provides the IRDIS transmission for a given neutral density and
    filter combination

    Description
    -----------
    This function works for all the IRDIS broad-band (BB), dual-band
    (DB) and narrow-band (NB) filters. The names of the filter combinations
    are provided below.

    Broad-band filters combinations:
     - BB_Y
     - BB_J
     - BB_H
     - BB_Ks

    Dual-band filters combinations:
     - DB_Y23
     - DB_J23
     - DB_H23
     - DB_H32
     - DB_NDH23
     - DB_NDH32
     - DB_H34
     - DB_K12

    Narrow-band filters:
     - NB_BrG
     - NB_CO
     - NB_CntH
     - NB_CntJ
     - NB_CntK1
     - NB_CntK2
     - NB_FeII
     - NB_H2
     - NB_Hel
     - NB_PaB

    The transmission is calculated as the ratio between the
    transmission of the filters with and without multiplication by the
    transmission curve of the CPI neutral density filter. The
    transmission of each filter is read from text files stored on
    disk.

    Parameters
    ----------
    combination : str
        Name of the filter combination. This parameter can be read
        directly from the header of any SPHERE/IRDIS raw file with the
        keyword 'HIERARCH ESO INS COMB IFLT'

    nd_filter : str
        Name of the neutral density filter. This parameter can be read
        directly from the header of any SPHERE/IRDIS raw file with the
        keyword 'HIERARCH ESO INS4 FILT2 NAME'

    Returns
    -------
    tr : array_like
        The transmission of the instrument on the left and right IRDIS
        fields of view

    '''

    # check if combination exists
    setup = combinations.get(combination)
    if setup is None:
        raise ValueError('Unknown filter combination {0}'.format(combination))

    # check if ND filter exists
    if nd_filter not in ['OPEN', 'ND_1.0', 'ND_2.0', 'ND_3.5']:
        raise ValueError('Unknown neutral density filter {0}'.format(nd_filter))

    # setup
    ndf = nd_filter
    cfw = setup['CFW']
    dfw = setup['DFW']

    # transmissions for NDF
    ndf_tr = transmissions.get(ndf)
    if ndf_tr is None:
        ndf_tr = _load('ndf', ndf)

    # transmissions for CFW
    cfw_tr = transmissions.get(cfw)
    if cfw_tr is None:
        cfw_tr = _load('cfw', cfw)

    # transmission for IRDIS ND
    if cfw == 'B_ND-H':
        ird_ndf_tr = transmissions.get('IRD-ND')
        if ird_ndf_tr is None:
            ird_ndf_tr = _load('ird_ndf', None)
    else:
        ird_ndf_tr = np.ones(wave_grid.size)

    # get transmissions for DFW
    if dfw is not None:
        dfw_tr = transmissions.get(dfw)
        if dfw_tr is None:
            dfw_tr = _load('dfw', dfw)
    else:
        # if BB or NB mode, just use 1 for DBF
        dfw_tr = np.ones((2, wave_grid.size))

    # integrated transmission value
    tr_0 = np.nansum(cfw_tr * dfw_tr[0])
    tr_1 = np.nansum(cfw_tr * dfw_tr[1])

    tr_nd_0 = np.nansum(cfw_tr * dfw_tr[0] * ndf_tr * ird_ndf_tr)
    tr_nd_1 = np.nansum(cfw_tr * dfw_tr[1] * ndf_tr * ird_ndf_tr)

    tr = (tr_nd_0 / tr_0, tr_nd_1 / tr_1)

    return tr


def transmission_nd(nd_filter, wave=None):
    '''
    Provides the transmission for a given neutral density

    Description
    -----------

    The function provides the transmission curve of a given CPI
    neutral density filter. The user can provide an array of
    wavelengths.  It works for both IRDIS and IFS, within the range of
    wavalength covered by SPHERE (0.95-2.3 microns).

    Parameters
    ----------
    nd_filter : str
        Name of the neutral density filter. This parameter can be read
        directly from the header of any SPHERE/IRDIS raw file with the
        keyword 'HIERARCH ESO INS4 FILT2 NAME'

    wave : array_like, optional
        Wavelengths at which the transmission is needed, in
        nanometers. Default is None

    Returns
    -------
    wave : array_like
        Wavelength, in nanometers

    tr : array_like
        Transmission of the neutral density filter

    '''

    # check if ND filter exists
    if nd_filter not in ['OPEN', 'ND_1.0', 'ND_2.0', 'ND_3.5']:
        raise ValueError('Unknown neutral density filter {0}'.format(nd_filter))

    ndf = nd_filter

    # get transmissions for NDF
    ndf_tr = transmissions.get(ndf)
    if ndf_tr is None:
        ndf_tr = _load('ndf', ndf)

    if wave is None:
        wave = wave_grid
    else:
        ndf_tr = _reinterpolate(ndf_tr, wave_grid, wave)

    return wave, ndf_tr


def transmission_filter(combination):
    '''
    Provides the IRDIS transmission curve for a given filter combination

    Description
    -----------
    This function works for all the IRDIS broad-band (BB), dual-band
    (DB) and narrow-band (NB) filters. The names of the filter combinations
    are provided below.

    Broad-band filters combinations:
     - BB_Y
     - BB_J
     - BB_H
     - BB_Ks

    Dual-band filters combinations:
     - DB_Y23
     - DB_J23
     - DB_H23
     - DB_H32
     - DB_NDH23
     - DB_NDH32
     - DB_H34
     - DB_K12

    Narrow-band filters:
     - NB_BrG
     - NB_CO
     - NB_CntH
     - NB_CntJ
     - NB_CntK1
     - NB_CntK2
     - NB_FeII
     - NB_H2
     - NB_Hel
     - NB_PaB


    Parameters
    ----------
    combination : str
        Name of the filter combination. This parameter can be read
        directly from the header of any SPHERE/IRDIS raw file with the
        keyword 'HIERARCH ESO INS COMB IFLT'

    Returns
    -------
    wave : array_like
        Wavelength, in nanometers

    tr_0, tr_1 : array_like
        The transmissions of the instrument on the left and right IRDIS
        fields of view

    '''

    # check if combination exists
    setup = combinations.get(combination)
    if setup is None:
        raise ValueError('Unknown filter combination {0}'.format(combination))

    # setup
    cfw = setup['CFW']
    dfw = setup['DFW']

    # transmissions for CFW
    cfw_tr = transmissions.get(cfw)
    if cfw_tr is None:
        cfw_tr = _load('cfw', cfw)

    # transmission for IRDIS ND
    if cfw == 'B_ND-H':
        ird_ndf_tr = transmissions.get('IRD-ND')
        if ird_ndf_tr is None:
            ird_ndf_tr = _load('ird_ndf', None)
    else:
        ird_ndf_tr = np.ones(wave_grid.size)

    # get transmissions for DFW
    if dfw is not None:
        dfw_tr = transmissions.get(dfw)
        if dfw_tr is None:
            dfw_tr = _load('dfw', dfw)
    else:
        # if BB or NB mode, just use 1 for DBF
        dfw_tr = np.ones((2, wave_grid.size))

    # integrated transmission value
    tr_0 = cfw_tr * dfw_tr[0] * ird_ndf_tr
    tr_1 = cfw_tr * dfw_tr[1] * ird_ndf_tr

    return wave_grid, tr_0, tr_1


def wavelength_bandwidth_filter(combination):
    '''
    Provides the wavelength and bandwidth of the filters for a given
    combination

    Description
    -----------
    This function works for all the IRDIS broad-band (BB), dual-band
    (DB) and narrow-band (NB) filters. The names of the filter combinations
    are provided below.

    Broad-band filters combinations:
     - BB_Y
     - BB_J
     - BB_H
     - BB_Ks

    Dual-band filters combinations:
     - DB_Y23
     - DB_J23
     - DB_H23
     - DB_H32
     - DB_NDH23
     - DB_NDH32
     - DB_H34
     - DB_K12

    Narrow-band filters:
     - NB_BrG
     - NB_CO
     - NB_CntH
     - NB_CntJ
     - NB_CntK1
     - NB_CntK2
     - NB_FeII
     - NB_H2
     - NB_Hel
     - NB_PaB


    Parameters
    ----------
    combination : str
        Name of the filter combination. This parameter can be read
        directly from the header of any SPHERE/IRDIS raw file with the
        keyword 'HIERARCH ESO INS COMB IFLT'

    Returns
    -------
    wave : array_like
        Tuple of central wavelengths, in nanometers

    bandwidth : array_like
        Tuple of bandwidth, in nanometers
    '''

    setup = combinations.get(combination)
    if setup is None:
        raise ValueError('Unknown filter combination {0}'.format(combination))

    wave = setup['Wavelength']
    bandwidth = setup['Bandwidth']

    return wave, bandwidth
