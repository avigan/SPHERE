#!/bin/sh

#+
# NAME:
#
#  data_reduction_ifs
#
# PURPOSE:
#
#  Performs the data reduction of SPHERE/IFS data
#
# CALLING SEQUENCE:
#
#  The script is called from the command-line as a normal shell script:
#
#    ~> ./data_reduction_ifs.sh
#
# DESCRIPTION:
#
#  This script provides a basic framework for the data reduction of
#  SPHERE/IFS data. It is based on the official pipeline provided by
#  ESO, plus some additionnal custom routines written in IDL.
#
#  It has been designed to be simple to use, with as little as
#  possible user intervention. For a given observing sequence, there
#  are only two inputs to provide:
#
#    * the full path to the analysis directory ==> ROOT variable
#    * the IFS mode: YJ or YJH                 ==> MODE variable
#
#  The analysis directory MUST contain a sub-directory called "raw",
#  in which you have to put together all the FITS files (science and
#  calibrations). Once all the data is gather in the raw/
#  subdirectory, you have to configure the ROOT and MODE variables at
#  the beginning of the script. The data should *not* be compressed
#  (.Z or .gz extension), and the file name should be in the standard
#  archive form (e.g. SPHER.2015-02-03T06:47:27.430.fits).
#
#  Then the first step is to do a sanity check (DO_SANITY_CHECK=1). It
#  will list the FITS files and will check that all required files (in
#  particular calibrations) are available. For each expected file type
#  it will display the available file names and some warning/error
#  messages if it notices potential problems. Any IFS sequence will
#  require the following files and calibrations:
#
#    * science files of type O, C or F (not all are required)
#    * dark/background with DITs corresponding to the science
#    * sky background with DITs corresponding to the science (optional)
#    * white lamp flat: 2 files
#    * 1020 nm filter flat: 2 files (YJ and YJH modes)
#    * 1230 nm filter flat: 2 files (YJ and YJH modes)
#    * 1300 nm filter flat: 2 files (YJ and YJH modes)
#    * 1550 nm filter flat: 2 files (YJH mode only)
#    * spectra positions: 1 file (YJ and YJH modes)
#    * wavelength calibration: 1 file (YJ and YJH modes)
#    * IFU flat: 1 file (YJ and YJH modes)
#    * dark for the calibration: 1 file, usually taken at minimum DIT (1.6 sec)
#
#  If all these files are available in the raw/ sub-directory, the
#  sanity check should display no errors (but possibly some
#  warnings). If you obtain errors, usually because some files are
#  missing, read the explanation carefully. After the sanity check,
#  the script automatically exits and you have to set
#  DO_SANITY_CHECK=0 to proceed with the data reduction.
#
#  The second step is to create the basic calibrations:
#  dark/backgrounds, detector flats, spectra position, wavelength
#  calibration and IFU flat. To enable the calibration step, you have
#  to set DO_CALIB=1. Then you can select which calibration you will
#  perform by switching the corresponding keywords to 1
#  (e.g. DO_DET_DARK=1) and run the script. If this is your first
#  reduction, I advise to run them one by one. For most of the basic
#  calibrations, the ESO pipeline is used by running the appropriate
#  recipe with esorex. The only exception is the detector flats, which
#  are created with a custom IDL routine. All the calibrations are
#  saved in the calib/ sub-directory. Don't forget to switch all basic
#  calibration keywords to 0 once they are done to avoid recreating
#  them when running the script again. Once all calibrations have been
#  performed, you can set DO_CALIB=0.
#
#  The following step is to run the pre-processing of the raw science
#  files by setting DO_PREPROC=1. The pre-processing is performed with
#  the custom IDL routine sph_ifs_preprocess(). There are several
#  options for this pre-processing, which are triggered by the
#  DO_COLLAPSE, DO_BKG_SUB, DO_BADPIXEL and DO_CROSSTALK keywords. The
#  effect of these options are documented in the header of the
#  sph_ifs_preprocess() routine. Brief summary:
#
#    * DO_COLLAPSE: performs a temporal binning of the data
#    * DO_BKG_SUB: performs the background subtraction
#    * DO_BADPIXEL: performs the bad pixels correction
#    * DO_CROSSTALK: performs the spectral crosstalk correction
#
#  The pre-processed raw science files are then saved in the interm/
#  subdirectory along with "sidecar" FITS files that contain some
#  useful information for each of the science frames (time, hour
#  angle, position angle, ...). For a standard reduction when you want
#  to keep all frames, disable the collapse, and enable all the others
#  (background subtraction, bad pixel correction and crosstalk
#  correction).
#
#  For the pre-processing step, there is the additional USE_SKY
#  keyword. This keyword forces the use of sky backgrounds if both
#  instrumental dark/backgrounds and sky backgrounds are
#  available. This usually provides slightly better cosmetics and
#  should be left enabled by default.
#
#  Once pre-processing is done, you can set DO_PREPROC=0 before going
#  to the next step.
#
#  The last step is to produce the (x,y,lambda) science cubes using as
#  input the pre-processed raw files. This is done by setting
#  DO_SCI_CORO=1. In this step, the ESO pipeline is run to generate
#  the cubes, which are copied into the products/ sub-directory at the
#  end. The keyword DO_POSTPROC=1 enables a simple post-processing
#  step of the science cubes that removes some unnecessary FITS
#  extensions that make the files much smaller (~16 MB instead of ~200
#  MB).
#
#  The wavelength calibration performed by the ESO pipeline is not
#  optimal. To be able to recalibrate properly the wavelength, the
#  script provides an additionnal step that is triggered by the
#  DO_SCI_LMBD keyword. In this procedure, the wavelength calibration
#  file is processed as a science frame and then used in the analysis
#  in combination with a star center frame to perform a clean
#  recalibration of the wavelength of each IFS spectral channels. The
#  complete procedure is described in Vigan et al. 2015 (MNRAS, 454,
#  129) and is included in the small data_reduction_ifs.pro IDL
#  program that comes with this script.
#
#  After the execution of all the steps from the data reduction
#  script, the products/ subdirectory should contain clean
#  (x,y,lambda) data cubes that can be used for science. They will
#  most likely still contain a small number of bad pixels, which can
#  easily be cleaned with a sigma-clipping procedure. Note that the
#  centring of the science cubes is not (yet) handled at the level of
#  the ESO pipeline, so additional code is required to perform an
#  accurate centring using star center (C) frames. I provide an
#  additionnal IDL code data_reduction_ifs.pro that does these
#  additional few steps, including the wavelength recalibration.
#
#  The steps required to obtain clean SPHERE/IFS data cubes are not
#  all straighforward. I hope that this script will simplify some of
#  these steps and will allow to quickly obtain results. If you
#  encounter any problem using this tool, please do not hesitate to
#  contact me.
#
# REQUIRED PROGRAMS AND LIBRARIES:
#
#  To be able to work properly, the reduction pipeline relies on a
#  small number of programs and IDL libraries that are listed below.
#
#  Command line tools:
#
#    * SPHERE pipeline: official pipeline available on the ESO website
#         http://www.eso.org/sci/software/pipelines/
#         Make sure that the esorex tool is working properly before
#         attempting any IFS reduction using this script
#
#    * dfits and fitsort: essential command line tools to display FITS
#         headers and values of specific keywords. They are part of
#         the ESO eclipse library available on the ESO website
#         http://www.eso.org/sci/software/eclipse/distrib/index.html
#
#  IDL libraries:
#
#    * astronomy user's library: http://idlastro.gsfc.nasa.gov/
#    * maskinterp: http://astro.vigan.fr/tools/maskinterp-1.2.tar.gz
#
# MODIFICATION HISTORY:
#
#  arthur.vigan - 08/2016 - v1.2
#                 small improvements. Added support for latest version
#                 of the ESO pipeline
#
#  arthur.vigan - 08/2015 - v1.1
#                 added sanity check, commented for distribution
#
#  arthur.vigan - 07/2015
#                 added keyword for waffle observations, added background
#                 specific for calibrations to improve cosmetics
#
#  arthur.vigan - 05/2015
#                 updated for latest version of ESO pipeline
#
#  arthur.vigan - 11/2014
#                 several small improvements including automatic
#                 classification of files
#
#  arthur.vigan - 10/2014
#                first usable version with main features
#
# AUTHOR:
# 
#   Arthur Vigan
#   Laboratoire d'Astrophysique de Marseille
#   arthur.vigan@lam.fr
#
# LICENSE:
#
#   This code is release under the MIT license. The full text of the
#   license is included in a separate file LICENSE.txt.
#
#   The developement of the SPHERE instrument has demanded a
#   tremendous effort from many scientists, who have devoted several
#   years of their life to design, build, test and commission this new
#   instrument. To recognize this work, we kindly ask you to cite the
#   relevant papers in your scientific work. More specifically,
#   because this script is the core of our public SPHERE/IFS reduction
#   pipeline, we would be grateful if you could cite the following
#   paper in any publication making use of it:
#
#     Vigan et al., 2015, MNRAS, 454, 129
#
#   We thank you for your effort, and hope that this tool will
#   contribute to your scientific work and discoveries. Please feel
#   free to report any bug or possible improvement to the author
#
#-

ROOT=/data/TargetDirectory/
MODE=YJ

#
# analysis parameters
#

# first step: sanity check!
DO_SANITY_CHECK=1       # sanity checks

# basic calibrations
DO_CALIB=0              # perform basic calibrations

DO_DET_DARK=1           # dark/background or sky
DO_DET_FLAT=1           # detector flats
DO_SPEC_POS=1           # IFU spectra position
DO_WAVE_CAL=1           # wavelength calibration
DO_IFU_FLAT=1           # IFU lenslet flat

# pre-processing
DO_PREPROC=0            # perform pre-processing

DO_COLLAPSE=1           # collapse cubes
DO_BKG_SUB=1            # subtract dark/background or sky
DO_BADPIXEL=1           # correct bad pixels
DO_CROSSTALK=1          # correct spectral cross-talk

COLLAPSE_TYPE=coadd     # collapse type: mean, angle, coadd
COLLAPSE_VAL=5          # collapse parameter (for mean, angle and coadd)
COLLAPSE_TOL=0.05       # collapse parameter tolerance (for mean and angle)

USE_SKY=1               # if available, use preferably a sky instead of a dark/background

# science data reduction
DO_SCIENCE=0            # perform science

DO_SCI_CORO=1           # process science data
DO_SCI_LMBD=1           # process wave cal cube

DO_POSTPROC=1           # perform post-processing on (x,y,lambda) cubes

# additional suffix
SUFFIX=
SUFFIX_PSF=

####################################################################

# create products directory if necessary
if [ ! -d ${ROOT}sof      ]; then mkdir ${ROOT}sof;      fi
if [ ! -d ${ROOT}calib    ]; then mkdir ${ROOT}calib;    fi
if [ ! -d ${ROOT}interm   ]; then mkdir ${ROOT}interm;   fi
if [ ! -d ${ROOT}products ]; then mkdir ${ROOT}products; fi

# build suffix
if [ $DO_COLLAPSE -ne 0 ]; then
    SUFFIX=${SUFFIX}_col_${COLLAPSE_TYPE}
fi
SUFFIX_PSF=${SUFFIX_PSF}_col_mean

if [ $DO_BKG_SUB -ne 0 ]; then
    SUFFIX=${SUFFIX}_bkg
    SUFFIX_PSF=${SUFFIX_PSF}_bkg
fi

if [ $DO_BADPIXEL -ne 0 ]; then
    SUFFIX=${SUFFIX}_bp
    SUFFIX_PSF=${SUFFIX_PSF}_bp
fi

if [ $DO_CROSSTALK -ne 0 ]; then
    SUFFIX=${SUFFIX}_ct
    SUFFIX_PSF=${SUFFIX_PSF}_ct
fi

#
# SORT FILES
#
echo 'List of files in raw directory:'
cd ${ROOT}raw/
dfits *.fits | fitsort dpr.type ins2.comb.ifs det.seq1.dit det.ndit tel.parang.start
cd ${ROOT}
echo
echo

#
# flats
#
echo '* White flat files:'
flat_whit_files=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit ins2.comb.ifs | grep FLAT | grep CAL_BB_2_${MODE}  | awk '{print $1}'`)
for (( i = 0 ; i < ${#flat_whit_files[@]} ; i++ )); do flat_whit_files[$i]=`basename ${flat_whit_files[$i]} .fits`; echo "  + ${flat_whit_files[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#flat_whit_files[@]} -ne 2 ]; then
        echo '  Error: there should be 2 flat files for white lamp!'
    fi
fi
echo

echo '* 1020 nm flat files:'
flat_1020_files=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit ins2.comb.ifs | grep FLAT | grep CAL_NB1_1_${MODE} | awk '{print $1}'`)
for (( i = 0 ; i < ${#flat_1020_files[@]} ; i++ )); do flat_1020_files[$i]=`basename ${flat_1020_files[$i]} .fits`; echo "  + ${flat_1020_files[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#flat_1020_files[@]} -ne 2 ]; then
        echo '  Error: there should be 2 flat files for 1020 nm filter!'
    fi
fi
echo

echo '* 1230 nm flat files:'
flat_1230_files=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit ins2.comb.ifs | grep FLAT | grep CAL_NB2_1_${MODE} | awk '{print $1}'`)
for (( i = 0 ; i < ${#flat_1230_files[@]} ; i++ )); do flat_1230_files[$i]=`basename ${flat_1230_files[$i]} .fits`; echo "  + ${flat_1230_files[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#flat_1230_files[@]} -ne 2 ]; then
        echo '  Error: there should be 2 flat files for 1230 nm filter!'
    fi
fi
echo

echo '* 1300 nm flat files:'
flat_1300_files=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit ins2.comb.ifs | grep FLAT | grep CAL_NB3_1_${MODE} | awk '{print $1}'`)
for (( i = 0 ; i < ${#flat_1300_files[@]} ; i++ )); do flat_1300_files[$i]=`basename ${flat_1300_files[$i]} .fits`; echo "  + ${flat_1300_files[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#flat_1300_files[@]} -ne 2 ]; then
        echo '  Error: there should be 2 flat files for 1300 nm filter!'
    fi
fi
echo

echo '* 1550 nm flat files:'
flat_1550_files=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit ins2.comb.ifs | grep FLAT | grep CAL_NB4_2_${MODE} | awk '{print $1}'`)
for (( i = 0 ; i < ${#flat_1550_files[@]} ; i++ )); do flat_1550_files[$i]=`basename ${flat_1550_files[$i]} .fits`; echo "  + ${flat_1550_files[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#flat_1550_files[@]} -ne 2 ] && [ ${MODE} = 'YJH' ]; then
        echo '  Error: there should be 2 flat files for 1550 nm filter in YJH mode!'
    fi
fi
echo

#
# spec pos
#    
if [ ${MODE} = 'YJ' ]; then OBS=OBS_YJ; else OBS=OBS_H; fi
echo '* Spectra position files:'
spec_files=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit ins2.comb.ifs | grep SPECPOS,LAMP | grep ${OBS} | awk '{print $1}'`)    
for (( i = 0 ; i < ${#spec_files[@]} ; i++ )); do spec_files[$i]=`basename ${spec_files[$i]} .fits`; echo "  + ${spec_files[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#spec_files[@]} -ne 1 ]; then
        echo '  Error: there should be 1 spectra position file!'
    fi
fi
echo

#
# wavelength
#
echo '* Wavelength calibration files:'
wave_files=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit ins2.comb.ifs | grep WAVE,LAMP | grep ${OBS} | awk '{print $1}'`)
for (( i = 0 ; i < ${#wave_files[@]} ; i++ )); do wave_files[$i]=`basename ${wave_files[$i]} .fits`; echo "  + ${wave_files[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#wave_files[@]} -ne 1 ]; then
        echo '  Error: there should be 1 wavelength calibration file!'
    fi
fi
echo

#
# IFU flat
#
echo '* IFU flat files:'
ifu_files=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit ins2.comb.ifs | grep FLAT,LAMP | grep ${OBS} | awk '{print $1}'`)
for (( i = 0 ; i < ${#ifu_files[@]} ; i++ )); do ifu_files[$i]=`basename ${ifu_files[$i]} .fits`; echo "  + ${ifu_files[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#ifu_files[@]} -ne 1 ]; then
        echo '  Error: there should be 1 IFU flat file!'
    fi
fi
echo

#
# PSF
#
echo '* PSF files:'
sci_files_psf=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit ins2.comb.ifs | grep OBJECT,FLUX | awk '{print $1}'`)
for (( i = 0 ; i < ${#sci_files_psf[@]} ; i++ )); do sci_files_psf[$i]=`basename ${sci_files_psf[$i]} .fits`; echo "  + ${sci_files_psf[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#sci_files_psf[@]} -eq 0 ]; then
        echo '  Warning: there are no PSF (F) files'
    fi
fi
echo

#
# star center
#
echo '* Star center files:'
sci_files_cen=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit ins2.comb.ifs | grep OBJECT | grep OBJECT,CENTER | awk '{print $1}'`)
for (( i = 0 ; i < ${#sci_files_cen[@]} ; i++ )); do sci_files_cen[$i]=`basename ${sci_files_cen[$i]} .fits`; echo "  + ${sci_files_cen[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#sci_files_cen[@]} -eq 0 ]; then
        echo '  Warning: there are no star center (C) files. You will not be able to perform a recalibration of the wavelength'
    fi
fi
echo

#
# corono
#
echo '* Corono data files:'
sci_files_cor=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit ins2.comb.ifs | grep OBJECT | grep -v OBJECT,FLUX | grep -v OBJECT,CENTER | awk '{print $1}'`)
for (( i = 0 ; i < ${#sci_files_cor[@]} ; i++ )); do sci_files_cor[$i]=`basename ${sci_files_cor[$i]} .fits`; echo "  + ${sci_files_cor[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#sci_files_cor[@]} -eq 0 ]; then
        echo '  Warning: there are no science (O) files. This is OK if your sequence is a continuous waffle observation.'
    fi
fi
echo
    
#
# darks
#    
echo '* Calib data dark files:'
DIT=1.6
dark_files_cal=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit | grep DARK | grep "${DIT}" | awk '{print $1}'`)
for (( i = 0 ; i < ${#dark_files_cal[@]} ; i++ )); do dark_files_cal[$i]=`basename ${dark_files_cal[$i]} .fits`; echo "  + ${dark_files_cal[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#dark_files_cal[@]} -eq 0 ]; then
        echo '  Warning: there is no dark/background for the basic calibrations (DIT=1.6 sec). It is *highly recommended* to include one to obtain the best data reduction. A single dark/background file is sufficient, and it can easily be downloaded from the ESO archive'
    fi
fi
echo

echo '* PSF data dark files:'
DIT=`dfits ${ROOT}raw/${sci_files_psf[0]}.fits | fitsort det.seq1.dit | grep -v FILE | awk '{print $2}'`
dark_files_psf=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit | grep DARK | grep "${DIT}" | awk '{print $1}'`)
for (( i = 0 ; i < ${#dark_files_psf[@]} ; i++ )); do dark_files_psf[$i]=`basename ${dark_files_psf[$i]} .fits`; echo "  + ${dark_files_psf[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#dark_files_psf[@]} -eq 0 ]; then
        echo '  Warning: there is no dark/background for the PSF (F) files. If you have PSF (F) files, it is *highly recommended* to include one to obtain the best data reduction. A single dark/background file is sufficient, and it can easily be downloaded from the ESO archive'
    fi
fi
echo

echo '* Star center dark files:'
DIT=`dfits ${ROOT}raw/${sci_files_cen[0]}.fits | fitsort det.seq1.dit | grep -v FILE | awk '{print $2}'`
dark_files_cen=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit | grep DARK | grep "${DIT}" | awk '{print $1}'`)
for (( i = 0 ; i < ${#dark_files_cen[@]} ; i++ )); do dark_files_cen[$i]=`basename ${dark_files_cen[$i]} .fits`; echo "  + ${dark_files_cen[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#dark_files_cen[@]} -eq 0 ]; then
        echo '  Warning: there is no dark/background for the star center (C) files. If you have star center (C) files, it is *highly recommended* to include one to obtain the best data reduction. A single dark/background file is sufficient, and it can easily be downloaded from the ESO archive'
    fi
fi
echo

echo '* Corono data dark files:'
DIT=`dfits ${ROOT}raw/${sci_files_cor[0]}.fits | fitsort det.seq1.dit | grep -v FILE | awk '{print $2}'`
dark_files_cor=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit | grep DARK | grep "${DIT}" | awk '{print $1}'`)
for (( i = 0 ; i < ${#dark_files_cor[@]} ; i++ )); do dark_files_cor[$i]=`basename ${dark_files_cor[$i]} .fits`; echo "  + ${dark_files_cor[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#dark_files_cor[@]} -eq 0 ]; then
        echo '  Warning: there is no dark/background for the science (O) files. If you have science (O) files, it is *highly recommended* to include one to obtain the best data reduction. A single dark/background file is sufficient, and it can easily be downloaded from the ESO archive'
    fi
fi
echo

#
# sky    
#
echo '* PSF data sky files:'
DIT=`dfits ${ROOT}raw/${sci_files_psf[0]}.fits | fitsort det.seq1.dit | grep -v FILE | awk '{print $2}'`
sky_files_psf=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit | grep SKY | grep "${DIT}" | awk '{print $1}'`)
for (( i = 0 ; i < ${#sky_files_psf[@]} ; i++ )); do sky_files_psf[$i]=`basename ${sky_files_psf[$i]} .fits`; echo "  + ${sky_files_psf[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#sky_files_psf[@]} -eq 0 ]; then
        echo '  Warning: there is no sky background for the PSF (F) files. If you have PSF (F) files, using a sky background instead of an internal instrumental background can sometimes provide a cleaner data reduction'
    fi
fi
echo

echo '* Star center sky files:'
DIT=`dfits ${ROOT}raw/${sci_files_cen[0]}.fits | fitsort det.seq1.dit | grep -v FILE | awk '{print $2}'`
sky_files_cen=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit | grep SKY | grep "${DIT}" | awk '{print $1}'`)
for (( i = 0 ; i < ${#sky_files_cen[@]} ; i++ )); do sky_files_cen[$i]=`basename ${sky_files_cen[$i]} .fits`; echo "  + ${sky_files_cen[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#sky_files_cen[@]} -eq 0 ]; then
        echo '  Warning: there is no sky background for the star center (C) files. If you have star center (C) files, using a sky background instead of an internal instrumental background can sometimes provide a cleaner data reduction'
    fi
fi
echo

echo '* Corono data sky files:'
DIT=`dfits ${ROOT}raw/${sci_files_cor[0]}.fits | fitsort det.seq1.dit | grep -v FILE | awk '{print $2}'`
sky_files_cor=(`dfits ${ROOT}raw/SPHER*[0-9].fits | fitsort dpr.type det.seq1.dit | grep SKY | grep "${DIT}" | awk '{print $1}'`)
for (( i = 0 ; i < ${#sky_files_cor[@]} ; i++ )); do sky_files_cor[$i]=`basename ${sky_files_cor[$i]} .fits`; echo "  + ${sky_files_cor[$i]}"; done
if [ $DO_SANITY_CHECK -ne 0 ]; then
    if [ ${#sky_files_cor[@]} -eq 0 ]; then
        echo '  Warning: there is no sky background for the science (O) files. If you have science (O) files, using a sky background instead of an internal instrumental background can sometimes provide a cleaner data reduction'
    fi
fi
echo

# exit if we are doing the sanity check
if [ $DO_SANITY_CHECK -ne 0 ]; then exit; fi

export ESO_BASE=ESO
export DYLD_LIBRARY_PATH=$HOME/${ESO_BASE}/lib:$DYLD_LIBRARY_PATH

# Waffle-only observation?
if [ ${#sci_files_cor[@]} -eq 0 ] && [ ${#sci_files_cen[@]} -ne 0 ]; then
    WAFFLE_OBS=1
    echo 'Information: this sequence is a continuous waffle observation'
else
    WAFFLE_OBS=0
fi 
echo

#
# basic calibrations
#
if [ $DO_CALIB -ne 0 ]; then
    #
    # DARK + SKY
    #
    if [ $DO_DET_DARK -ne 0 ]; then
        #
        # calibs dark
        #
        if [ -n "$dark_files_cal" ]; then
    	    SOF=${ROOT}sof/dark.sof

    	    if [ -f ${SOF} ]; then rm ${SOF}; fi
    	    for f in ${dark_files_cal[*]}
    	    do
    	        echo "${ROOT}raw/${f}.fits   IFS_DARK_RAW" >> ${SOF}
    	    done

    	    esorex --msg-level=debug sph_ifs_master_dark \
    	           --ifs.master_dark.coll_alg=2 \
    	           --ifs.master_dark.sigma_clip=5.0 \
    	           --ifs.master_dark.smoothing=5 \
    	           --ifs.master_dark.min_acceptable=0.0 \
    	           --ifs.master_dark.max_acceptable=2000.0 \
    	           --ifs.master_dark.outfilename=${ROOT}calib/dark_cal.fits \
    	           --ifs.master_dark.badpixfilename=${ROOT}calib/dark_bpm_cal.fits \
    	           ${SOF}
        fi
        
        #
        # PSF dark
        #
        if [ -n "$dark_files_psf" ]; then
    	    SOF=${ROOT}sof/dark.sof

    	    if [ -f ${SOF} ]; then rm ${SOF}; fi
    	    for f in ${dark_files_psf[*]}
    	    do
    	        echo "${ROOT}raw/${f}.fits   IFS_DARK_RAW" >> ${SOF}
    	    done
	    
    	    esorex --msg-level=debug sph_ifs_master_dark \
    	           --ifs.master_dark.coll_alg=2 \
    	           --ifs.master_dark.sigma_clip=3.0 \
    	           --ifs.master_dark.smoothing=5 \
    	           --ifs.master_dark.min_acceptable=0.0 \
    	           --ifs.master_dark.max_acceptable=2000.0 \
    	           --ifs.master_dark.outfilename=${ROOT}calib/dark_psf.fits \
    	           --ifs.master_dark.badpixfilename=${ROOT}calib/dark_bpm_psf.fits \
    	           ${SOF}
        fi

        #
        # star center dark
        #
        if [ -n "$dark_files_cen" ]; then
    	    SOF=${ROOT}sof/dark.sof
	    
    	    if [ -f ${SOF} ]; then rm ${SOF}; fi
    	    for f in ${dark_files_cen[*]}
    	    do
    	        echo "${ROOT}raw/${f}.fits   IFS_DARK_RAW" >> ${SOF}
    	    done
	    
    	    esorex --msg-level=debug sph_ifs_master_dark \
    	           --ifs.master_dark.coll_alg=2 \
    	           --ifs.master_dark.sigma_clip=3.0 \
    	           --ifs.master_dark.smoothing=5 \
    	           --ifs.master_dark.min_acceptable=0.0 \
    	           --ifs.master_dark.max_acceptable=2000.0 \
    	           --ifs.master_dark.outfilename=${ROOT}calib/dark_cen.fits \
    	           --ifs.master_dark.badpixfilename=${ROOT}calib/dark_bpm_cen.fits \
    	           ${SOF}
        fi

        #
        # corono dark
        #
        if [ -n "$dark_files_cor" ]; then
    	    SOF=${ROOT}sof/dark.sof
	    
    	    if [ -f ${SOF} ]; then rm ${SOF}; fi
    	    for f in ${dark_files_cor[*]}
    	    do
    	        echo "${ROOT}raw/${f}.fits   IFS_DARK_RAW" >> ${SOF}
    	    done
	    
    	    esorex --msg-level=debug sph_ifs_master_dark \
    	           --ifs.master_dark.coll_alg=2 \
    	           --ifs.master_dark.sigma_clip=3.0 \
    	           --ifs.master_dark.smoothing=5 \
    	           --ifs.master_dark.min_acceptable=0.0 \
    	           --ifs.master_dark.max_acceptable=2000.0 \
    	           --ifs.master_dark.outfilename=${ROOT}calib/dark_cor.fits \
    	           --ifs.master_dark.badpixfilename=${ROOT}calib/dark_bpm_cor.fits \
    	           ${SOF}
        fi
        
        #
        # PSF sky
        #
        if [ -n "$sky_files_psf" ]; then
    	    SOF=${ROOT}sof/sky.sof
	    
    	    if [ -f ${SOF} ]; then rm ${SOF}; fi
    	    for f in ${sky_files_psf[*]}
    	    do
    	        echo "${ROOT}raw/${f}.fits   IFS_DARK_RAW" >> ${SOF}
    	    done
	    
    	    esorex --msg-level=debug sph_ifs_master_dark \
    	           --ifs.master_dark.coll_alg=2 \
    	           --ifs.master_dark.sigma_clip=3.0 \
    	           --ifs.master_dark.smoothing=5 \
    	           --ifs.master_dark.min_acceptable=0.0 \
    	           --ifs.master_dark.max_acceptable=2000.0 \
    	           --ifs.master_dark.outfilename=${ROOT}calib/sky_psf.fits \
    	           --ifs.master_dark.badpixfilename=${ROOT}calib/sky_bpm_psf.fits \
    	           ${SOF}
        fi
        
        #
        # star center sky
        #
        if [ -n "$sky_files_cen" ]; then
    	    SOF=${ROOT}sof/sky.sof
	    
    	    if [ -f ${SOF} ]; then rm ${SOF}; fi
    	    for f in ${sky_files_cen[*]}
    	    do
    	        echo "${ROOT}raw/${f}.fits   IFS_DARK_RAW" >> ${SOF}
    	    done
	    
    	    esorex --msg-level=debug sph_ifs_master_dark \
    	           --ifs.master_dark.coll_alg=2 \
    	           --ifs.master_dark.sigma_clip=3.0 \
    	           --ifs.master_dark.smoothing=5 \
    	           --ifs.master_dark.min_acceptable=0.0 \
    	           --ifs.master_dark.max_acceptable=2000.0 \
    	           --ifs.master_dark.outfilename=${ROOT}calib/sky_cen.fits \
    	           --ifs.master_dark.badpixfilename=${ROOT}calib/sky_bpm_cen.fits \
    	           ${SOF}
        fi
        
        #
        # corono sky
        #
        if [ -n "$sky_files_cor" ]; then
	    SOF=${ROOT}sof/sky.sof
	    
	    if [ -f ${SOF} ]; then rm ${SOF}; fi
	    for f in ${sky_files_cor[*]}
	    do
	        echo "${ROOT}raw/${f}.fits   IFS_DARK_RAW" >> ${SOF}
	    done
	    
	    esorex --msg-level=debug sph_ifs_master_dark \
	           --ifs.master_dark.sigma_clip=3.0 \
	           --ifs.master_dark.min_acceptable=0.0 \
	           --ifs.master_dark.max_acceptable=2000.0 \
	           --ifs.master_dark.outfilename=${ROOT}calib/sky_cor.fits \
	           --ifs.master_dark.badpixfilename=${ROOT}calib/sky_bpm_cor.fits \
	           ${SOF}
        fi
    fi

    #
    # FLAT
    #
    if [ $DO_DET_FLAT -ne 0 ]; then
        #
        # white flat
        #
        SOF=${ROOT}sof/flat_white.sof

        if [ -f ${SOF} ]; then rm ${SOF}; fi
        for f in ${flat_whit_files[*]}
        do
	    echo "${ROOT}raw/${f}.fits         IFS_DETECTOR_FLAT_FIELD_RAW" >> ${SOF}
        done
        echo "${ROOT}calib/dark_bpm_cal.fits   IFS_STATIC_BADPIXELMAP"      >> ${SOF}
        echo "${ROOT}calib/dark_bpm_cen.fits   IFS_STATIC_BADPIXELMAP"      >> ${SOF}
        echo "${ROOT}calib/dark_bpm_cor.fits   IFS_STATIC_BADPIXELMAP"      >> ${SOF}
        
        # esorex --msg-level=debug sph_ifs_master_detector_flat \
        # 	   --ifs.master_detector_flat.save_addprod=TRUE \
        # 	   --ifs.master_detector_flat.outfilename=${ROOT}calib/master_detector_flat_white_drh.fits \
        # 	   --ifs.master_detector_flat.lss_outfilename=${ROOT}calib/large_scale_flat_white_drh.fits \
        # 	   --ifs.master_detector_flat.preamp_outfilename=${ROOT}calib/preamp_flat_white_drh.fits \
        # 	   --ifs.master_detector_flat.badpixfilename=${ROOT}calib/dff_badpixelname_white_drh.fits \
        # 	   --ifs.master_detector_flat.lambda=-1.0 \
        # 	   --ifs.master_detector_flat.smoothing_length=10.0 \
        # 	   --ifs.master_detector_flat.smoothing_method=1 \
        # 	   ${SOF}
        
        FFNAME=${ROOT}calib/master_detector_flat_white.fits
        BPNAME=${ROOT}calib/dff_badpixelname_white.fits
        idl -e sph_ifs_detector_flat_manual -args ${SOF} ${FFNAME} ${BPNAME}
        
        if [ $? -ne 0 ]; then
    	    echo "IDL returned an error value. Stopping here..."
    	    exit 1
        fi
        
        #
        # 1020 nm flat
        #
        SOF=${ROOT}sof/flat_1020.sof

        if [ -f ${SOF} ]; then rm ${SOF}; fi
        for f in ${flat_1020_files[*]}
        do
	    echo "${ROOT}raw/${f}.fits         IFS_DETECTOR_FLAT_FIELD_RAW" >> ${SOF}
        done
        echo "${ROOT}calib/dark_bpm_cal.fits   IFS_STATIC_BADPIXELMAP"      >> ${SOF}
        echo "${ROOT}calib/dark_bpm_cen.fits   IFS_STATIC_BADPIXELMAP"      >> ${SOF}
        echo "${ROOT}calib/dark_bpm_cor.fits   IFS_STATIC_BADPIXELMAP"      >> ${SOF}
        
        # esorex --msg-level=debug sph_ifs_master_detector_flat \
        # 	   --ifs.master_detector_flat.save_addprod=TRUE \
        # 	   --ifs.master_detector_flat.outfilename=${ROOT}calib/master_detector_flat_1020.fits \
        # 	   --ifs.master_detector_flat.lss_outfilename=${ROOT}calib/large_scale_flat_1020.fits \
        # 	   --ifs.master_detector_flat.preamp_outfilename=${ROOT}calib/preamp_flat_1020.fits \
        # 	   --ifs.master_detector_flat.badpixfilename=${ROOT}calib/dff_badpixelname_1020.fits \
        # 	   --ifs.master_detector_flat.lambda=1.020 \
        # 	   --ifs.master_detector_flat.smoothing_length=10.0 \
        # 	   --ifs.master_detector_flat.smoothing_method=1 \
        # 	   ${SOF}

        FFNAME=${ROOT}calib/master_detector_flat_1020.fits
        BPNAME=${ROOT}calib/dff_badpixelname_1020.fits
        idl -e sph_ifs_detector_flat_manual -args ${SOF} ${FFNAME} ${BPNAME}
        if [ $? -ne 0 ]; then
    	    echo "IDL returned an error value. Stopping here..."
    	    exit 1
        fi
        
        #
        # 1230 nm flat
        #
        SOF=${ROOT}sof/flat_1230.sof

        if [ -f ${SOF} ]; then rm ${SOF}; fi
        for f in ${flat_1230_files[*]}
        do
	    echo "${ROOT}raw/${f}.fits         IFS_DETECTOR_FLAT_FIELD_RAW" >> ${SOF}
        done
        echo "${ROOT}calib/dark_bpm_cal.fits   IFS_STATIC_BADPIXELMAP"      >> ${SOF}
        echo "${ROOT}calib/dark_bpm_cen.fits   IFS_STATIC_BADPIXELMAP"      >> ${SOF}
        echo "${ROOT}calib/dark_bpm_cor.fits   IFS_STATIC_BADPIXELMAP"      >> ${SOF}
        
        # esorex --msg-level=debug sph_ifs_master_detector_flat \
        # 	   --ifs.master_detector_flat.save_addprod=TRUE \
        # 	   --ifs.master_detector_flat.outfilename=${ROOT}calib/master_detector_flat_1230.fits \
        # 	   --ifs.master_detector_flat.lss_outfilename=${ROOT}calib/large_scale_flat_1230.fits \
        # 	   --ifs.master_detector_flat.preamp_outfilename=${ROOT}calib/preamp_flat_1230.fits \
        # 	   --ifs.master_detector_flat.badpixfilename=${ROOT}calib/dff_badpixelname_1230.fits \
        # 	   --ifs.master_detector_flat.lambda=1.230 \
        # 	   --ifs.master_detector_flat.smoothing_length=10.0 \
        # 	   --ifs.master_detector_flat.smoothing_method=1 \
        # 	   ${SOF}

        FFNAME=${ROOT}calib/master_detector_flat_1230.fits
        BPNAME=${ROOT}calib/dff_badpixelname_1230.fits
        idl -e sph_ifs_detector_flat_manual -args ${SOF} ${FFNAME} ${BPNAME}
        if [ $? -ne 0 ]; then
    	    echo "IDL returned an error value. Stopping here..."
    	    exit 1
        fi
        
        #
        # 1300 nm flat
        #
        SOF=${ROOT}sof/flat_1300.sof

        if [ -f ${SOF} ]; then rm ${SOF}; fi
        for f in ${flat_1300_files[*]}
        do
	    echo "${ROOT}raw/${f}.fits         IFS_DETECTOR_FLAT_FIELD_RAW" >> ${SOF}
        done
        echo "${ROOT}calib/dark_bpm_cal.fits   IFS_STATIC_BADPIXELMAP"      >> ${SOF}
        echo "${ROOT}calib/dark_bpm_cen.fits   IFS_STATIC_BADPIXELMAP"      >> ${SOF}
        echo "${ROOT}calib/dark_bpm_cor.fits   IFS_STATIC_BADPIXELMAP"      >> ${SOF}
        
        # esorex --msg-level=debug sph_ifs_master_detector_flat \
        # 	   --ifs.master_detector_flat.save_addprod=TRUE \
        # 	   --ifs.master_detector_flat.outfilename=${ROOT}calib/master_detector_flat_1300.fits \
        # 	   --ifs.master_detector_flat.lss_outfilename=${ROOT}calib/large_scale_flat_1300.fits \
        # 	   --ifs.master_detector_flat.preamp_outfilename=${ROOT}calib/preamp_flat_1300.fits \
        # 	   --ifs.master_detector_flat.badpixfilename=${ROOT}calib/dff_badpixelname_1300.fits \
        # 	   --ifs.master_detector_flat.lambda=1.300 \
        # 	   --ifs.master_detector_flat.smoothing_length=10.0 \
        # 	   --ifs.master_detector_flat.smoothing_method=1 \
        # 	   ${SOF}

        FFNAME=${ROOT}calib/master_detector_flat_1300.fits
        BPNAME=${ROOT}calib/dff_badpixelname_1300.fits
        idl -e sph_ifs_detector_flat_manual -args ${SOF} ${FFNAME} ${BPNAME}
        if [ $? -ne 0 ]; then
    	    echo "IDL returned an error value. Stopping here..."
    	    exit 1
        fi
        
        #
        # 1550 nm flat
        #
        if [ ${MODE} = 'YJH' ]; then
	    SOF=${ROOT}sof/flat_1550.sof
	    
	    if [ -f ${SOF} ]; then rm ${SOF}; fi
	    for f in ${flat_1300_files[*]}
	    do
	        echo "${ROOT}raw/${f}.fits         IFS_DETECTOR_FLAT_FIELD_RAW" >> ${SOF}
	    done
	    echo "${ROOT}calib/dark_bpm_cal.fits   IFS_STATIC_BADPIXELMAP"      >> ${SOF}
	    echo "${ROOT}calib/dark_bpm_cen.fits   IFS_STATIC_BADPIXELMAP"      >> ${SOF}
	    echo "${ROOT}calib/dark_bpm_cor.fits   IFS_STATIC_BADPIXELMAP"      >> ${SOF}
	    
	    # esorex --msg-level=debug sph_ifs_master_detector_flat \
	    #        --ifs.master_detector_flat.save_addprod=TRUE \
	    #        --ifs.master_detector_flat.outfilename=${ROOT}calib/master_detector_flat_1550.fits \
	    #        --ifs.master_detector_flat.lss_outfilename=${ROOT}calib/large_scale_flat_1550.fits \
	    #        --ifs.master_detector_flat.preamp_outfilename=${ROOT}calib/preamp_flat_1550.fits \
	    #        --ifs.master_detector_flat.badpixfilename=${ROOT}calib/dff_badpixelname_1550.fits \
	    #        --ifs.master_detector_flat.lambda=1.550 \
	    #        --ifs.master_detector_flat.smoothing_length=10.0 \
	    #        --ifs.master_detector_flat.smoothing_method=1 \
	    #        ${SOF}
	    
	    FFNAME=${ROOT}calib/master_detector_flat_1550.fits
	    BPNAME=${ROOT}calib/dff_badpixelname_1550.fits
	    idl -e sph_ifs_detector_flat_manual -args ${SOF} ${FFNAME} ${BPNAME}
	    if [ $? -ne 0 ]; then
    	        echo "IDL returned an error value. Stopping here..."
    	        exit 1
	    fi
	    
        fi
        
    fi

    #
    # SPEC POS
    #
    if [ $DO_SPEC_POS -ne 0 ]; then
        SOF=${ROOT}sof/specpos.sof
        
        if [ -f ${SOF} ]; then rm ${SOF}; fi
        for f in ${spec_files[*]}
        do
	    echo "${ROOT}raw/${f}.fits         IFS_SPECPOS_RAW"        >> ${SOF}
        done
        echo "${ROOT}calib/dark_cal.fits       IFS_MASTER_DARK"        >> ${SOF}
        
        if [ ${MODE} = 'YJH' ]; then Hmode=TRUE; else Hmode=FALSE; fi
        
        esorex --msg-level=debug sph_ifs_spectra_positions \
	       --ifs.spectra_positions.outfilename=${ROOT}calib/spectra_positions.fits \
	       --ifs.spectra_positions.hmode=${Hmode} \
	       ${SOF}
        
    fi

    #
    # WAVE CAL
    #
    if [ $DO_WAVE_CAL -ne 0 ]; then
        SOF=${ROOT}sof/wave.sof
        
        if [ -f ${SOF} ]; then rm ${SOF}; fi
        for f in ${wave_files[*]}
        do
	    echo "${ROOT}raw/${f}.fits               IFS_WAVECALIB_RAW"      >> ${SOF}
        done
        echo "${ROOT}calib/spectra_positions.fits    IFS_SPECPOS"            >> ${SOF}
        echo "${ROOT}calib/dark_cal.fits             IFS_MASTER_DARK"        >> ${SOF}
        
        if [ ${MODE} = 'YJ' ]; then
	    esorex --msg-level=debug sph_ifs_wave_calib \
	           --ifs.wave_calib.number_lines=3 \
	           --ifs.wave_calib.outfilename=${ROOT}calib/wave_calib.fits \
	           --ifs.wave_calib.wavelength_line1=0.9877 \
	           --ifs.wave_calib.wavelength_line2=1.1237 \
	           --ifs.wave_calib.wavelength_line3=1.3094 \
	           ${SOF}
        else	
	    esorex --msg-level=debug sph_ifs_wave_calib \
	           --ifs.wave_calib.number_lines=4 \
	           --ifs.wave_calib.outfilename=${ROOT}calib/wave_calib.fits \
	           --ifs.wave_calib.wavelength_line1=0.9877 \
	           --ifs.wave_calib.wavelength_line2=1.1237 \
	           --ifs.wave_calib.wavelength_line3=1.3094 \
	           --ifs.wave_calib.wavelength_line4=1.5451 \
	           ${SOF}
        fi
    fi

    #
    # IFU FLAT
    #
    if [ $DO_IFU_FLAT -ne 0 ]; then
        SOF=${ROOT}sof/ifu_flat.sof
        
        if [ -f ${SOF} ]; then rm ${SOF}; fi
        for f in ${ifu_files[*]}
        do
	    echo "${ROOT}raw/${f}.fits                            IFS_FLAT_FIELD_RAW"     >> ${SOF}
        done
        echo "${ROOT}calib/wave_calib.fits                        IFS_WAVECALIB"          >> ${SOF}    
        echo "${ROOT}calib/master_detector_flat_1020_l1.fits      IFS_MASTER_DFF_LONG1"   >> ${SOF}
        echo "${ROOT}calib/master_detector_flat_1230_l2.fits      IFS_MASTER_DFF_LONG2"   >> ${SOF}
        echo "${ROOT}calib/master_detector_flat_1300_l3.fits      IFS_MASTER_DFF_LONG3"   >> ${SOF}
        if [ ${MODE} = 'YJH' ]; then
    	    echo "${ROOT}calib/master_detector_flat_1550_l3.fits  IFS_MASTER_DFF_LONG4"   >> ${SOF}
        fi
        echo "${ROOT}calib/master_detector_flat_white_l5.fits     IFS_MASTER_DFF_LONGBB"  >> ${SOF}
        echo "${ROOT}calib/master_detector_flat_white_l5.fits     IFS_MASTER_DFF_SHORT"   >> ${SOF}
        echo "${ROOT}calib/dark_cal.fits                          IFS_MASTER_DARK"        >> ${SOF}
        
        esorex sph_ifs_instrument_flat \
	       --ifs.instrument_flat.ifu_filename=${ROOT}calib/ifu_flat.fits \
	       --ifs.instrument_flat.nofit=TRUE \
	       ${SOF}

    fi
    
fi
    
#
# PRE-PROCESSING
#
if [ $DO_PREPROC -ne 0 ]; then
    SOF=${ROOT}sof/preproc.sof

    #
    # PSF
    #
    if [ -n "$sci_files_psf" ]; then
	if [ -f ${SOF} ]; then rm ${SOF}; fi
	for f in ${sci_files_psf[*]}
	do
    	    echo "${ROOT}raw/${f}.fits            IFS_RAW"                >> ${SOF}
	done
	if [ -f ${ROOT}calib/sky_psf.fits ] && [ $USE_SKY -ne 0 ]; then
    	    echo "${ROOT}calib/sky_psf.fits       IFS_MASTER_DARK"        >> ${SOF}
    	    echo "${ROOT}calib/sky_bpm_psf.fits   IFS_STATIC_BADPIXELMAP" >> ${SOF}
	else 
    	    echo "${ROOT}calib/dark_psf.fits      IFS_MASTER_DARK"        >> ${SOF}
    	    echo "${ROOT}calib/dark_bpm_psf.fits  IFS_STATIC_BADPIXELMAP" >> ${SOF}
	fi
	
	# always collapse off-axis PSF
	idl -e sph_ifs_preprocess -args ${SOF} 1 ${DO_BKG_SUB} ${DO_BADPIXEL} ${DO_CROSSTALK} MEAN 0.0 0.0
	if [ $? -ne 0 ]; then
    	    echo "IDL returned an error value. Stopping here..."
    	    exit 1
	fi
    fi
    
    #
    # star center
    #
    if [ -n "$sci_files_cen" ]; then
	if [ -f ${SOF} ]; then rm ${SOF}; fi
	for f in ${sci_files_cen[*]}
	do
    	    echo "${ROOT}raw/${f}.fits            IFS_RAW"                >> ${SOF}
	done
	if [ -f ${ROOT}calib/sky_cen.fits ] && [ $USE_SKY -ne 0 ]; then
    	    echo "${ROOT}calib/sky_cen.fits       IFS_MASTER_DARK"        >> ${SOF}
    	    echo "${ROOT}calib/sky_bpm_cen.fits   IFS_STATIC_BADPIXELMAP" >> ${SOF}
	else       
    	    echo "${ROOT}calib/dark_cen.fits      IFS_MASTER_DARK"        >> ${SOF}
    	    echo "${ROOT}calib/dark_bpm_cen.fits  IFS_STATIC_BADPIXELMAP" >> ${SOF}
	fi

	if [ $WAFFLE_OBS -ne 0 ]; then
	    # if waffle was used on all frames on purpose, then process the data as normal science
	    idl -e sph_ifs_preprocess -args ${SOF} ${DO_COLLAPSE} ${DO_BKG_SUB} ${DO_BADPIXEL} ${DO_CROSSTALK} ${COLLAPSE_TYPE} ${COLLAPSE_VAL} ${COLLAPSE_TOL}
	else
	    # otherwise collapse cube for better SNR
	    idl -e sph_ifs_preprocess -args ${SOF} 1 ${DO_BKG_SUB} ${DO_BADPIXEL} ${DO_CROSSTALK} MEAN 0.0 0.0
	fi
	
	if [ $? -ne 0 ]; then
    	    echo "IDL returned an error value. Stopping here..."
    	    exit 1
	fi
    fi
	
    #
    # corono data
    #
    if [ -n "$sci_files_cor" ]; then
	if [ -f ${SOF} ]; then rm ${SOF}; fi
	for f in ${sci_files_cor[*]}
	do
    	    echo "${ROOT}raw/${f}.fits            IFS_RAW"                >> ${SOF}
	done
	if [ -f ${ROOT}calib/sky_cor.fits ] && [ $USE_SKY -ne 0 ]; then
    	    echo "${ROOT}calib/sky_cor.fits       IFS_MASTER_DARK"        >> ${SOF}
    	    echo "${ROOT}calib/sky_bpm_cor.fits   IFS_STATIC_BADPIXELMAP" >> ${SOF}
	else
    	    echo "${ROOT}calib/dark_cor.fits      IFS_MASTER_DARK"        >> ${SOF}
    	    echo "${ROOT}calib/dark_bpm_cor.fits  IFS_STATIC_BADPIXELMAP" >> ${SOF}
	fi
	
	idl -e sph_ifs_preprocess -args ${SOF} ${DO_COLLAPSE} ${DO_BKG_SUB} ${DO_BADPIXEL} ${DO_CROSSTALK} ${COLLAPSE_TYPE} ${COLLAPSE_VAL} ${COLLAPSE_TOL}
	if [ $? -ne 0 ]; then
    	    echo "IDL returned an error value. Stopping here..."
    	    exit 1
	fi
    fi
    
    #
    # move to preprocessed directory
    #
    mv ${ROOT}raw/*preproc*.fits ${ROOT}interm/    
fi

#
# science
#
if [ $DO_SCIENCE -ne 0 ]; then

    #
    # science target
    #
    if [ $DO_SCI_CORO -ne 0 ]; then
        SOF=${ROOT}sof/science.sof
        
        #
        # PSF
        #
        if [ -n "$sci_files_psf" ]; then
            if [ -f ${SOF} ]; then rm ${SOF}; fi
            for f in ${sci_files_psf[*]}
            do
    	        echo "${ROOT}interm/${f}_preproc${SUFFIX_PSF}.fits      IFS_SCIENCE_DR_RAW"     >> ${SOF}
            done
            echo "${ROOT}calib/master_detector_flat_1020_l1.fits        IFS_MASTER_DFF_LONG1"   >> ${SOF}
            echo "${ROOT}calib/master_detector_flat_1230_l2.fits        IFS_MASTER_DFF_LONG2"   >> ${SOF}
            echo "${ROOT}calib/master_detector_flat_1300_l3.fits        IFS_MASTER_DFF_LONG3"   >> ${SOF}
            if [ ${MODE} = 'YJH' ]; then
    	        echo "${ROOT}calib/master_detector_flat_1550_l3.fits    IFS_MASTER_DFF_LONG4"   >> ${SOF}
            fi
            echo "${ROOT}calib/master_detector_flat_white_l5.fits       IFS_MASTER_DFF_LONGBB"  >> ${SOF}
            echo "${ROOT}calib/master_detector_flat_white_l5.fits       IFS_MASTER_DFF_SHORT"   >> ${SOF}
            echo "${ROOT}calib/dff_badpixelname_white_l5.fits           IFS_STATIC_BADPIXELMAP"  >> ${SOF}
            echo "${ROOT}calib/wave_calib.fits                          IFS_WAVECALIB"          >> ${SOF}

            #new method: IFU flat done later
            echo "${ROOT}calib/ifu_flat.fits                            IFS_IFU_FLAT_FIELD"     >> ${SOF}

            esorex --msg-level=debug sph_ifs_science_dr \
    	           --ifs.science_dr.use_adi=0 \
    	           --ifs.science_dr.spec_deconv=FALSE \
    	           ${SOF}
	    
            if [ $DO_POSTPROC -ne 0 ]; then
    	        idl -e sph_ifs_postprocess -args ${SOF}
    	        if [ $? -ne 0 ]; then
    		    echo "IDL returned an error value. Stopping here..."
    		    exit 1
    	        fi
            fi
        fi
        
        #
        # corono + star center data
        #
        if [ -n "$sci_files_cen" ] || [ -n sci_files_cor ]; then
	    if [ -f ${SOF} ]; then rm ${SOF}; fi
	    for f in ${sci_files_cor[*]}
	    do
    	        echo "${ROOT}interm/${f}_preproc${SUFFIX}.fits          IFS_SCIENCE_DR_RAW"     >> ${SOF}
	    done
	    for f in ${sci_files_cen[*]}
	    do
	        if [ $WAFFLE_OBS -ne 0 ]; then
		    # if waffle was used on all frames on purpose, use science suffix
		    echo "${ROOT}interm/${f}_preproc${SUFFIX}.fits      IFS_SCIENCE_DR_RAW"     >> ${SOF}
	        else
		    # otherwise use the suffix for collapsed data
		    echo "${ROOT}interm/${f}_preproc${SUFFIX_PSF}.fits  IFS_SCIENCE_DR_RAW"     >> ${SOF}
	        fi
	    done
	    echo "${ROOT}calib/master_detector_flat_1020_l1.fits        IFS_MASTER_DFF_LONG1"   >> ${SOF}
	    echo "${ROOT}calib/master_detector_flat_1230_l2.fits        IFS_MASTER_DFF_LONG2"   >> ${SOF}
	    echo "${ROOT}calib/master_detector_flat_1300_l3.fits        IFS_MASTER_DFF_LONG3"   >> ${SOF}
	    if [ ${MODE} = 'YJH' ]; then
    	        echo "${ROOT}calib/master_detector_flat_1550_l3.fits    IFS_MASTER_DFF_LONG4"   >> ${SOF}
	    fi
	    echo "${ROOT}calib/master_detector_flat_white_l5.fits       IFS_MASTER_DFF_LONGBB"  >> ${SOF}
	    echo "${ROOT}calib/master_detector_flat_white_l5.fits       IFS_MASTER_DFF_SHORT"   >> ${SOF}
	    echo "${ROOT}calib/dff_badpixelname_white_l5.fits           IFS_STATIC_BADPIXELMAP" >> ${SOF}
	    echo "${ROOT}calib/wave_calib.fits                          IFS_WAVECALIB"          >> ${SOF}

	    #new method: IFU flat done later
	    echo "${ROOT}calib/ifu_flat.fits                            IFS_IFU_FLAT_FIELD"     >> ${SOF}
	    
	    esorex --msg-level=debug sph_ifs_science_dr \
    	           --ifs.science_dr.use_adi=0 \
    	           --ifs.science_dr.spec_deconv=FALSE \
    	           ${SOF}

	    if [ $DO_POSTPROC -ne 0 ]; then
    	        idl -e sph_ifs_postprocess -args ${SOF}
    	        if [ $? -ne 0 ]; then
    		    echo "IDL returned an error value. Stopping here..."
    		    exit 1
    	        fi
	    fi
        fi
        
        #
        # move products to final directory
        #
        mv ${ROOT}SPHER*${SUFFIX}*.fits ${ROOT}products/
        mv ${ROOT}SPHER*${SUFFIX_PSF}*.fits ${ROOT}products/
        cp ${ROOT}interm/SPHER*info.fits ${ROOT}products/
    fi

    #
    # wavelength calibration
    #
    if [ $DO_SCI_LMBD -ne 0 ]; then
        #
        # wavelength calib collapse
        #
        SOF=${ROOT}sof/preproc.sof

        if [ -f ${SOF} ]; then rm ${SOF}; fi
        for f in ${wave_files[*]}
        do
    	    echo "${ROOT}raw/${f}.fits          IFS_RAW"                >> ${SOF}
        done
        # for f in ${ifu_files[*]}
        # do
        #     echo "${ROOT}raw/${f}.fits          IFS_RAW"                >> ${SOF}
        # done
        echo "${ROOT}calib/dark_cal.fits        IFS_MASTER_DARK"        >> ${SOF}
        echo "${ROOT}calib/dark_bpm_psf.fits    IFS_STATIC_BADPIXELMAP" >> ${SOF}

        # collapse, no background subtraction, bad pixel correction, crosstalk correction
        idl -e sph_ifs_preprocess -args ${SOF} 1 0 1 1 MEAN 0.0 0.0
        if [ $? -ne 0 ]; then
    	    echo "IDL returned an error value. Stopping here..."
    	    exit 1
        fi
        
        #
        # move to preprocessed directory
        #
        mv ${ROOT}raw/*preproc*.fits ${ROOT}interm/

        #
        # wavelength calib <==> science
        #
        SOF=${ROOT}sof/science_lmbd.sof
        
        if [ -f ${SOF} ]; then rm ${SOF}; fi
        for f in ${wave_files[*]}
        do
    	    echo "${ROOT}interm/${f}_preproc_col_mean_bp_ct.fits    IFS_SCIENCE_DR_RAW"     >> ${SOF}
        done
        # for f in ${ifu_files[*]}
        # do
        #     echo "${ROOT}interm/${f}_preproc_col_mean_bp_ct.fits    IFS_SCIENCE_DR_RAW"     >> ${SOF}
        # done
        echo "${ROOT}calib/master_detector_flat_1020_l1.fits        IFS_MASTER_DFF_LONG1"   >> ${SOF}
        echo "${ROOT}calib/master_detector_flat_1230_l2.fits        IFS_MASTER_DFF_LONG2"   >> ${SOF}
        echo "${ROOT}calib/master_detector_flat_1300_l3.fits        IFS_MASTER_DFF_LONG3"   >> ${SOF}
        if [ ${MODE} = 'YJH' ]; then
	    echo "${ROOT}calib/master_detector_flat_1550_l3.fits    IFS_MASTER_DFF_LONG4"   >> ${SOF}
        fi
        echo "${ROOT}calib/master_detector_flat_white_l5.fits       IFS_MASTER_DFF_LONGBB"  >> ${SOF}
        echo "${ROOT}calib/master_detector_flat_white_l5.fits       IFS_MASTER_DFF_SHORT"   >> ${SOF}
        echo "${ROOT}calib/dff_badpixelname_white_l5.fits           IFS_STATIC_BADPIXELMAP" >> ${SOF}
        echo "${ROOT}calib/wave_calib.fits                          IFS_WAVECALIB"          >> ${SOF}
        echo "${ROOT}calib/ifu_flat.fits                            IFS_IFU_FLAT_FIELD"     >> ${SOF}

        esorex --msg-level=debug sph_ifs_science_dr \
    	       --ifs.science_dr.use_adi=0 \
    	       --ifs.science_dr.spec_deconv=FALSE \
    	       ${SOF}

        if [ $DO_POSTPROC -ne 0 ]; then
    	    idl -e sph_ifs_postprocess -args ${SOF}
    	    if [ $? -ne 0 ]; then
    	        echo "IDL returned an error value. Stopping here..."
    	        exit 1
    	    fi
        fi

        #
        # move products to final directory
        #
        mv ${ROOT}SPHER*_mean_bp_ct*.fits ${ROOT}products/
    fi

fi
