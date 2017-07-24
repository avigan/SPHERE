;+
; NAME:
;
;   SPH_IRD_LSS
;
; PURPOSE:
;
;   Analysis package for IRDIS/LSS data
;
; CALLING SEQUENCE:
;
;  From the IDL command prompt:
;
;    IDL> .run sph_ird_lss
;
; DESCRIPTION:
;
;  This package is dedicated to the data reduction and analysis of
;  IRDIS long-slit spectroscopy data. It is made of several modules
;  that have to be executed sequentially, starting from the creation
;  of the static calibrations, the calibration and cleaning of the
;  data, the speckle subtraction, and finally the spectral extraction.
;
;  *
;  * Setting-up the analysis
;  *
;
;  You will find below a typical setup for an analysis. The essential
;  parameters are:
;
;    - aroot: path to the working directory. In this directory, the
;      raw data must be in a "raw/" subdirectory
;    - resol: specify LRS or MRS
;    - suffix: optional suffix that will be appended to the files created by the pipeline
;    - pla_sep: separation of the planet in the data in as. This value
;      should be negative if in the OB you entered the on-sky PA of
;      the companion. For example, if the separation is 0.5" and the
;      on-sky PA is 100°, and that you have entered 100° in the POSANG
;      keyword of the OB, then you should entered -0.5 for pla_sep
;    - files: enter the raw file names (without the .fits
;      extension). The wavelength calibration expects only one file,
;      while the others can be a list of files 
;
;  Once you have set up the analysis, you start the analysis. 
;
;  *
;  * static calibrations
;  *
;
;  You execute this part by setting the keyword do_calib to 1 and then
;  executing the code.
;
;  This module calls esorex and specific recipes to create all the
;  static calibrations: dark/background, flat, wavelength calibration,
;  bad pixel maps. They are all saved in the "calib/"
;  subdirectory. You should check the quality of the created
;  calibrations before going further.
;  
;  Note for Mac OS X 10.11: for esorex to be able to find the recipes,
;  you have to set properly the DYLD_LIBRARY_PATH environment variable
;  in the program. This is done near line 460: 
;
;     env = "export DYLD_LIBRARY_PATH=$HOME/ESO/lib:$DYLD_LIBRARY_PATH && "
;
;  You have to replace $HOME/ESO/lib with the path to the ESO lib directory.
;
;  *
;  * Calibration and cleaning of the data
;  *
;
;  You execute this part by setting the keyword do_clean to 1 and then
;  executing the code. 
;
;  This module calibrates the science data, cleans it and save it for
;  future analysis. During the step, you have the choice to collapse
;  each science data cube into a single image (faster analysis) or to
;  keep all the individual DITs. To keep all the DITs, set the all_DIT
;  keyword to 1.
;
;  *
;  * Registration and combination of the data
;  *
;
;  You execute this part by setting the keyword do_combine to 1 and
;  then executing the code.
;
;  This module reads the calibrated data created in the previous step,
;  aligns them and combines them into a single file that is saved into
;  the "products/" subdirectory. The off-axis PSF, coronagraphic data
;  and telluric standard are saved independently. The 2D wavelength
;  map and 1D wavelength vector are also saved into files.
;
;  The off-axis PSF data (science target and telluric standard) are
;  aligned by performing a Gaussian fit to find the position of the
;  star inside the frame. The coronagraphic data are aligned using an
;  FFT-based cross-correlation to find a possible shift between the
;  different files. This procedure can be biased if the companion is
;  really bright and you will need to mask it manually (search for the
;  warning message around lines 1090 and 1110).
;  
;  This module also allows measuring and compensating a chromatic
;  effect in IRDIS that significantly affects the MRS data (but not
;  the LRS data). The visible result is that the PSF is shifted with
;  wavelength. This effect can be corrected by setting the keyword
;  corr_chroma to 1 before executing this module. During this step, a
;  window showing the position of the PSF as a function if it location
;  on the detector (for both IRDIS fields) will appear. It will also
;  show a fit (4th order in LRS, linear in MRS) that will be used to
;  correcte the effect.
;
;  *
;  * Speckle subtraction
;  *
;
;  You execute this part by setting the keyword do_specsub to 1 and
;  then executing the code.

;  This module performs the speckles subtraction. Several algorithms
;  are available at the moment (others are in development but not
;  robust enough to be distributed):
;  
;    - simplesub: subtraction of the speckles using the symmetric with
;      respect to the star (without spatial rescaling) 
;    - symmetry: subtraction of the speckles using the symmetric with
;      respect to the star (with spatial rescaling) 
;    - vigan2008: subtraction of the speckles using the method
;      described in Vigan et al. (2008) (with spatial rescaling) 
;    - locs: LOCI-type subtraction of the speckles (with spatial rescaling)
;    
;  The pipeline is very flexible for the addition of new
;  speckle-subtraction methods. If you want to implement a new method
;  you can get in touch with me or look at the methods already
;  implemented.

;  For this module you have to specify which of the two IRDIS fields
;  you are considering. This is done by setting the field variable at
;  the beginning of the program to 'F0' or 'F1' (left or right
;  field).
;  
;  There are some user parameters for the vigan2008 and locs
;  methods. The main one is the planet_size parameter, which the size
;  of a mask that will be applied at the position of the planet to
;  avoid fitting its signal. It is given in unite of lambda/D, and its
;  value is defined in the "parameters for speckle subtraction" part
;  at the beginning of the code. For locs, the locs_crit parameter
;  defines an exclusion zone for the fit of the speckles, also
;  expressed in lambda/D. It equivalent to Nδ in the LOCI reference
;  paper of Lafrenière et al. (2007). A value between 2 and 3 usually
;  gives good results. The scale parameter corresponds to an
;  over-sampling factor along the spatial dimension: tests have showed
;  that a value of 2 gives better results than 1 (i.e. than no
;  oversampling).
;  
;  At the end of this module, a final dataset is saved in the
;  "products/". It contains the off-axis PSF, the collapsed
;  coronagraphic data and the speckle-subtracted data. 
;
;  *
;  * Spectral extraction
;  *
;
;  You execute this part by setting the keyword do_extract to 1 and
;  then executing the code.
;  
;  This module reads the final products and extracts the spectra of
;  the star and planet using aperture photometry in each spectral
;  channel. Each spectrum is corrected for the effect of any neutral
;  density filter that was in the optical path.
;
;  At the end it produces two plots that show to the extracted spectra
;  (PSF, halo, halo+planet, planet, noise) and the spectrum of the
;  planet calibrated in contrast. A text file containing the spectra
;  extracted for the PSF, the planet and the noise is also saved on
;  disk.
;
;  *
;  * Fake planet injection
;  *
;
;  The code enables the possibility to inject the spectrum of a fake
;  planet inside the data before performing the analysis. This is
;  enabled by the PLA_FAKE keyword. The spectrum is based on the
;  stellar spectrum obtained from the off-axis reference PSF. The user
;  can choose the peak contrast in H-band and the H-J color.
;
;  Unless directly modified by the user, the fake planet is injected
;  at a position symmetric to the real companion with respect to the
;  star, i.e. if pla_sep = -0.35", the fake planet will the injected
;  on the other side of the star at a separation of 0.35". Note that
;  in the case where a fake planet is injected, the speckle
;  subtraction will be optimised for the fake, and NOT for the real
;  companion. The files for the fake planet are saved independently to
;  allow separate analysis. The spectral extraction will allow you to
;  retrieved the spectrum of the fake planet and compare it to the
;  input spectrum. 
;
;
; INPUTS:
;
;  General inputs (independent of the target):
;
;    DO_CALIB - keyword to create static calibrations
;    
;    DO_CLEAN - keyword to calibrate and clean the science data
;    
;    DO_COMBINE - keyword to align and combine the science data
;    
;    DO_WAVECOR - keyword to recalibrate the wavelength using sky lines
;    
;    DO_SPECSUB - keyword to perform speckle subtraction
;    
;    DO_EXTRACT - keyword to perform spectral extraction and save
;                 final data
;
;    ALL_DIT - keyword to save and analyse all DITs independently
;
;    CORR_CHROMA - keyword to measured and correct chromatic
;                  distortion (essential in MRS data)
;
;    FIELD - name of the IRDIS field to be analysed (F0 or F1)
;
;    METHOD - name of the speckle subtraction method
;
;    EPS - keyword to save all figures in EPS files
;
;    PLA_FAKE - keyword to enable injection of a fake planet (see
;               documentation for details)
;
;    PLA_CNT_FAKE - contrast of the fake planet in H-band, in mag
;
;    PLA_COL_FAKE - H-J color of the fake planet
;
;
;  Inputs for the analysis of a specific target:
;
;    AROOT - path to the working directory
;
;    SUFFIX - suffix that will be added to the produced files
;
;    PLA_SEP - angular separation of the planet, in as
;
;    SCI_FILES - science data files
;
;    SCI_BKG_FILES - background files for science data
;
;    PSF_FILES - off-axis PSF files
;
;    PSF_BKG_FILES - background files for off-axis PSF
;
;    CAL_FILES - telluric standard files
;
;    CAL_BKG_FILES - background files for telluric standard
;
;    WAVE_FILES - wavelength calibration file
;
;    WAVE_BKG_FILE - background file for wavelength calibration
;
;    FLAT_FILES - flat field files
;
; OUTPUTS:
;
;  Many outputs produced at different steps by the pipeline. See
;  documentation for details.
;
; REQUIRED PROGRAMS AND LIBRARIES:
;
;  The following IDL generic libraries are necessary:
;
;    * astronomy user's library: http://idlastro.gsfc.nasa.gov/
;    * MPFIT: https://www.physics.wisc.edu/~craigm/idl/fitting.html
;
;  The following SPHERE-specific library is also necessary:
;
;    * SPHERE transmission: http://astro.vigan.fr/tools.html
;
;  To create the standard calibrations, it is also required that you
;  install the officiel sphere pipeline from ESO:
;
;    * http://www.eso.org/sci/software/pipelines/
;
; AUTHOR:
; 
;   Arthur Vigan
;   Laboratoire d'Astrophysique de Marseille
;   arthur.vigan@lam.fr
;
; LICENSE:
;
;   This code is release under the MIT license. The full text of the
;   license is included in a separate file LICENSE.txt.
;
;   The developement of the SPHERE instrument has demanded a
;   tremendous effort from many scientists, who have devoted several
;   years of their life to design, build, test and commission this new
;   instrument. To recognize this work, we kindly ask you to cite the
;   relevant papers in your scientific work. More specifically,
;   because this code is dedicated to the SPHERE/IRDIS subsystem,
;   please cite the two following papers:
;
;    * IRDIS general descripton: Dohlen et al., 2008, SPIE, 7014
;    * Long-Slit Spectroscopy mode: Vigan et al., 2008, A&A, 489, 1345
;
;   In addition, this analysis code is available on the Astrophysics
;   Source Code Library and should be cited explicitely:
;
;    * 
;
;   We are grateful for your effort, and hope that this tool will
;   contribute to your scientific work and discoveries. Please feel
;   free to report any bug or possible improvement to the author(s)
;
; MODIFICATION HISTORY:
;
;   arthur.vigan - v1.2 - 02/2016
;                  fix calibrations for LRS mode, improved
;                  documentation, added support for fake planet
;                  simulation
;
;   arthur.vigan - v1.1 - 01/2016
;                  bug fix for chromaticity correction, bug fix for
;                  all DIT analysis, improved cosmetics, moved
;                  parameters for speckle subtraction at the
;                  beginning, save cleaned science data in new
;                  sub-directory "interm", improved cosmetics
;
;   arthur.vigan - v1.0 - 01/2016
;                  first public version based on AIT and commissioning
;                  tools
;                            
;-

do_calib    = 0     ;; create static calibrations
do_clean    = 0     ;; clean/calibrate data
do_combine  = 0     ;; combine data
do_specsub  = 0     ;; speckle subtraction
do_extract  = 0     ;; spectral extraction

all_DIT     = 0     ;; save and analyze all DITs
corr_chroma = 0     ;; measure and correct chromaticity (in do_combine step)

field = 'F0'        ;; field to be analyzed (F0, F1)

;; speckle subtraction method
;method = 'simplesub'
;method = 'vigan2008'
;method = 'symmetry'
method  = 'locs'

;; parameters for speckle subtraction
pla_size   = 2.5
locs_crit  = 1.5
nmodes_pca = 1

eps = 0

;; fake planet parameters
pla_fake = 0
pla_cnt_fake = 12.5       ;; contrast in H-band [mag]
pla_col_fake = 0.2        ;; H-J color [mag]


;; ----------------------------------------------------------
;; 2014-02-03 - GTO - BetaPic - MR_WL
;;
aroot = '~/data/SPHERE/LSS/TargetName/'

resol   = 'LRS'     ;; resolution: LRS or MRS
suffix  = ''        ;; optional file suffix
pla_sep = -0.335    ;; angular separation of the planet in as

sci_files = ['...']
sci_bkg_files = ['...']

psf_files = ['...']
psf_bkg_files = ['...']

cal_files = ['...']
cal_bkg_files = ['...']

wave_files = '...'
wave_bkg_file = '...'

flat_files = ['...']

;; -----------------------------------------------------------------------------------
;; -----------------------------------------------------------------------------------
;; -----------------------------------------------------------------------------------
postfix = ''

;; subfields
if (resol eq 'LRS') then begin
   c0 = [484+15,496]
   c1 = [1512+15,486]
   extx = 400;250
   exty = 200;90
endif else if (resol eq 'MRS') then begin
   c0 = [476+15,519]
   c1 = [1503+15,510]
   extx = 250
   exty = 390
endif else message,'Unknown resolution "'+resol+'"'

x0_l = c0[0]-extx
y0_l = c0[1]-exty
x1_l = c0[0]+extx-1
y1_l = c0[1]+exty-1

x0_r = c1[0]-extx
y0_r = c1[1]-exty
x1_r = c1[0]+extx-1
y1_r = c1[1]+exty-1

taillex = 2*extx
tailley = 2*exty

file_mkdir,aroot+'tmp/'
file_mkdir,aroot+'calib/'
file_mkdir,aroot+'interm/'
file_mkdir,aroot+'products/'

;; fake planet
if keyword_set(pla_fake) then begin
   pla_sep = -pla_sep
   postfix = 'fake'
endif

pixel = 0.01225
if (suffix  ne '') then suffix  = '_'+suffix
if (postfix ne '') then postfix = '_'+postfix

if keyword_set(do_calib) then begin
   ;;
   ;; create basic calibrations
   ;;

   ;; Mac OS X > 10.11
   env = "export DYLD_LIBRARY_PATH=$HOME/ESO/lib:$DYLD_LIBRARY_PATH && "

   ;; flat
   openw,lun,aroot+'tmp/tmp.sof',/get_lun,width=1000
   printf,lun,transpose(aroot+'raw/'+flat_files+'.fits   IRD_FLAT_FIELD_RAW')
   printf,lun
   free_lun,lun

   esorex = "esorex sph_ird_instrument_flat "+ $
            "--ird.instrument_flat.save_addprod=TRUE "+ $
            "--ird.instrument_flat.outfilename='"+aroot+"calib/flat.fits' "+ $
            "--ird.instrument_flat.badpixfilename='"+aroot+"calib/flat_bpm.fits' "+ $
            aroot+"tmp/tmp.sof"
   cmd = env+"cd "+aroot+"tmp/ && "+esorex
   
   spawn,cmd
   
   data = readfits(aroot+'calib/flat.fits',hdr)
   writefits,aroot+'calib/flat.fits',data,hdr

   ;; wavelength calibration
   if (resol eq 'MRS') then begin
      ;; hides the second order in the data
      data = readfits(aroot+'raw/'+wave_files+'.fits',hdr)
      data[*,0:60] = 0
      writefits,aroot+'raw/'+wave_files+'_modified.fits',data,hdr
      wave_files = wave_files+'_modified'

      grism  = 'TRUE'
      thresh = '1000'
      nlines = '5'
   endif else if (resol eq 'LRS') then begin
      grism  = 'FALSE'
      thresh = '2000'
      nlines = '6'
   endif else stop
   
   openw,lun,aroot+'tmp/tmp.sof',/get_lun,width=1000
   printf,lun,aroot+'raw/'+wave_files+'.fits   IRD_WAVECALIB_RAW'
   printf,lun,aroot+'raw/'+wave_bkg_file+'.fits   IRD_MASTER_DARK'
   printf,lun,aroot+'calib/flat.fits   IRD_FLAT_FIELD'
   printf,lun,aroot+'calib/flat_bpm.fits   IRD_STATIC_BADPIXELMAP'
   printf,lun
   free_lun,lun

   esorex = "esorex sph_ird_wave_calib "+ $
            "--ird.wave_calib.column_width=200 "+ $
            "--ird.wave_calib.grism_mode="+grism+" "+ $
            "--ird.wave_calib.threshold="+thresh+" "+ $
            "--ird.wave_calib.number_lines="+nlines+" "+ $
            "--ird.wave_calib.outfilename='"+aroot+"calib/wave_cal.fits' "+ $
            aroot+"tmp/tmp.sof"
   cmd = env+"cd "+aroot+"tmp/ && "+esorex
   
   spawn,cmd

   data = readfits(aroot+'calib/wave_cal.fits',hdr)
   writefits,aroot+'calib/wave_cal.fits',data,hdr
   
   ;; science background/dark
   openw,lun,aroot+'tmp/tmp.sof',/get_lun,width=1000
   printf,lun,transpose(aroot+'raw/'+sci_bkg_files+'.fits   IRD_DARK_RAW')
   printf,lun
   free_lun,lun

   esorex = "esorex sph_ird_master_dark "+ $
            "--ird.master_dark.save_addprod=TRUE "+ $
            "--ird.master_dark.outfilename='"+aroot+"calib/dark_sci.fits' "+ $
            "--ird.master_dark.badpixfilename='"+aroot+"calib/dark_sci_bpm.fits' "+ $
            aroot+"tmp/tmp.sof"
   cmd = env+"cd "+aroot+"tmp/ && "+esorex
   
   spawn,cmd

   data = readfits(aroot+'calib/dark_sci.fits',hdr)
   writefits,aroot+'calib/dark_sci.fits',data,hdr
   
   ;; off-axis PSF background/dark
   openw,lun,aroot+'tmp/tmp.sof',/get_lun,width=1000
   printf,lun,transpose(aroot+'raw/'+psf_bkg_files+'.fits   IRD_DARK_RAW')
   printf,lun
   free_lun,lun

   esorex = "esorex sph_ird_master_dark "+ $
            "--ird.master_dark.save_addprod=TRUE "+ $
            "--ird.master_dark.outfilename='"+aroot+"calib/dark_psf.fits' "+ $
            "--ird.master_dark.badpixfilename='"+aroot+"calib/dark_psf_bpm.fits' "+ $
            aroot+"tmp/tmp.sof"
   cmd = env+"cd "+aroot+"tmp/ && "+esorex
   
   spawn,cmd

   data = readfits(aroot+'calib/dark_psf.fits',hdr)
   writefits,aroot+'calib/dark_psf.fits',data,hdr
   
   ;; telluric standard background/dark
   if (cal_files[0] ne '') then begin   
      openw,lun,aroot+'tmp/tmp.sof',/get_lun,width=1000
      printf,lun,transpose(aroot+'raw/'+cal_bkg_files+'.fits   IRD_DARK_RAW')
      printf,lun
      free_lun,lun

      esorex = "esorex sph_ird_master_dark "+ $
               "--ird.master_dark.save_addprod=TRUE "+ $
               "--ird.master_dark.outfilename='"+aroot+"calib/dark_cal.fits' "+ $
               "--ird.master_dark.badpixfilename='"+aroot+"calib/dark_cal_bpm.fits' "+ $
               aroot+"tmp/tmp.sof"
      cmd = env+"cd "+aroot+"tmp/ && "+esorex
      
      spawn,cmd

      data = readfits(aroot+'calib/dark_cal.fits',hdr)
      writefits,aroot+'calib/dark_cal.fits',data,hdr
   endif
endif

if keyword_set(do_clean) then begin
   ;; ------------------------------
   ;; read calibrations
   ;;
   print,'+ calibrations'
   flat = readfits(aroot+'calib/flat.fits')
   bkg  = readfits(aroot+'calib/dark_sci.fits')
   pbkg = readfits(aroot+'calib/dark_psf.fits')
   cbkg = readfits(aroot+'calib/dark_cal.fits')
   
   bpm_files = file_search(aroot+'calib/*bpm*')
   bpm  = fltarr(size(flat,/dim))
   for f=0,n_elements(bpm_files)-1 do begin
      tmp = readfits(bpm_files[f])
      bpm = bpm or tmp
   endfor
   bpm = byte(abs(1-bpm))

   wave = readfits(aroot+'calib/wave_cal.fits')

   ;; ------------------------------
   ;; extract sub-fields
   ;;
   wave0 = wave[x0_l:x1_l,y0_l:y1_l]
   wave1 = wave[x0_r:x1_r,y0_r:y1_r]   
   
   bpm0 = bpm[x0_l:x1_l,y0_l:y1_l]
   bpm1 = bpm[x0_r:x1_r,y0_r:y1_r]

   ;; corono flat
   ;; cflat = readfits(aroot+'calib/coro_flat.fits')
   ;; dim = size(coro_flat0,/dim)
   ;; coro_flat = fltarr(dim[0],dim[1],2)
   ;; coro_flat[*,*,0] = cflat[x0_l:x1_l,y0_l:y1_l]
   ;; coro_flat[*,*,1] = cflat[x0_r:x1_r,y0_r:y1_r]
   ;; writefits,aroot+'interm/coro_flat.fits',coro_flat,/compress
   
   ;; ------------------------------
   ;; reference PSF
   ;;
   print,'+ reference PSF data'

   ncube = n_elements(psf_files)
   for n=0,ncube-1 do begin
      print,' + cube '+numformat(n+1)+' / '+numformat(ncube)

      cube = readfits(aroot+'raw/'+psf_files[n]+'.fits',hdr)
      dit  = sxpar_eso(hdr,'HIERARCH ESO DET SEQ1 DIT')
      ndit = sxpar_eso(hdr,'HIERARCH ESO DET NDIT')
      rom  = sxpar_eso(hdr,'HIERARCH ESO DET READ CURNAME')      

      ;; collapse cube if needed
      s = size(cube)
      if (s[0] eq 3) then cube = total(cube,3) / s[3]
      
      ;; filters
      BB_filter = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS1 FILT NAME'),2)
      lyot_stop = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS1 OPTI1 NAME'),2)
      DB_filter = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS1 OPTI2 NAME'),2)
      
      ;; neutral densities
      ampli_cal = 0.0
      ampli_cpi = 0.0
      tmp = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS4 FILT1 NAME'),2)
      if (tmp ne 'OPEN') then ampli_cal = float(strmid(tmp,3,3))
      tmp = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS4 FILT2 NAME'),2)
      if (tmp ne 'OPEN') then nd_cpi = float(strmid(tmp,3,3)) else nd_cpi = 0.0
      
      print,' --> CAL ND = '+numformat(ampli_cal,deci=1)
      print,' --> CPI ND = '+numformat(nd_cpi,deci=1)
      print,' --> DIT    = '+numformat(dit,deci=2)
      print,' --> BB     = '+BB_filter
      print,' --> STOP   = '+lyot_stop
      print,' --> DB     = '+DB_filter
      
      ;; calibrate
      final_cube = fltarr(taillex,tailley,2)      
      
      img = cube
      img = img - pbkg
      img = img / flat

      ;; fix pixels
      apstep = 4
      maxap  = 10
      
      ;; left field
      img0_ini = img[x0_l:x1_l,y0_l:y1_l]
      img0 = maskinterp(img0_ini,bpm0,1,6,'csplinterp',gpix=10,gpoints=4,cdis=2)      
      if (resol eq 'MRS') then begin
         ;; cluster of bad pixels
         cluster_bpm = sph_ird_lss_cluster_bpm()
         cluster_bpm = abs(1-cluster_bpm)
         cluster_bpm0 = cluster_bpm[x0_l:x1_l,y0_l:y1_l]
         
         fixpix,img0,cluster_bpm0,img0_out,npix=12,/weight,/silent
         img0 = temporary(img0_out)
      endif
      img0 = sigma_filter(img0,7.,n_sigma=5,n_change=nchange)
      img0 = sigma_filter(img0,7.,N_sigma=3,radius=2,/iterate,n_change=n_change)
      final_cube[*,*,0] = img0 / dit
      
      ;; right field
      img1_ini = img[x0_r:x1_r,y0_r:y1_r]         
      img1 = maskinterp(img1_ini,bpm1,1,6,'csplinterp',gpix=10,gpoints=4,cdis=2)      
      if (resol eq 'MRS') then begin
         ;; cluster of bad pixels
         cluster_bpm = sph_ird_lss_cluster_bpm()
         cluster_bpm = abs(1-cluster_bpm)
         cluster_bpm0 = cluster_bpm[x0_r:x1_r,y0_r:y1_r]

         fixpix,img1,cluster_bpm0,img1_out,npix=12,/weight,/silent
         img1 = temporary(img1_out)
      endif
      img1 = sigma_filter(img1,7.,n_sigma=5,n_change=nchange)
      img1 = sigma_filter(img1,7.,N_sigma=3,radius=2,/iterate,n_change=n_change)
      final_cube[*,*,1] = img1 / dit
      
      writefits,aroot+'interm/'+psf_files[n]+'_clean'+suffix+'.fits.gz',final_cube,/compress
   endfor

   ;; ------------------------------
   ;; corono data
   ;;
   print,'+ corono data'
   
   ncube = n_elements(sci_files)
   for n=0,ncube-1 do begin
      print,' + cube '+numformat(n+1)+' / '+numformat(ncube)
      
      cube = readfits(aroot+'raw/'+sci_files[n]+'.fits',hdr)
      dit  = sxpar_eso(hdr,'HIERARCH ESO DET SEQ1 DIT')
      ndit = sxpar_eso(hdr,'HIERARCH ESO DET NDIT')
      rom  = sxpar_eso(hdr,'HIERARCH ESO DET READ CURNAME')

      if ~keyword_set(all_DIT) then begin
         ;; collpase cube
         s = size(cube)
         if (s[0] eq 3) then cube = total(cube,3) / s[3]
         ndit = 1
      endif

      final_cube = fltarr(taillex,tailley,2,ndit)         
      for d=0,ndit-1 do begin
         ;; filters
         BB_filter = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS1 FILT NAME'),2)
         lyot_stop = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS1 OPTI1 NAME'),2)
         DB_filter = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS1 OPTI2 NAME'),2)
         
         ;; neutral densities
         ampli_cal = 0.0
         ampli_cpi = 0.0
         tmp = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS4 FILT1 NAME'),2)
         if (tmp ne 'OPEN') then ampli_cal = float(strmid(tmp,3,3))
         tmp = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS4 FILT2 NAME'),2)
         if (tmp ne 'OPEN') then nd_cpi = float(strmid(tmp,3,3)) else nd_cpi = 0.0
         
         print,' --> CAL ND = '+numformat(ampli_cal,deci=1)
         print,' --> CPI ND = '+numformat(nd_cpi,deci=1)
         print,' --> DIT    = '+numformat(dit,deci=2)
         print,' --> BB     = '+BB_filter
         print,' --> STOP   = '+lyot_stop
         print,' --> DB     = '+DB_filter
         print
         
         ;; calibrate
         img = cube[*,*,d]
         img = img - bkg
         img = img / flat
         
         ;; fix pixels
         apstep = 4
         maxap  = 10
         
         ;; left field
         img0_ini = img[x0_l:x1_l,y0_l:y1_l]
         img0 = maskinterp(img0_ini,bpm0,1,6,'csplinterp',gpix=10,gpoints=4,cdis=2)      
         if (resol eq 'MRS') then begin
            ;; cluster of bad pixels
            cluster_bpm = sph_ird_lss_cluster_bpm()
            cluster_bpm = abs(1-cluster_bpm)
            cluster_bpm0 = cluster_bpm[x0_l:x1_l,y0_l:y1_l]

            fixpix,img0,cluster_bpm0,img0_out,npix=12,/weight,/silent
            img0 = temporary(img0_out)
         endif
         img0 = sigma_filter(img0,7.,n_sigma=5,n_change=nchange)
         img0 = sigma_filter(img0,7.,N_sigma=3,radius=2,/iterate,n_change=n_change)
         final_cube[*,*,0,d] = img0 / dit
         
         
         ;; right field
         img1_ini = img[x0_r:x1_r,y0_r:y1_r]         
         img1 = maskinterp(img1_ini,bpm0,1,6,'csplinterp',gpix=10,gpoints=4,cdis=2)      
         if (resol eq 'MRS') then begin
            ;; cluster of bad pixels
            cluster_bpm = sph_ird_lss_cluster_bpm()
            cluster_bpm = abs(1-cluster_bpm)
            cluster_bpm0 = cluster_bpm[x0_r:x1_r,y0_r:y1_r]

            fixpix,img1,cluster_bpm0,img1_out,npix=12,/weight,/silent
            img1 = temporary(img1_out)
         endif
         img1 = sigma_filter(img1,7.,n_sigma=5,n_change=nchange)
         img1 = sigma_filter(img1,7.,N_sigma=3,radius=2,/iterate,n_change=n_change)
         final_cube[*,*,1,d] = img1 / dit
      endfor

      writefits,aroot+'interm/'+sci_files[n]+'_clean'+suffix+'.fits.gz',final_cube,/compress
   endfor

   ;; ------------------------------
   ;; spectro calib PSF
   ;;
   print,'+ spectro calib PSF data'

   ncube = (cal_files[0] eq '') ? 0 : n_elements(cal_files)
   for n=0,ncube-1 do begin
      print,' + cube '+numformat(n+1)+' / '+numformat(ncube)

      cube = readfits(aroot+'raw/'+cal_files[n]+'.fits',hdr)
      dit  = sxpar_eso(hdr,'HIERARCH ESO DET SEQ1 DIT')
      ndit = sxpar_eso(hdr,'HIERARCH ESO DET NDIT')
      rom  = sxpar_eso(hdr,'HIERARCH ESO DET READ CURNAME')      

      ;; collpase cube if needed
      s = size(cube)
      if (s[0] eq 3) then cube = total(cube,3) / s[3]
      
      ;; filters
      BB_filter = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS1 FILT NAME'),2)
      lyot_stop = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS1 OPTI1 NAME'),2)
      DB_filter = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS1 OPTI2 NAME'),2)
      
      ;; neutral densities
      ampli_cal = 0.0
      ampli_cpi = 0.0
      tmp = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS4 FILT1 NAME'),2)
      if (tmp ne 'OPEN') then ampli_cal = float(strmid(tmp,3,3))
      tmp = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS4 FILT2 NAME'),2)
      if (tmp ne 'OPEN') then nd_cpi = float(strmid(tmp,3,3)) else nd_cpi = 0.0
      
      print,' --> CAL ND = '+numformat(ampli_cal,deci=1)
      print,' --> CPI ND = '+numformat(nd_cpi,deci=1)
      print,' --> DIT    = '+numformat(dit,deci=2)
      print,' --> BB     = '+BB_filter
      print,' --> STOP   = '+lyot_stop
      print,' --> DB     = '+DB_filter
      
      ;; calibrate
      final_cube = fltarr(taillex,tailley,2)
      
      img = cube
      img = img - cbkg
      img = img / flat

      ;; fix pixels
      apstep = 4
      maxap  = 10
      
      ;; left field
      img0_ini = img[x0_l:x1_l,y0_l:y1_l]
      img0 = maskinterp(img0_ini,bpm0,1,6,'csplinterp',gpix=10,gpoints=4,cdis=2)      
      if (resol eq 'MRS') then begin
         ;; cluster of bad pixels
         cluster_bpm = sph_ird_lss_cluster_bpm()
         cluster_bpm = abs(1-cluster_bpm)
         cluster_bpm0 = cluster_bpm[x0_l:x1_l,y0_l:y1_l]
         
         fixpix,img0,cluster_bpm0,img0_out,npix=12,/weight,/silent
         img0 = temporary(img0_out)
      endif
      img0 = sigma_filter(img0,7.,n_sigma=5,n_change=nchange)
      img0 = sigma_filter(img0,7.,N_sigma=3,radius=2,/iterate,n_change=n_change)
      final_cube[*,*,0] = img0 / dit
      
      ;; right field
      img1_ini = img[x0_r:x1_r,y0_r:y1_r]         
      img1 = maskinterp(img1_ini,bpm1,1,6,'csplinterp',gpix=10,gpoints=4,cdis=2)      
      if (resol eq 'MRS') then begin
         ;; cluster of bad pixels
         cluster_bpm = sph_ird_lss_cluster_bpm()
         cluster_bpm = abs(1-cluster_bpm)
         cluster_bpm0 = cluster_bpm[x0_r:x1_r,y0_r:y1_r]

         fixpix,img1,cluster_bpm0,img1_out,npix=12,/weight,/silent
         img1 = temporary(img1_out)
      endif
      img1 = sigma_filter(img1,7.,n_sigma=5,n_change=nchange)
      img1 = sigma_filter(img1,7.,N_sigma=3,radius=2,/iterate,n_change=n_change)
      final_cube[*,*,1] = img1 / dit
      
      writefits,aroot+'interm/'+cal_files[n]+'_clean'+suffix+'.fits.gz',final_cube,/compress
   endfor
endif

if keyword_set(do_combine) then begin
   ;;------------------------------
   ;; reference PSF
   ;;
   print,'Combining reference PSFs'
   
   ;; wavelength calibration
   wave = readfits(aroot+'calib/wave_cal.fits')
   wave0 = wave[x0_l:x1_l,y0_l:y1_l]
   wave1 = wave[x0_r:x1_r,y0_r:y1_r]   
   
   ;; read data
   ncube = n_elements(psf_files)
   cube = fltarr(taillex,tailley,2,ncube)
   for n=0,ncube-1 do begin
      data = readfits(aroot+'interm/'+psf_files[n]+'_clean'+suffix+'.fits.gz')
      cube[*,*,*,n] = data
   endfor   

   ;; combine data
   if ((size(cube))[0] eq 4) then data = median(cube,dim=4) $
   else data = cube
   
   pix  = findgen(taillex)
   ext  = 10
   
   ;; left field
   img0 = data[*,*,0] - median(data[*,*,0])
   tot0 = total(img0,2)
   fit0 = mpfitpeak(pix,tot0,a0,/gauss)
   
   img0 = translate(img0,(taillex-1)/2.-a0[1],0.0)
   lam0_1D = total(wave0[fix(a0[1])-ext:fix(a0[1])+ext,*],1)/(2*ext+1)
   lam0_2D = replicate(1,taillex) # lam0_1D
   
   ;; right field
   img1 = data[*,*,1] - median(data[*,*,1])
   tot1 = total(img1,2)
   fit1 = mpfitpeak(pix,tot1,a1,/gauss)
   
   img1 = translate(img1,(taillex-1)/2.-a1[1],0.0)
   lam1_1D = total(wave1[fix(a1[1])-ext:fix(a1[1])+ext,*],1)/(2*ext+1)
   lam1_2D = replicate(1,taillex) # lam1_1D

   ;; combined data
   psf = fltarr(taillex,tailley,2)
   psf[*,*,0] = img0
   psf[*,*,1] = img1

   lam = fltarr(tailley,2)
   lam[*,0] = lam0_1D
   lam[*,1] = lam1_1D

   lam_2D = fltarr(taillex,tailley,2)
   lam_2D[*,*,0] = lam0_2D
   lam_2D[*,*,1] = lam1_2D

   ;;------------------------------
   ;; chromatic effects
   ;;   
   if keyword_set(corr_chroma) then begin
      print,' * measure chromaticity'
      
      nlambda = n_elements(lam[*,0])
      pix = findgen(nlambda)
      
      center = fltarr(nlambda,2)
      chroma = fltarr(nlambda,2)
      for f=0,1 do begin
         ;tmp = sigma_filter(psf[*,*,f],5,n_sigma=3,/iter)

         for l=0,nlambda-1 do begin
            ;; s = reform(tmp[*,l])
            
            s = reform(psf[taillex/2-75:taillex/2+75,l,f])
            fit = mpfitpeak(findgen(151),s,a,/gauss)
            center[l,f] = a[1]+taillex/2-75
         endfor
      endfor

      if (resol eq 'LRS') then begin
         ii  = where(((1000 le lam0_1D) and (lam0_1D le 1300)) or $
                     ((1450 le lam0_1D) and (lam0_1D le 2300)))
         fit = poly_fit(pix[ii],center[ii,0],3)
         chroma[*,0] = fit[0]+pix*fit[1]+pix^2*fit[2]+pix^3*fit[3];+pix^4*fit[4]
         
         ii  = where(((1000 le lam1_1D) and (lam1_1D le 1300)) or $
                     ((1650 le lam1_1D) and (lam1_1D le 2300)))
         fit = poly_fit(pix[ii],center[ii,1],3)
         chroma[*,1] = fit[0]+pix*fit[1]+pix^2*fit[2]+pix^3*fit[3];+pix^4*fit[4]
      endif else if (resol eq 'MRS') then begin
         ;; ii  = where((1000 le lam0_1D) and (lam0_1D le 1800))
         ii  = where(((1000 le lam0_1D) and (lam0_1D le 1300)) or $
                     ((1450 le lam0_1D) and (lam0_1D le 1650)))
         fit = poly_fit(pix[ii],center[ii,0],1)
         chroma[*,0] = fit[0]+pix*fit[1];+pix^2*fit[2]
         
         ;; ii  = where((1000 le lam1_1D) and (lam1_1D le 1800))
         ii  = where(((1000 le lam1_1D) and (lam1_1D le 1300)) or $
                     ((1650 le lam1_1D) and (lam1_1D le 1800)))
         fit = poly_fit(pix[ii],center[ii,1],1)
         chroma[*,1] = fit[0]+pix*fit[1];+pix^2*fit[2]
      endif
         
      ;; mplot,/open,eps=1,x=6.4,y=4.8,file='~/Desktop/ADC_test.eps'
      loadct,39,/silent
      wopen,0,xs=1200,ys=500
      !p.multi = [0,2,1]
      
      plot,lam0_1D,center,/nodata,xs=1,psym=1,/yno,ys=1,yr=taillex/2.-1+[-3,3], $
           xtitle='Wavelength [nm]',ytitle='Center position [pix]',title='Field 0'
      plots,!x.crange,(taillex-1)/2.,linestyle=2
      
      oplot,lam0_1D,center[*,0],color=50
      oplot,lam0_1D,chroma[*,0],color=100,linestyle=1

      plot,lam0_1D,center,/nodata,xs=1,psym=1,/yno,ys=1,yr=taillex/2.-1+[-3,3], $
           xtitle='Wavelength [nm]',ytitle='Center position [pix]',title='Field 1'
      plots,!x.crange,(taillex-1)/2.,linestyle=2

      oplot,lam1_1D,center[*,1],color=250
      oplot,lam1_1D,chroma[*,1],color=150,linestyle=1

      ;; print,transpose([[lam0_1D],[center[*,0]],[center[*,1]]])
      ;; mplot,/close,/disp

      chroma[*,0] = (taillex-1)/2. - chroma[*,0]
      chroma[*,1] = (taillex-1)/2. - chroma[*,1]

      !p.multi = 0
      hak
   endif
   
   ;; correct effect for PSF
   if keyword_set(corr_chroma) then begin
      print,' * correct chromaticity'
      
      ;; PSF
      for l=0,nlambda-1 do begin
         psf[*,l,0] = subpixel_shift_1d(psf[*,l,0],pshift=chroma[l,0])
         psf[*,l,1] = subpixel_shift_1d(psf[*,l,1],pshift=chroma[l,1])
      endfor      
   endif
   
   ;;------------------------------
   ;; spectro calib PSF
   ;;
   if (cal_files[0] ne '') then begin   
      print,'Combining telluric standard'
      
      ;; read data
      ncube = n_elements(cal_files)
      cube = fltarr(taillex,tailley,2,ncube)
      for n=0,ncube-1 do begin
         data = readfits(aroot+'interm/'+cal_files[n]+'_clean'+suffix+'.fits.gz')
         cube[*,*,*,n] = data
      endfor   

      ;; combine data
      if ((size(data))[0] eq 4) then data = median(cube,dim=4) $
      else data = cube
      
      pix  = findgen(taillex)
      ext  = 10
      
      ;; left field
      img0 = data[*,*,0]
      tot0 = total(img0,2)
      fit0 = mpfitpeak(pix,tot0,a0,/gauss)
      
      img0 = translate(img0,(taillex-1)/2.-a0[1],0.0)
      lam0_1D = total(wave0[fix(a0[1])-ext:fix(a0[1])+ext,*],1)/(2*ext+1)
      lam0_2D = replicate(1,taillex) # lam0_1D

      ;; right field
      img1 = data[*,*,1]
      tot1 = total(img1,2)
      fit1 = mpfitpeak(pix,tot1,a1,/gauss)
      
      img1 = translate(img1,(taillex-1)/2.-a1[1],0.0)
      lam1_1D = total(wave1[fix(a1[1])-ext:fix(a1[1])+ext,*],1)/(2*ext+1)
      lam1_2D = replicate(1,taillex) # lam1_1D

      ;; combined data
      cal = fltarr(taillex,tailley,2)
      cal[*,*,0] = img0
      cal[*,*,1] = img1

      ;; correct chromatic effect
      if keyword_set(corr_chroma) then begin
         print,' * correct chromaticity'
         for l=0,nlambda-1 do begin
            cal[*,l,0] = subpixel_shift_1d(cal[*,l,0],pshift=chroma[l,0])
            cal[*,l,1] = subpixel_shift_1d(cal[*,l,1],pshift=chroma[l,1])
         endfor
      endif
   endif
      
   ;; ------------------------------
   ;; corono data
   ;;
   print,'Combining coronagraphic data'

   ;; number of DITs
   ncube = n_elements(sci_files)
   ndit  = 0
   for n=0,ncube-1 do begin
      hdr = headfits(aroot+'interm/'+sci_files[n]+'_clean'+suffix+'.fits.gz')

      naxis = sxpar(hdr,'NAXIS')
      if (naxis eq 3) then ndit += 1 $
      else if (naxis eq 4) then ndit += sxpar(hdr,'NAXIS4')
   endfor      
   
   ;; read data
   ncube = n_elements(sci_files)
   cube = fltarr(taillex,tailley,2,ndit)
   idx = 0
   for n=0,ncube-1 do begin
      data = readfits(aroot+'interm/'+sci_files[n]+'_clean'+suffix+'.fits.gz',hdr)

      naxis = sxpar(hdr,'NAXIS')
      if (naxis eq 3) then inc = 1 $
      else if (naxis eq 4) then inc = sxpar(hdr,'NAXIS4')
      
      cube[*,*,*,idx:idx+inc-1] = data
      idx += inc
   endfor   
   
   ;;
   ;; center collapsed frame
   ;;
   maxdim = max([taillex,tailley])
   if ((size(cube))[0] eq 4) then coll = median(cube,dim=4)
   
   ;; correct chromatic effect for collapsed image
   if keyword_set(corr_chroma) then begin
      print,' * correct chromaticity'
      for l=0,nlambda-1 do begin
         coll[*,l,0] = subpixel_shift_1d(coll[*,l,0],pshift=chroma[l,0])
         coll[*,l,1] = subpixel_shift_1d(coll[*,l,1],pshift=chroma[l,1])
      endfor
   endif

   img0 = double(coll[*,*,0])
   img1 = double(coll[*,*,1])

   ;; filtering can improve result if remaining large scale structure
   ;; img0 = img0 - median(img0,15)
   ;; img1 = img1 - median(img1,15)
   
   ref0 = img0
   ref1 = img1         
   
   fftref0 = fft(ref0,-1,dim=1)
   fftref1 = fft(ref1,-1,dim=1)
   
   ;; center left field
   ref0 = img0
   ;; message,'!!!!!!!!!!!!! Warning !!!!!!!!!!!!!',/inform
   ;; ref0[197-10:197+10,*] = 0
   inv0 = reverse(ref0,1)
   
   fftref0 = fft(ref0,-1,dim=1)
   fftinv0 = fft(inv0,-1,dim=1)
   
   cps0 = (fftinv0*conj(fftref0)) / abs(fftinv0*conj(fftref0))
   res0 = abs(fft(cps0,1,dim=1))^2         
   tot0 = total(res0,2,/nan)
   tot0[50:*] = 0.0
   ext  = 10
   max0 = max(tot0,imax)
   fit0 = mpfitpeak(pix[imax-ext:imax+ext],tot0[imax-ext:imax+ext],a0,/gauss)
   sh0  = a0[1]/2

   tmp0 = fltarr(maxdim,maxdim)
   tmp0[0:taillex-1,0:tailley-1] = img0
   tmp0 = subpixel_shift(tmp0,xs=sh0,ys=0)
   img0 = tmp0[0:taillex-1,0:tailley-1]
   
   ;; center right field
   ref1 = img1
   ;; message,'!!!!!!!!!!!!! Warning !!!!!!!!!!!!!',/inform
   ;; ref1[197-10:197+10,*] = 0   
   inv1 = reverse(ref1,1)

   fftref1 = fft(ref1,-1,dim=1)
   fftinv1 = fft(inv1,-1,dim=1)

   cps1 = (fftinv1*conj(fftref1)) / abs(fftinv1*conj(fftref1))
   res1 = abs(fft(cps1,1,dim=1))^2
   tot1 = total(res1,2,/nan)
   tot1[taillex/2:*] = 0.0
   tot1[50:*] = 0.0
   ext  = 10
   max1 = max(tot1,imax)         
   fit1 = mpfitpeak(pix[imax-ext:imax+ext],tot1[imax-ext:imax+ext],a1,/gauss)
   sh1  = a1[1]/2

   tmp1 = fltarr(maxdim,maxdim)
   tmp1[0:taillex-1,0:tailley-1] = img1
   tmp1 = subpixel_shift(tmp1,xs=sh1,ys=0)
   img1 = tmp1[0:taillex-1,0:tailley-1]
   
   coll[*,*,0] = img0
   coll[*,*,1] = img1

   print,'Centering shift:',sh0,sh1
   
   if ~keyword_set(all_DIT) then begin
      coro = temporary(coll)
   endif else begin
      for d=0,ndit-1 do begin
         img0 = cube[*,*,0,d]
         img1 = cube[*,*,1,d]

         tmp0 = fltarr(maxdim,maxdim)
         tmp0[0:taillex-1,0:tailley-1] = img0
         if ~keyword_set(corr_chroma) then begin
            ;; shift only if not correcting chromaticity
            ;; otherwise, shift is done below
            tmp0 = subpixel_shift(tmp0,xs=sh0,ys=0)
         endif
         tmp0 = tmp0[0:taillex-1,0:tailley-1]

         tmp1 = fltarr(maxdim,maxdim)
         tmp1[0:taillex-1,0:tailley-1] = img1
         if ~keyword_set(corr_chroma) then begin
            ;; shift only if not correcting chromaticity
            ;; otherwise, shift is done below
            tmp1 = subpixel_shift(tmp1,xs=sh1,ys=0)
         endif
         tmp1 = tmp1[0:taillex-1,0:tailley-1]
         
         cube[*,*,0,d] = tmp0
         cube[*,*,1,d] = tmp1
      endfor
      coro = temporary(cube)

      ;; correct chromatic effect
      if keyword_set(corr_chroma) then begin
         print,' * correct chromaticity'
         for l=0,nlambda-1 do begin
            for d=0,ndit-1 do begin
               coro[*,l,0,d] = subpixel_shift_1d(coro[*,l,0,d],pshift=chroma[l,0]+sh0)
               coro[*,l,1,d] = subpixel_shift_1d(coro[*,l,1,d],pshift=chroma[l,1]+sh1)
            endfor
         endfor
      endif      
   endelse
   
   ;;------------------------------
   ;; final data
   ;;   
   
   ;; cut signal to usable area, choose field
   cut = 20
   if (resol eq 'LRS') then begin
      lmin = 920
      lmax = 2400
   endif else if (resol eq 'MRS') then begin
      lmin = 920
      lmax = 1900
   endif else message,'Resolution '+resol+' is undefined!'
   
   for ifield=0,1 do begin
      dim = (size(psf,/dim))[0]
      ll  = where((lmin le lam[*,ifield]) and (lam[*,ifield] le lmax),nll0)
      
      if (nll0 mod 2) then ll = ll[1:nll0-1]
      
      lam_1D_s = reverse(lam[ll,ifield])
      lam_2D_s = lam_2D[cut:dim-1-cut,ll,ifield]
      psf_s    = psf[cut:dim-1-cut,ll,ifield]
      coro_s   = reform(coro[cut:dim-1-cut,ll,ifield,*])
      
      writefits,aroot+'products/reference_psf'+suffix+'_F'+numformat(ifield)+'.fits.gz',psf_s,/compress
      writefits,aroot+'products/wavelength'+suffix+'_F'+numformat(ifield)+'.fits.gz',lam_1D_s,/compress
      writefits,aroot+'products/wavelength_2D'+suffix+'_F'+numformat(ifield)+'.fits.gz',lam_2D_s,/compress
      
      if ~keyword_set(all_DIT) then begin
         writefits,aroot+'products/combined_corono'+suffix+'_F'+numformat(ifield)+'.fits.gz',coro_s,/compress
      endif else begin
         coro_s_combined = total(coro_s,3) / (size(coro_s))[3]
         writefits,aroot+'products/combined_corono'+suffix+'_F'+numformat(ifield)+'.fits.gz',coro_s_combined,/compress
         writefits,aroot+'products/individual_corono'+suffix+'_F'+numformat(ifield)+'.fits.gz',coro_s,/compress
      endelse
      
      if (cal_files[0] ne '') then begin
         cal_s = cal[cut:dim-1-cut,ll,ifield]
         writefits,aroot+'products/spectro_cal'+suffix+'_F'+numformat(ifield)+'.fits.gz',cal_s,/compress
      endif

      ;; cflat = readfits(aroot+'interm/coro_flat.fits.gz')
      ;; coro_flat = cflat[cut:dim-1-cut,ll,ifield]
      ;; coro_flat_lin = reverse(total(coro_flat,1) / (size(coro_flat))[1])
      ;; writefits,aroot+'products/coro_flat'+suffix+'_F'+numformat(ifield)+'.fits.gz',coro_flat_lin,/compress
      
      dim = (size(psf_s,/dim))
      ds9
      !v->im,coro_s[*,*,*,0]
      !v->box,(dim[0]-1)/2.,(dim[1]-1)/2.,30,dim[1],0
      ;; !v->im,psf_s
      ;; !v->line,(dim[0]-1)/2.,0,(dim[0]-1)/2.,dim[1]
   endfor
endif

if keyword_set(do_specsub) then begin
   if file_test(aroot+'products/wavelength_corrected'+suffix+'_'+field+'.fits.gz') then begin
      lambda = readfits(aroot+'products/wavelength_corrected'+suffix+'_'+field+'.fits.gz')
   endif else begin
      lambda = readfits(aroot+'products/wavelength'+suffix+'_'+field+'.fits.gz')
   endelse
   if keyword_set(all_DIT) then begin
      coro = readfits(aroot+'products/individual_corono'+suffix+'_'+field+'.fits.gz')
      ndit = (size(coro,/dim))[2]
      nsig = fltarr((size(coro,/dim))[1],(size(coro,/dim))[0],(size(coro,/dim))[2])
      for d=0,ndit-1 do nsig[*,*,d] = rotate(coro[*,*,d],1)
      sig  = temporary(nsig)
   endif else begin
      coro = readfits(aroot+'products/combined_corono'+suffix+'_'+field+'.fits.gz')
      ndit = 1
      sig  = rotate(coro,1)
   endelse

   psf = readfits(aroot+'products/reference_psf'+suffix+'_'+field+'.fits.gz')
   psf = rotate(psf,1)
   
   nlambda = n_elements(lambda)
   nsep    = (size(coro,/dim))[0]
   center  = (nsep-1)/2.0
      
   w     = nlambda
   h     = nsep
   bigh  = h*2
   scale = float(bigh)/ h

   lsurD = lambda*1d-9/8*180/!pi*3600/pixel

   ;; fake planet
   if keyword_set(pla_fake) then begin
      ;; neutral density
      hdr_psf   = headfits(aroot+'raw/'+psf_files[0]+'.fits')
      ampli_psf = 0.0
      nd_psf = strtrim(sxpar_eso(hdr_psf,'HIERARCH ESO INS4 FILT2 NAME'),2)
      ampli_psf = sph_cpi_nd_transmission(nd_psf,lambda)

      a = (0-pla_col_fake) / (1600-1200)
      b = 0-1600*a
      correction = 10^(-(lambda*a + b)/2.5)

      npsf = psf * 10^(replicate(1,taillex) ## (ampli_psf*correction))
      bad = where(finite(npsf) ne 1,nbad)
      if (nbad ne 0) then npsf[bad] = 0.0

      pla = npsf * 10^(-pla_cnt_fake/2.5)
      pla = translate(pla,0,pla_sep/pixel)

      sig_nopla = sig
      for d=0,ndit-1 do sig[*,*,d] = sig[*,*,d] + pla

      writefits,aroot+'products/fake_planet_'+method+suffix+'_'+field+'.fits.gz',pla
   endif else begin
      npsf = psf
   endelse
   bad = where(finite(sig) eq 0,nbad)
   if (nbad gt 0) then sig[bad] = 0.0

   sig_final = fltarr(w,h,ndit)
   t_start = systime(/sec)
   for d=0,ndit-1 do begin
      t_left = ((systime(/sec)-t_start)/d*(ndit-d)) / 60D
      print,d+1,ndit,t_left,format='(" * DIT ",I0,"/",I0," - ",D0.1," min left")'

      ;; radius of the coronagraphic mask in as
      ;; currently only the 0.20" mask is available
      lyot_size = 0.20
            
      ;; rescaling and oversampling      
      sig_over = sph_ird_lss_rescale(sig[*,*,d],lambda, $
                                     oversample=scale,lambda_ref=lambda[0], $
                                     pixel=pixel,coro_rad=lyot_size, $
                                     planet_sep=pla_sep,planet_size=pla_size, $
                                     coro_mask=use_mask,use_mask=use_mask_over, $
                                     planet_mask=pla_mask_over)

      if (method eq 'simplesub') then begin
         synthetic = sig[*,*,d]
         synthetic[*,0:h/2-1] = reverse(synthetic[*,h/2:*],2)
      endif else begin         
         ;; remove speckles
         if (method eq 'vigan2008') then begin
            synthetic_over = sph_ird_lss_synthetic('synthetic_vigan2008',sig_over,lambda, $
                                                   use_mask_over,pla_mask_over,fill=1, $
                                                   oversample=scale,lambda_ref=lambda[0], $
                                                   pixel=pixel,lyot_rad=lyot_size, $
                                                   limit=80,tr_lim=0.75)
         endif else if (method eq 'symmetry') then begin
            synthetic_over = sph_ird_lss_synthetic('synthetic_symmetry',sig_over,lambda, $
                                                   use_mask_over,pla_mask_over,fill=0)
         endif else if (method eq 'locs') then begin
            synthetic_over = sph_ird_lss_synthetic('synthetic_locs',sig_over,lambda, $
                                                   use_mask_over,pla_mask_over,fill=1, $
                                                   pixel=pixel,lambda_ref=lambda[0], $
                                                   planet_sep=pla_sep,planet_size=pla_size, $
                                                   locs_crit=locs_crit,normalize=0, $
                                                   oversample=scale,silent=0,spatial=0)
         endif
         final_over = (sig_over-synthetic_over)*use_mask_over
         
         ;; rescaling
         synthetic = sph_ird_lss_rescale(synthetic_over,lambda,/reverse, $
                                         oversample=scale,lambda_ref=lambda[0])
      endelse
   
      ;; final signal
      sig_final[*,*,d] = (sig[*,*,d] - synthetic)*use_mask

      ;; save in temporary file
      writefits,aroot+'products/temporary_analyzed_corono_'+method+ $
                suffix+'_'+field+'.fits.gz',sig_final[*,*,0:d]
      print
   endfor

   ;; save
   final = fltarr(w,h,2+ndit)
   if keyword_set(all_DIT) then begin
      final[*,*,0] = psf
      final[*,*,1] = total(sig,3)
      final[*,*,2:*] = sig_final
      writefits,aroot+'products/analyzed_corono_'+method+suffix+postfix+'_'+field+'_individual.fits.gz',final
   endif else begin
      final[*,*,0] = psf
      final[*,*,1] = sig
      final[*,*,2:*] = sig_final
      writefits,aroot+'products/analyzed_corono_'+method+suffix+postfix+'_'+field+'.fits.gz',final
   endelse

   ;; remove temporary file
   file_delete,aroot+'products/temporary_analyzed_corono_'+method+ $
               suffix+'_'+field+'.fits.gz'
   
   ;; show final product   
   !v->im,sig_final[*,*,0]
   !v->box,nlambda/2,(h-1)/2.+pla_sep/pixel,nlambda,10,0
endif

if keyword_set(do_extract) then begin
   if keyword_set(all_DIT) then begin
      final = readfits(aroot+'products/analyzed_corono_'+method+suffix+postfix+'_'+field+'_individual.fits.gz')
      dim   = size(final,/dim)
      tmp   = fltarr(dim[0],dim[1],3)
      tmp[*,*,0] = final[*,*,0]
      tmp[*,*,1] = final[*,*,1]
      tmp[*,*,2] = median(final[*,*,2:*],dim=3)
      final = tmp
   endif else begin
      final = readfits(aroot+'products/analyzed_corono_'+method+suffix+postfix+'_'+field+'.fits.gz')
   endelse
   if file_test(aroot+'products/wavelength_corrected'+suffix+'_'+field+'.fits.gz') then begin
      lambda = readfits(aroot+'products/wavelength_corrected'+suffix+'_'+field+'.fits.gz')
   endif else begin
      lambda = readfits(aroot+'products/wavelength'+suffix+'_'+field+'.fits.gz')
   endelse

   if keyword_set(pla_fake) then begin
      fake_pla = readfits(aroot+'products/fake_planet_'+method+suffix+'_'+field+'.fits.gz')
   endif
   
   psf       = final[*,*,0]
   sig       = final[*,*,1]
   sig_final = final[*,*,2]

   nlambda = n_elements(lambda)
   nsep    = (size(sig,/dim))[0]
   center  = (nsep-1)/2.0
   
   ds9
   !v->im,[[[sig_final]],[[sig]]]
   
   aper_size = 1.0
   aper_pla = sph_ird_lss_aperture(sig_final,lambda,pla_sep,aper_diam=aper_size,pixel=pixel, $
                                   lambda_ref=lambda[0],aper_rescaled=aper_pla_r)   
   aper_pla_sym = sph_ird_lss_aperture(sig_final,lambda,-pla_sep,aper_diam=aper_size, $
                                       pixel=pixel,lambda_ref=lambda[0], $
                                       aper_rescaled=aper_pla_r)
   aper_psf = sph_ird_lss_aperture(sig_final,lambda,0.0,aper_diam=aper_size,pixel=pixel)

   lin_psf = total(aper_psf*psf,2)
   lin_hal = total(aper_pla_sym*sig,2)
   lin_sig = total(aper_pla*sig,2)
   lin_ext = total(aper_pla*sig_final,2)
   lin_sym = total(aper_pla_sym*sig_final,2)
   if keyword_set(pla_fake) then lin_fak = total(aper_pla*fake_pla,2)
   
   coronoisered = lss_noise(sig_final,lambda,pixel,nlsd=5)
   lin_noise_sym = total(aper_pla_sym*coronoisered,2)

   ;;
   ;; correct for ND
   ;;
   hdr_psf   = headfits(aroot+'raw/'+psf_files[0]+'.fits')
   ampli_psf = 0.0
   nd_psf = strtrim(sxpar_eso(hdr_psf,'HIERARCH ESO INS4 FILT2 NAME'),2)
   ampli_psf = sph_cpi_nd_transmission(nd_psf,lambda)

   hdr_sig   = headfits(aroot+'raw/'+sci_files[0]+'.fits')
   ampli_sig = 0.0
   nd_sig = strtrim(sxpar_eso(hdr_sig,'HIERARCH ESO INS4 FILT2 NAME'),2)
   ampli_sig = sph_cpi_nd_transmission(nd_sig,lambda)

   lin_psf = lin_psf * 10^ampli_psf
   lin_sig = lin_sig * 10^ampli_sig
   lin_hal = lin_hal * 10^ampli_sig
   lin_ext = lin_ext * 10^ampli_sig
   lin_sym = lin_sym * 10^ampli_sig
   lin_noise_sym = lin_noise_sym * 10^ampli_sig
   if keyword_set(pla_fake) then lin_fak = lin_fak * 10^ampli_sig
   
   ;;
   ;; spectro calibrator
   ;;
   if (cal_files[0] ne '') then begin
      cal = readfits(aroot+'products/spectro_cal'+suffix+'_'+field+'.fits.gz')
      cal = rotate(cal,1)

      hdr_cal   = headfits(aroot+'raw/'+cal_files[0]+'.fits')
      ampli_cal = 0.0
      nd_cal = strtrim(sxpar_eso(hdr_cal,'HIERARCH ESO INS4 FILT2 NAME'),2)
      ampli_cal = sph_cpi_nd_transmission(nd_cal,lambda)
      
      lin_cal = total(aper_psf*cal,2)
      lin_cal = lin_cal * 10^ampli_cal
   endif else lin_cal = replicate(1D,nlambda)

   ;;
   ;; save
   ;;
   if keyword_set(all_DIT) then suffix = suffix+'_individual'
   if keyword_set(pla_fake) then begin
      openw,lun,aroot+'spectra_'+method+'_'+field+suffix+postfix+'.dat',/get_lun
      printf,lun,'       lambda          PSF          fake_pla_out     fake_pla_in'
      printf,lun,'----------------------------------------------------------------'
      printf,lun,transpose([[lambda],[lin_psf],[lin_ext],[lin_fak]])
      free_lun,lun      
   endif else begin
      openw,lun,aroot+'spectra_'+method+'_'+field+suffix+'.dat',/get_lun
      printf,lun,'       lambda         PSF       planet     planet_sym'
      printf,lun,'-----------------------------------------------------'
      printf,lun,transpose([[lambda],[lin_psf],[lin_ext],[lin_noise_sym]])
      free_lun,lun
   endelse
   
   ;;
   ;; plot
   ;;
   if (resol eq 'LRS') then begin
      xmin = 900
      xmax = 2300
   endif else if (resol eq 'MRS') then begin
      xmin = 900
      xmax = 1900
   endif
   
   ;;
   ;; extracted spectra
   ;;
   mplot,/open,eps=eps,xs=6.4,ys=5,zoom=1.2,thick=4,win=0, $
         file=aroot+'spectra_comparison.eps'
   loadct,39,/silent

   ymin = 5d-8
   ymax = 2

   norm = max(lin_psf,/nan)
   plot,lambda,[0],/nodata,ylog=1,ys=1,yr=[ymin,ymax],xs=1,xr=[xmin,xmax], $
        xtitle='Wavelength [nm]',ytitle='Signal [normalized]',xmargin=[8,2], $
        ymargin=[3.5,1]
   oplot,lambda,lin_psf/norm,linestyle=1
   oplot,lambda,lin_hal/norm,color=150,linestyle=2
   oplot,lambda,lin_sig/norm,color=250
   oplot,lambda,lin_ext/norm,color=50
   oplot,lambda,lin_noise_sym/norm,color=200

   if keyword_set(pla_fake) then oplot,lambda,lin_fak/norm,color=80

   al_legend,['Star (peak)','Halo+planet','Halo (symmetric)','Planet (extracted)','Noise (symmetric)'], $
             color=[0,250,150,50,200],linestyle=[1,0,2,0,0],/bottom,/right,/clear
   
   loadct,0,/silent
   mplot,/close,/disp

   ;;
   ;; contrast
   ;;
   mplot,/open,eps=eps,xs=6.4,ys=5,zoom=1.2,thick=4,win=1, $
         file=aroot+'planet_spectrum_contrast.eps'
   loadct,39,/silent

   planet = lin_ext/lin_psf

   ymin = 1e-6
   ymax = 1e-3
   plot,lambda,[0],/nodata,ylog=1,ys=9,yr=[ymin,ymax],xs=1,xr=[xmin,xmax], $
        xtitle='Wavelength [nm]',ytitle='Contrast',xmargin=[8,8], $
        ymargin=[3.5,1]
   !y.type = 0
   axis,/yaxis,ys=1,yr=-2.5*alog10(10^!y.crange),ytitle='Contrast [mag]'
   !y.type = 1
   
   oplot,lambda,planet,color=250

   if keyword_set(pla_fake) then oplot,lambda,lin_fak/norm,color=200
   
   loadct,0,/silent
   mplot,/close,/disp
endif

fin:

end

