;+
; NAME:
;
;  data_reduction_ifs
;
; PURPOSE:
;
;  Performs the final steps of the SPHERE/IFS data reduction
;
; CALLING SEQUENCE:
;
;  From the IDL command prompt:
;
;    IDL> .run data_reduction_ifs
;
; DESCRIPTION:
;
;  This IDL program is designed to perform the final steps of the
;  SPHERE/IFS data reduction, mainly clean the few remaining bad
;  pixels, recalibrate the wavelength using a star center (C) frame
;  and combine all the data inside a single data cube.
;
;  For the procedure to work, you have to enter some basic information
;  concerning your reduction:
;
;    ROOT - path to the data reduction directory. It should be the
;           same at the one you have used in the data_reduction_ifs.sh
;           program
;
;    WAVE_FILE - name of the raw wavelength calibration file
;
;    CENT_FILE - name of the raw star center (C) calibration file. If
;                you have several in your observing sequence, you
;                should normally be able to use any of them
;
;    FLUX_FILE - name of the raw off-axis PSF (F) files. The variable
;                is a vector that can contain several file names
;
;    CORO_FILE - name of the science (O or C) files. The variable is a
;                vector that can contain several file names
;
;  Once the path and file names have been specified, you can execute
;  the three calibration steps sequentially. They are controled by the
;  keywords at the beginning: DO_CLEAN, DO_CENTER and DO_COMBINE.
;
;  The DO_CLEAN step uses sigma-clipping to identify the remaining
;  bad pixels, and then the maskinterp() routine to correct for
;  them. The sigma-clipping level can be controled by setting the
;  value of the NSIGMA variable below the control keywords. The
;  cleaned cubes are saved into temporary files with the *_clean
;  suffix in the products/ subdirectory. Note that you can skip the
;  cleaning by setting NSIGMA = 0.0. It will not clean the data but it
;  will still copy the data into new files that are required for the
;  subsequent steps.
;
;  The DO_CENTER step performs two distinct actions. The first action
;  is to determine the position of the star behind the coronagraph
;  using the dedicated calibration frames (C). The satellite spots are
;  fitted using a 2D gaussian. The output of this step is (1) the star
;  center in each spectral channel and (2) the wavelength dependent
;  scaling factor that is used to recalibrate the wavelength (see
;  below) and can be used in the data analysis to perform accurate
;  spatial rescaling for spectral differential imaging (SDI). The
;  second action is the recalibration of the wavelength using the
;  wavelength calibration frame processed as a science frame in the
;  data_reduction_IFS.sh script and the scaling factor described
;  above. This procedure is detailed in Vigan et al. (2015, MNRAS,
;  454, 129). The outputs of the DO_CENTER step are saved as FITS
;  files into the products/ subdirectory:
;
;    - data_centers.fits: vector giving the (x,y) centers for all
;                         spectral channels (2x39 elements)
;
;    - data_wavelength.fits: vector giving the recalibrated wavelength
;                            for all spectral channels (39 elements)
;
;    - data_scaling.fits: vector giving the rescaling factor with
;                         respect to the first spectral channel
;                         (39 elements)
;
;  Finally, the DO_COMBINE step will combine all the science data
;  cubes into a single master cube. The frames are centered using
;  either a FFT-based or interpolation-based subpixel shift
;  procedures. The procedure to be used can be specified with the
;  FFT_SHIFT variable variable below the control keywords. Note that
;  FFT-based routine is NOT provided with the pipeline. Replace it by
;  your own or use the provided interpolation-based routine. This step
;  has three outputs that are saved as FITS files in the products/
;  subdirectory:
;
;     - data_cube_psf.fits: data cube containing all the off-axis PSFs 
;
;     - data_cube_coro.fits: data cube containing all the science data
;
;     - data_cube_info.fits: binary FITS table containing the basic
;                           information for all science frames (see
;                           details below)
;
;  IMPORTANT: during this combination step, each frame is normalized
;  to a DIT of 1 second, and the effect of any neutral density filter
;  that was in the optical path during the observation (e.g. for the
;  off-axis PSF) is compensated for.
;
;  The data_cube_info.fits contains a binary FITS table with as many
;  entries as the number of science frames (DITs). Each entry contains
;  the basic information about each frame: time (at start/middle/end
;  of the exposure), hour angle (start/middle/end), parallactic angle
;  (start/middle/end), pupil orientation, DIT value. The content of
;  the table can be read and used as follow:
;
;    IDL> info = mrdfits('data_info.fits',1)   ;; read data
;    IDL> plot,info[*].time,info[*].pa         ;; plot pa values
;    IDL> pa = info[*].pa                      ;; extract pa vector
;
;  The pupil offset value corresponds to a fixed angular offset when
;  the instrument is used in pupil-tracking (for angular differential
;  imaging). It is fixed for all observations, although there is an
;  entry in the table for each science frame. The value is hard-coded
;  in the reduction pipeline and has been calibrated during various
;  commissionings. However, it does not take into account any fine
;  true North orientation correction that would be required for
;  accurate astrometry. Below is a simple IDL example that shows how
;  to use the table values to derotate your frames at the end of an
;  ADI analysis for all cubes and wavelengths:
;
;    for c=0,ncubes-1 do begin
;       for l=0,nlambda-1 do begin
;          frame = cube[*,*,l,c]
;          cube[*,*,l,c] = rot(frame,-parang[c]-pupoff,1.0,cx,cy)
;       endfor
;    endfor
;
;  Finally, the DO_ERASE keyword offers the possibility to erase the
;  temporary files that were created by the DO_CLEAN at the
;  beginning. Confirmation at the prompt will be asked before erasing
;  the files.
;
; REQUIRED PROGRAMS AND LIBRARIES:
;
;  The following IDL generic libraries are necessary:
;
;    * astronomy user's library: http://idlastro.gsfc.nasa.gov/
;    * MPFIT: https://www.physics.wisc.edu/~craigm/idl/fitting.html
;    * maskinterp: http://astro.vigan.fr/tools/maskinterp-1.2.tar.gz
;
;  The following SPHERE-specific library is also necessary:
;
;    * SPHERE transmission: http://astro.vigan.fr/tools.html
;
; MODIFICATION HISTORY:
;
;  arthur.vigan - 08/2015 - public version of personnal tools
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
;   because this script is the core of our public SPHERE/IFS reduction
;   pipeline, we would be grateful if you could cite the following
;   paper in any publication making use of it:
;
;     Vigan et al., 2015, MNRAS, 454, 129
;
;   We thank you for your effort, and hope that this tool will
;   contribute to your scientific work and discoveries. Please feel
;   free to report any bug or possible improvement to the author
;
;-

do_clean   = 1
do_center  = 1
do_combine = 1
do_erase   = 1

nsigma    = 3.5  ;; sigma-clipping threshold for additional cleaning
fft_shift = 0    ;; 0: interpolation-based shift; 1: FFT-based shift

;;
;; Target data configuration
;;
root = '/data/TargetDirectory/'

wave_file = 'SPHER.2015-01-01T00:00:00.000'
cent_file = 'SPHER.2015-01-01T00:00:01.000'
flux_file = ['SPHER.2015-01-01T00:00:02.000']
coro_file = ['SPHER.2015-01-01T00:00:03.000','...']

;; #####################################################################################

;; fixed parameters
nlambda = 39       ;; number of wavelength channels (no do modify!)
ntaille = 280      ;; final size of the data cubes (do not modify!)
pixel   = 7.3      ;; pixel size [mas]

data_file = [cent_file,flux_file,coro_file]

;;
;; additional cleaning of the (x,y,lambda) science cubes
;;
if keyword_set(do_clean) then begin
   print,'Final cleaning of the data'
   nfiles = n_elements(data_file)
   for f=0,nfiles-1 do begin
      ;; find all the products associated to the current raw science file
      dits  = file_search(root+'/products/'+data_file[f]+'*[0-9].fits',count=ndits)
      infos = file_search(root+'/products/'+data_file[f]+'*info.fits',count=ninfos)

      ;; check for potential problems
      if (ninfos gt 1) then $
         message,'More than one info file for '+data_file[f]

      ;; read info file
      info = mrdfits(infos[0],1,/silent)
      if (n_elements(info) ne ndits) then $
         message,'Number of info does not correspond to the number ' + $
                 'of DITs for '+data_file[f]+'!'

      ;; analyze all DITs
      for d=0,ndits-1 do begin
         cube = readfits(dits[d],hdr)

         ;; size and wavelength
         taille  = (size(cube,/dim))[0]
         nlambda = (size(cube,/dim))[2]
         
         print,' * '+file_basename(dits[d])

         ;; loop over all wavelengths
         ncube = cube
         for l=0,nlambda-1 do begin
            frame = cube[*,*,l]

            if (nsigma gt 0) then begin
               ;; use sigma-clipping to identify bad pixels
               frame_out = sigma_filter(frame,7,n_sigma=nsigma,/iter)
               frame_out[where(frame eq 0)] = 0
               
               ;; correct them using maskinterp, which usually provides
               ;; better correction 
               bpm = frame eq frame_out
               frame_out = maskinterp(frame,bpm,1,5,'csplinterp')
            endif else frame_out = frame
            
            ncube[*,*,l] = frame_out
         endfor         

         ;; save cleaned data
         ditnum = string(d,format='(I05)')
         writefits,root+'/products/'+data_file[f]+'_'+ditnum+'_clean.fits',ncube,hdr
         mwrfits,reform(info[d]),root+'/products/'+data_file[f]+'_'+ditnum+'_params.fits',/create
      endfor
   endfor
endif

;;
;; frame centers determination and wavelength recalibration
;;
if keyword_set(do_center) then begin
   ;;
   ;; frame centers determination
   ;;
   print,'Determination of the star center and wavelength scaling'
   cube = readfits(root+'/products/'+cent_file+'_00000_clean.fits',hdr)
   
   ;; wavelength best guess
   wave_min = sxpar_eso(hdr,'HIERARCH ESO DRS IFS MIN LAMBDA')
   wave_max = sxpar_eso(hdr,'HIERARCH ESO DRS IFS MAX LAMBDA')
   wave_step = (wave_max-wave_min) / (nlambda-1)
   lambda = wave_min + wave_step * findgen(nlambda)   
   idx    = indgen(39)
   
   ;; useful parameters
   taille = (size(cube,/dim))[0]
   lsurD  = lambda*1d-6/8D*180/!dpi*3600*1000/pixel

   ;; precise center
   spot_centers = fltarr(nlambda,4,2)
   centers  = fltarr(nlambda,2)
   distance = fltarr(nlambda,6)
   window,0,xs=taille,ys=taille
   for l=0,nlambda-1 do begin
      frame = cube[*,*,l]

      ;; orientation and separation
      waffle_orientation = strtrim(sxpar_eso(hdr,'HIERARCH ESO OCS WAFFLE ORIENT'),2)     
      if (waffle_orientation eq '+') then offset = !pi/4D else offset = 0D
      orient = 57 * !pi / 180 + offset
      freq   = 10 * sqrt(2) * 0.97
      ext    = 8

      ;; center
      if (l eq 0) then begin
         ;; initial guess at the frame center
         ;;   ==> can be modified manually if necessary
         cx_int = fix(taille / 2.)
         cy_int = fix(taille / 2.)
      endif else begin
         ;; subsequent guess is center of previous channel
         cx_int = fix(centers[l-1,0])
         cy_int = fix(centers[l-1,1])
      endelse
      
      ;;
      ;; sattelite spots fitting
      ;;
      
      ;; spot 0
      cx0_int = fix(cx_int+freq*lsurD[l]*cos(orient))
      cy0_int = fix(cy_int+freq*lsurD[l]*sin(orient))
      
      sub0 = frame[cx0_int-ext:cx0_int+ext,cy0_int-ext:cy0_int+ext]
      fit0 = mpfit2dpeak(sub0,a)
      
      cx0 = a[4]-ext+cx0_int
      cy0 = a[5]-ext+cy0_int

      ;; spot 1
      cx1_int = fix(cx_int+freq*lsurD[l]*cos(orient+!pi/2))
      cy1_int = fix(cy_int+freq*lsurD[l]*sin(orient+!pi/2))
      
      sub1 = frame[cx1_int-ext:cx1_int+ext,cy1_int-ext:cy1_int+ext]
      fit1 = mpfit2dpeak(sub1,a)
      
      cx1 = a[4]-ext+cx1_int
      cy1 = a[5]-ext+cy1_int
      
      ;; spot 2
      cx2_int = fix(cx_int+freq*lsurD[l]*cos(orient+!pi))
      cy2_int = fix(cy_int+freq*lsurD[l]*sin(orient+!pi))
         
      sub2 = frame[cx2_int-ext:cx2_int+ext,cy2_int-ext:cy2_int+ext]
      fit2 = mpfit2dpeak(sub2,a)
      
      cx2 = a[4]-ext+cx2_int
      cy2 = a[5]-ext+cy2_int
      
      ;; spot 3
      cx3_int = fix(cx_int+freq*lsurD[l]*cos(orient+3*!pi/2))
      cy3_int = fix(cy_int+freq*lsurD[l]*sin(orient+3*!pi/2))
      
      sub3 = frame[cx3_int-ext:cx3_int+ext,cy3_int-ext:cy3_int+ext]
      fit3 = mpfit2dpeak(sub3,a)
      
      cx3 = a[4]-ext+cx3_int
      cy3 = a[5]-ext+cy3_int

      spot_centers[l,0,*] = [cx0,cy0]
      spot_centers[l,1,*] = [cx1,cy1]
      spot_centers[l,2,*] = [cx2,cy2]
      spot_centers[l,3,*] = [cx3,cy3]
      
      ;;
      ;; final center
      ;;
      lint,[cx0,cy0],[cx2,cy2],[cx1,cy1],[cx3,cy3],center
      centers[l,*] = center

      ;;
      ;; final scaling
      ;;
      distance[l,0] = sqrt((cx0-cx2)^2+(cy0-cy2)^2)
      distance[l,1] = sqrt((cx1-cx3)^2+(cy1-cy3)^2)
      
      distance[l,2] = sqrt((cx0-cx1)^2+(cy0-cy1)^2)         
      distance[l,3] = sqrt((cx0-cx3)^2+(cy0-cy3)^2)         

      distance[l,4] = sqrt((cx1-cx2)^2+(cy1-cy2)^2)         

      distance[l,5] = sqrt((cx2-cx3)^2+(cy2-cy3)^2)         

      ;;
      ;; display
      ;;
      loadct,0,/silent
      tvscl,sqrt(frame)
      loadct,39,/silent
      plots,[cx0,cx2],[cy0,cy2],color=250,/device
      plots,[cx1,cx3],[cy1,cy3],color=250,/device
      ;;wait,0.1
   endfor
   print

   ;; average wavelength-scaling factor (used for spectral
   ;; differential imaging)
   scale = distance / (distance[0,*] ## replicate(1,nlambda))
   scale_average = total(scale,2)/6

   ;;
   ;; wavelength recalibration
   ;;
   print,'Wavelength recalibration'
   wave = readfits(root+'/products/'+wave_file+'_preproc_col_mean_bp_ct_00000.fits',hdr)
   
   ;; cube size
   taille = (size(wave,/dim))[0]

   ;; measure flux in all channels
   aper = fltarr(taille,taille)
   aper[(taille-150)/2:(taille+150)/2-1,(taille-150)/2:(taille+150)/2-1] = 1
   nidx  = where(aper ne 0)
   flux = fltarr(nlambda)
   for l=0,nlambda-1 do begin
      frame = wave[*,*,l]      
      flux[l] = median(frame[nidx])
   endfor

   ;; fitting procedure to improve wavelength calibration
   ;; see Vigan et al. (submitted) for details
   mode = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS2 COMB IFS'),2)
   lambda_drh = lambda
   if (mode eq 'OBS_YJ') then begin
      ;; peak 1
      sub_idx  = idx[0:10]
      sub_flux = flux[0:10]
      fit = mpfitpeak(sub_idx,sub_flux,par1,nterms=5)
      
      ;; peak 2
      sub_idx  = idx[10:26]
      sub_flux = flux[10:26]
      fit = mpfitpeak(sub_idx,sub_flux,par2,nterms=5)

      ;; peak 3
      sub_idx  = idx[26:*]
      sub_flux = flux[26:*]
      fit = mpfitpeak(sub_idx,sub_flux,par3,nterms=5)

      ;; lasers
      lam_laser = [0.9877,1.1237,1.3094]
      pos_laser = [par1[1],par2[1],par3[1]]
      
      ;; fit parameter
      ref_idx = 2
      ref_lam = 0.9

      parinfo = replicate({value:0.D,fixed:0},1)            
      parinfo[0].value = ref_lam

      cscale = scale_average
      cscale = min(scale,dim=2)
      fargs = {lam_laser:lam_laser,pos_laser:pos_laser, $
               ref_idx:ref_idx,nlambda:nlambda, $
               scale_average:cscale,flux:flux,disp:1}
      param = mpfit('sph_ifs_wavelength_optimisation',parinfo.value,functargs=fargs, $
                    parinfo=parinfo,status=status,quiet=1)
      ref_lam = param[0]

      ;; final wavalength
      lambda = replicate(ref_lam,nlambda) * scale_average/scale_average[ref_idx]
   endif else begin
      ;; peak 1
      sub_idx  = idx[0:7]
      sub_flux = flux[0:7]
      fit = mpfitpeak(sub_idx,sub_flux,par1,nterms=5)
      
      ;; peak 2
      sub_idx  = idx[5:16]
      sub_flux = flux[5:16]
      fit = mpfitpeak(sub_idx,sub_flux,par2,nterms=5)

      ;; peak 3
      sub_idx  = idx[14:25]
      sub_flux = flux[14:25]
      fit = mpfitpeak(sub_idx,sub_flux,par3,nterms=5)

      ;; peak 4
      sub_idx  = idx[25:*]
      sub_flux = flux[25:*]
      fit = mpfitpeak(sub_idx,sub_flux,par4,nterms=5)

      ;; lasers
      lam_laser = [0.9877,1.1237,1.3094,1.5451]
      pos_laser = [par1[1],par2[1],par3[1],par4[1]]
      
      ;; fit parameter
      ref_idx = 2
      ref_lam = 0.9

      parinfo = replicate({value:0.D,fixed:0},1)            
      parinfo[0].value = ref_lam

      cscale = scale_average
      fargs = {lam_laser:lam_laser,pos_laser:pos_laser, $
               ref_idx:ref_idx,nlambda:nlambda, $
               scale_average:cscale,flux:flux,disp:1}
      param = mpfit('sph_ifs_wavelength_optimisation',parinfo.value,functargs=fargs, $
                    parinfo=parinfo,status=status,quiet=1)
      ref_lam = param[0]

      ;; final wavalength
      lambda = replicate(ref_lam,nlambda) * scale_average/scale_average[ref_idx]
   endelse
      
   ;; saving
   writefits,root+'/products/'+'data_centers.fits',centers
   writefits,root+'/products/'+'data_wavelength.fits',lambda
   writefits,root+'/products/'+'data_scaling.fits',scale_average
   
   ;;
   ;; summary plot
   ;;
   window,1,xs=1200,ys=400
   loadct,39,/silent
   !p.multi = [0,3,1]
   
   plot,centers[*,0],centers[*,1],psym=1,/yno,/iso, $
        xs=1,xr=cx_int+[-5,5],ys=1,yr=cy_int+[-5,5], $
        xtitle='x center [pix]',ytitle='y center [pix]', $
        title=fname,/nodata,charsize=2
   for l=0,nlambda-1 do plots,centers[l,0],centers[l,1],psym=1,color=float(l)/nlambda*250
   
   plot,idx,scale[*,0],/nodata,/yno,title=fname, $
        xtitle='Spectral channel index',ytitle='Scaling factor', $
        ys=1,yr=[0.95,1.05*max(scale_average)],xs=1,xr=[-1,nlambda+1], $
        charsize=2
   for i=0,5 do oplot,idx,scale[*,i],psym=1,color=50+40*i
   oplot,idx,scale_average,linestyle=0

   plot,lambda,flux,psym=-4,xtitle='Wavelength [micron]',charsize=2
   oplot,lambda_drh,flux,linestyle=1
   plots,0.9877,!y.crange,color=250,linestyle=2
   plots,1.1237,!y.crange,color=250,linestyle=2
   plots,1.3094,!y.crange,color=250,linestyle=2
   plots,1.5451,!y.crange,color=250,linestyle=2
   
   !p.multi = 0
   loadct,0,/silent
endif

if keyword_set(do_combine) then begin
   ;; template for frame information
   template_frame = {template_frame_ifs,     $
                     file:'',                $
                     img:-1L,                $
                     dit:-1D,                $
                     time:-1D,               $
                     time_start:-1D,         $
                     time_end:-1D,           $
                     ha:!values.d_nan,       $
                     ha_start:!values.d_nan, $
                     ha_end:!values.d_nan,   $
                     pa:!values.d_nan,       $
                     pa_start:!values.d_nan, $
                     pa_end:!values.d_nan,   $
                     seeing:!values.d_nan,   $
                     pupoff:!values.d_nan,   $
                     filtname:'',            $
                     filtswap:0,             $
                     nwave:39,               $
                     nd:!values.f_nan}
   
   ;; centers, wavelength and scaling
   centers = mrdfits(root+'/products/'+'data_centers.fits',/silent)
   lambda  = mrdfits(root+'/products/'+'data_wavelength.fits',/silent)
   nlambda = n_elements(lambda)
   lsurD   = lambda*1d-6/8D*180/!dpi*3600*1000/pixel
      
   ;;
   ;; read and format psf data (FLUX)
   ;;
   print,'Off-axis PSF centering'
   nflux_file = n_elements(flux_file)
   cube_psf = fltarr(ntaille,ntaille,nlambda,nflux_file)
   for f=0,nflux_file-1 do begin
      print,' * '+flux_file[f]+'_00000_clean.fits'
      
      ;; read data
      cube = readfits(root+'/products/'+flux_file[f]+'_00000_clean.fits',hdr)
      
      ;; DIT value
      DIT_psf = sxpar_eso(hdr,'HIERARCH ESO DET SEQ1 DIT')
      
      ;; neutral density
      nd_psf = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS4 FILT2 NAME'),2)
      ampli_psf = sph_cpi_nd_transmission(nd_psf,lambda*1000.)
      
      ;; size and center
      taille = (size(cube,/dim))[0]
      
      ;; work with even-sized cubes
      taille = (taille mod 2) ? (taille-1) : taille
      cube = cube[0:taille-1,0:taille-1,*]
      
      ;;
      ;; recenter and PSF resize data
      ;;
      for l=0,nlambda-1 do begin
         frame = cube[*,*,l]
         
         ;; normalization
         frame = frame / DIT_psf * 10^ampli_psf[l]

         ;;
         ;; recenter and resize data
         ;;
         c = peak_center(frame,hide=50,ext=8)

         ttaille = ntaille+4
         cc = (ttaille-1)/2.
         
         shiftx = cc-c[0]
         shifty = cc-c[1]
         
         shiftx_int = fix(shiftx)
         shifty_int = fix(shifty)
         
         shiftx_rem = shiftx - shiftx_int
         shifty_rem = shifty - shifty_int
         
         nframe = frame
         nframe = shift(nframe,shiftx_int,shifty_int)
         nframe = nframe[0:ttaille-1,0:ttaille-1]

         if keyword_set(fft_shift) then begin
            ;; FFT-based shift
            ;;nframe = subpixel_shift(nframe,xs=shiftx_rem,ys=shifty_rem)
            message,'No FFT-based shifting routine is provided. Please use your own.'
         endif else begin
            ;; interpolation-based
            nframe = translate(nframe,shiftx_rem,shifty_rem)
         endelse
         
         nframe = nframe[2:ttaille-3,2:ttaille-3]

         cube_psf[*,*,l,f] = nframe
      endfor
   endfor
   print
      
   ;; save data
   writefits,root+'/products/'+'data_cube_psf.fits',cube_psf

   ;;
   ;; read and format corono data
   ;;
   print,'Science data centering'

   ;; first count the number of cubes
   nfiles = n_elements(coro_file)
   ncubes = 0
   for c=0,nfiles-1 do begin
      dits    = file_search(root+'/products/'+coro_file[c]+'*[0-9]_clean.fits',count=ndits)
      ncubes += ndits
   endfor

   ;; loop on all cubes
   frames = replicate(create_struct(name='template_frame_ifs'),ncubes)
   cube_noscl = fltarr(ntaille,ntaille,nlambda,ncubes)
   cubeidx    = 0
   for f=0,nfiles-1 do begin
      print,' * '+coro_file[f]

      ;; loop on all DITs in the cubes
      dits = file_search(root+'/products/'+coro_file[f]+'*[0-9]_clean.fits',count=ndits)
      for d=0,ndits-1 do begin
         print,'   ==> '+file_basename(dits[d])
         
         ;; read data
         cube = readfits(dits[d],hdr)
         info = mrdfits(strmid(dits[d],0,strlen(dits[d])-10)+'params.fits',1,/silent)

         ;; DIT value
         DIT_coro = sxpar_eso(hdr,'HIERARCH ESO DET SEQ1 DIT')

         ;; neutral density
         nd_coro = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS4 FILT2 NAME'),2)
         ampli_coro = sph_cpi_nd_transmission(nd_coro,lambda*1000.)
         
         ;; size and center
         taille = (size(cube,/dim))[0]
         
         ;; work with even-sized cubes
         taille = (taille mod 2) ? (taille-1) : taille
         cube = cube[0:taille-1,0:taille-1,*]
         
         ;; frame info
         frames[cubeidx] = info
         
         ;; spatial rescaling
         for l=0,nlambda-1 do begin
            frame = cube[*,*,l]
            
            ;; normalization
            frame = frame / DIT_coro * 10^ampli_coro[l]
            
            ;;
            ;; recenter data
            ;;
            ttaille = ntaille+4
            cc = (ttaille-1)/2.
            
            shiftx = cc-centers[l,0]
            shifty = cc-centers[l,1]
         
            shiftx_int = fix(shiftx)
            shifty_int = fix(shifty)
            
            shiftx_rem = shiftx - shiftx_int
            shifty_rem = shifty - shifty_int

            nframe = frame
            nframe = shift(nframe,shiftx_int,shifty_int)
            nframe = nframe[0:ttaille-1,0:ttaille-1]

            if keyword_set(fft_shift) then begin
               ;; FFT-based shift
               ;; nframe = subpixel_shift(nframe,xs=shiftx_rem,ys=shifty_rem)
               message,'No FFT-based shifting routine is provided. Please use your own.'
            endif else begin
               ;; interpolation-based
               nframe = translate(nframe,shiftx_rem,shifty_rem)
            endelse            
            
            nframe = nframe[2:ttaille-3,2:ttaille-3]
            
            ;; data without scaling
            cube_noscl[*,*,l,cubeidx] = nframe
         endfor     ;; lambda
         cubeidx++
      endfor        ;; DITs
   endfor           ;; file
   print
   
   ;; parallactic angle correction for VLT: correction is needed
   ;; because when some targets cross meridian, there is a
   ;; discontinuity (going from -180 to +180). However, even if the
   ;; discontinuity is not corrected, the derotation will be right
   deriv = (frames.pa-shift(frames.pa,1))[1:*]
   if (max(abs(deriv)) gt 180.) then frames[where(frames.pa lt 0)].pa += 360
   
   deriv = (frames.pa_start-shift(frames.pa_start,1))[1:*]
   if (max(abs(deriv)) gt 180.) then frames[where(frames.pa_start lt 0)].pa_start += 360
   
   deriv = (frames.pa_end-shift(frames.pa_end,1))[1:*]
   if (max(abs(deriv)) gt 180.) then frames[where(frames.pa_end lt 0)].pa_end += 360

   ;; save data
   writefits,root+'/products/'+'data_cube_coro.fits',cube_noscl
   mwrfits,frames,root+'/products/'+'data_info.fits',/create   
endif

if keyword_set(do_erase) then begin
   ;; ask whether the user is sure he wants to erase the data
   ok = 0
   er = 0
   while (ok eq 0) do begin
      print,'Are you sure you want to erase all the intermediary files? (*_clean.fits, *_parames.fits)'
      a = ''
      read,a,prompt='[y/n] > '
      case a of
         'Y': begin
            ok = 1
            er = 1
         end
         'y': begin
            ok = 1
            er = 1
         end
         'N': begin
            ok = 1
            er = 0
         end
         'n': begin
            ok = 1
            er = 0
         end
         else: begin
            print,'Please answer y or n'
            print
         end
      endcase
   endwhile

   ;; if yes, erase the data
   if (er ne 0) then begin
      print
      print,'Erasing files...'

      files = file_search(root+'/products/'+'*[0-9]_clean.fits',count=nfiles)
      if (nfiles ne 0) then file_delete,files

      files = file_search(root+'/products/'+'*[0-9]_params.fits',count=nfiles)
      if (nfiles ne 0) then file_delete,files
   endif
endif

fin:

end
