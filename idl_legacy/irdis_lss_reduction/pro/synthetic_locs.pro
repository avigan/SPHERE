;+
; NAME:
;
;  SYNTHETIC_LOCS
;
; AUTHOR:
;
;   Arthur Vigan
;   Laboratoire d'Astrophysique de Marseille
;   38, rue Frederic Joliot-Curie, 13388 Marseille Cedex 13, France
;   arthur.vigan@lam.fr
;
; PURPOSE:
;
;  Locally Optimized Combination of Spectra
;  Performs speckle estimation to create a "synthetic" spectrum using
;  an application of the LOCI for spectral deconvolution of IRDIS
;  spectra
;
; CALLING SEQUENCE:
;
;   synth = SYNTHETIC_LOCS(sig_over,lambda,use_mask,pla_mask, $
;                          PIXEL=pixel,LAMBDA_REF=lambda_ref, $
;                          LOCS_CRIT=locs_crit,SILENT=silent)
;
; DESCRIPTION:
;
;  This function is designed to estimate the speckles in the rescaled
;  spectrum. It will create a synthetic spectrum that represent the
;  speckles and stellar halo. It uses a method based on the Locally
;  Optimized Combination of Images from Lafreniere et al. (2007),
;  adapted to work for the spectral deconvolution of IRDIS spectra
;
; INPUTS:
;
;  SIG_OVER - 2D rescaled (and sometimes oversampled) spectrum
;
;  LAMBDA - wavelength vector; in nanometers
;
;  PLANET_MASK - binary mask for the planet, as created by the
;                SPH_IRD_LSS_RESCALE function
;
;  USE_MASK - binary mask of the "usable" area, as created by the
;             SPH_IRD_LSS_RESCALE function
;
;  PIXEL - pixel scale; in arcsecond/pixel
;
;  LAMBDA_REF - reference wavelength used for the rescaling, that is
;               here used to calculate the size of the exclusion
;               window arround the planet
;
;  SPECTRAL - performs a spectral fit for the speckles estimation
;
;  SPATIAL - performs a spatial fit for the speckles estimation
;
;  LOCS_CRIT - size of the exclusion window arround the planet
;
;  SILENT - suppress verbose outputs
;
; OUTPUTS:
;
;  A synthetic spectrum representing the speckles and stellar halo in
;  the input SIG_OVER spectrum.
;
; MODIFICATION HISTORY:
;
;  28/01/2016 - arthur.vigan - implemented a spatial version of the fit
;
;  22/01/2013 - arthur.vigan@lam.fr
;  imp: added descriptive header with proper documentation
;
;  08/11/2012 - arthur.vigan@lam.fr
;  add: SILENT keyword
;
;  24/07/2012 - arthur.vigan@lam.fr
;  Created
;
;-

function synthetic_locs_spectral,sig_over,lambda,use_mask,pla_mask, $
                                 PIXEL=pixel,LAMBDA_REF=lambda_ref, $
                                 PLANET_SEP=pla_sep,PLANET_SIZE=pla_size, $
                                 LOCS_CRIT=locs_crit,SILENT=silent, $
                                 OVERSAMPLE=scale,NORMALIZE=normalize
  if ~keyword_set(lambda_ref) then lambda_ref = lambda[0]  
  if ~keyword_set(pixel)      then pixel = 0.01225
  if ~keyword_set(locs_crit)  then locs_crit = 1.0

  nlambda = n_elements(lambda)
  w = (size(sig_over,/dim))[0]
  h = (size(sig_over,/dim))[1]
  bigh = (size(sig_over,/dim))[1]
  sep = indgen(bigh) - (bigh-1)/2.

  all_mask = use_mask*pla_mask

  signal = sig_over
  
  ;; normalization
  if keyword_set(normalize) then begin
     norm = synthetic_vigan2008(sig_over,lambda,use_mask,pla_mask,/SILENT)
     signal = signal / norm
  endif

  ;; remove speckles
  if ~keyword_set(silent) then print,'Speckles attenuation (LOCS)...'
  synthetic_over = fltarr(w,bigh)
  for l=0,bigh-1 do begin
;  for l=440,440 do begin
     if ~keyword_set(silent) then $
        print,cr()+' * separation '+numformat(l+1)+' / '+ $
              numformat(bigh),format='($,a)'
     spectrum = reform(signal[*,l],nlambda)
     mask     = reform(all_mask[*,l],nlambda)

     ;; ;; planet aperture
     ;; psep = pla_sep/pixel*lambda_ref/lambda*scale+(h-1)/2.
     ;; pla_r_lin = pla_size*lambda_ref*1d-9/8D*180/!pi*3600/pixel*scale
     ;; pla_r_min = psep-pla_r_lin/2.0
     ;; pla_r_max = psep+pla_r_lin/2.0

     ;; pmin = pla_r_min[l]
     ;; pmax = pla_r_max[l]

     ;; ;; usable areas
     ;; lmin = max(where(pla_r_max lt pmin))
     ;; lmax = min(where(pla_r_min gt pmax))

     ;; loadct,39,/silent
     ;; plot,lambda,psep,/yno,xs=1
     ;; oplot,lambda,pla_r_min,color=250
     ;; oplot,lambda,pla_r_max,color=50
     ;; plots,!x.crange,l,linestyle=2
     ;; plots,lambda[lmin],!y.crange,linestyle=2
     ;; plots,lambda[lmax],!y.crange,linestyle=2
     
     ;; stop
     ;; hak
     ;; continue
     
     ;; reference spectra --> at xx lambda/D away from the current wavelength
     lsd = ceil(lambda_ref*1d-9/8.0*180/!pi*3600/pixel/2*locs_crit)
     idx = indgen(bigh)
     tmp = total(all_mask,1)
     iref = where(((idx lt (l-lsd)) or ((l+lsd) lt idx)) and (tmp ge max(tmp)),nref)     
     ref  = signal[*,iref]

     good = where(mask eq 1,ngood,comp=bad,ncomp=nbad)
     if (ngood eq 0) then continue
     
     ref_pix = ref[good,*]
     obj_pix = spectrum[good]
     
     AA = ref_pix ## transpose(ref_pix)
     BB = reform(ref_pix ## transpose(obj_pix),nref)

     ;; solve linear system with positive coefficients only
     bnd = fltarr(2,nref)
     bnd[1,*] = 1e20
     bvls,AA,BB,bnd,coeff
     ;; coeff = la_linear_equation(AA,BB,status=status)
     
     ii = where(coeff ne 0,nii)
     print,nii,n_elements(coeff)
     
     ;; synthetic spectrum
     synthetic_over[*,l] = ref # coeff
  endfor
  if ~keyword_set(silent) then begin
     print,cr()+'  --> done!',format='($,a)'
     print
  endif

  ;; inverse normalization
  if keyword_set(normalize) then begin
     synthetic_over *= norm
  endif

  return,synthetic_over
end

function synthetic_locs_spatial,sig_over,lambda,use_mask,pla_mask, $
                                PIXEL=pixel,LAMBDA_REF=lambda_ref, $
                                PLANET_SEP=pla_sep,PLANET_SIZE=pla_size, $
                                LOCS_CRIT=locs_crit,SILENT=silent, $
                                OVERSAMPLE=scale
  if ~keyword_set(lambda_ref) then lambda_ref = lambda[0]  
  if ~keyword_set(pixel)      then pixel = 0.01225
  if ~keyword_set(locs_crit)  then locs_crit = 1.0

  nlambda = n_elements(lambda)
  w = (size(sig_over,/dim))[0]
  h = (size(sig_over,/dim))[1]
  bigh = (size(sig_over,/dim))[1]
  sep = indgen(bigh) - (bigh-1)/2.

  all_mask = use_mask*pla_mask

  signal = sig_over
  
  ;; normalization
  if keyword_set(normalize) then begin
     norm = synthetic_vigan2008(sig_over,lambda,use_mask,pla_mask,/SILENT)
     signal = signal / norm
  endif

  pix = findgen(bigh)
  
  ;; remove speckles
  if ~keyword_set(silent) then print,'Speckles attenuation (LOCS)...'
  synthetic_over = fltarr(w,bigh)
  for l=0,nlambda-1 do begin
     if ~keyword_set(silent) then $
        print,cr()+' * wavelength '+numformat(l+1)+' / '+ $
              numformat(nlambda),format='($,a)'
     spectrum = reform(signal[l,*],bigh)
     mask     = reform(all_mask[l,*],bigh)
     mask[bigh/2:*] = 0
     
     ;; planet aperture
     psep = pla_sep/pixel*lambda_ref/lambda*scale+(h-1)/2.
     pla_r_lin = pla_size*lambda_ref*1d-9/8D*180/!pi*3600/pixel*scale
     pla_r_min = psep-pla_r_lin/2.0
     pla_r_max = psep+pla_r_lin/2.0

     pmin = pla_r_min[l]
     pmax = pla_r_max[l]

     ;; usable areas
     lmin = max(where(pla_r_max lt pmin))
     lmax = min(where(pla_r_min gt pmax))

     ldown = lmin
     lup   = nlambda-lmax
     
     ;; plot
     ;; loadct,39,/silent
     ;; plot,lambda,psep,/yno,xs=1
     ;; oplot,lambda,pla_r_min,color=250
     ;; oplot,lambda,pla_r_max,color=50
     ;; plots,lambda[l],psep[l],linestyle=2,psym=4
     ;; if (lmin ne -1) then plots,lambda[lmin],!y.crange,linestyle=2
     ;; if (lmax ne -1) then plots,lambda[lmax],!y.crange,linestyle=2
     ;; hak
     
     if (ldown ge lup) or (lmax eq -1) then begin
        idx = indgen(lmin)
        ref = signal[idx,*]
        nref = n_elements(idx)
        good = pix[where((pix ge pmin) and (mask ne 0),ngood)]

        ref_pix = ref[*,good]
        obj_pix = spectrum[good]
        
        AA = transpose(ref_pix) ## ref_pix
        BB = reform(transpose(ref_pix) ## obj_pix,nref)

        ;; solve linear system with positive coefficients only
        bnd = fltarr(2,nref)
        bnd[1,*] = 1e20
        bvls,AA,BB,bnd,coeff,NSETP=NSETP
        
        plot,pix,spectrum,title='Down'
        plots,psep[l],!y.crange,linestyle=1
        oplot,pix[good],spectrum[good],color=250
        oplot,pix,coeff # ref,color=200

        print,'down',nsetp
        
        ;; synthetic spectrum
        synthetic_over[l,*] = coeff # ref
     endif else if (ldown lt lup) then begin
        idx = lmax+indgen(lup)
        ref = signal[idx,*]
        nref = n_elements(idx)
        good = pix[where((pix le pmax) and (mask ne 0),ngood)]

        ref_pix = ref[*,good]
        obj_pix = spectrum[good]
        
        AA = transpose(ref_pix) ## ref_pix
        BB = reform(transpose(ref_pix) ## obj_pix,nref)

        ;; solve linear system with positive coefficients only
        bnd = fltarr(2,nref)
        bnd[1,*] = 1e20
        bvls,AA,BB,bnd,coeff,NSETP=NSETP
        
        plot,pix,spectrum,title='Up'
        plots,psep[l],!y.crange,linestyle=1
        oplot,pix[good],spectrum[good],color=250
        oplot,pix,coeff # ref,color=200

        print,'up',nsetp
        
        ;; synthetic spectrum
        synthetic_over[l,*] = coeff # ref
     endif
     
     ;; stop
     ;; hak & continue     
  endfor
  if ~keyword_set(silent) then begin
     print,cr()+'  --> done!',format='($,a)'
     print
  endif

  ;; inverse normalization
  if keyword_set(normalize) then begin
     synthetic_over *= norm
  endif

  return,synthetic_over
end

function synthetic_locs,sig_over,lambda,use_mask,pla_mask, $
                        PIXEL=pixel,LAMBDA_REF=lambda_ref, $
                        PLANET_SEP=pla_sep,PLANET_SIZE=pla_size, $
                        LOCS_CRIT=locs_crit,SILENT=silent, $
                        OVERSAMPLE=scale,NORMALIZE=normalize, $
                        SPATIAL=spatial,SPECTRAL=spectral
  synthetic = 0
  
  ;; default is SPATIAL
  if ~keyword_set(spatial) and ~keyword_set(spectral) then spectral = 1

  if keyword_set(spatial) and keyword_set(spectral) then $
     message,'You cannot set both SPATIAL and SPECTRAL. You have to choose on (default=SPECTRAL).'

  if keyword_set(spatial) then begin
     synthetic = synthetic_locs_spatial(sig_over,lambda,use_mask,pla_mask, $
                                        PIXEL=pixel,LAMBDA_REF=lambda_ref, $
                                        PLANET_SEP=pla_sep,PLANET_SIZE=pla_size, $
                                        LOCS_CRIT=locs_crit,SILENT=silent, $
                                        OVERSAMPLE=scale)
  endif

  if keyword_set(spectral) then begin
     synthetic = synthetic_locs_spectral(sig_over,lambda,use_mask,pla_mask, $
                                         PIXEL=pixel,LAMBDA_REF=lambda_ref, $
                                         PLANET_SEP=pla_sep,PLANET_SIZE=pla_size, $
                                         LOCS_CRIT=locs_crit,SILENT=silent, $
                                         OVERSAMPLE=scale,NORMALIZE=normalize)
  endif

  return,synthetic
end
