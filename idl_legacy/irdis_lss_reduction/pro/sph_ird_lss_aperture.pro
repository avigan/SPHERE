;+
; NAME:
;
;  SPH_IRD_LSS_APERTURE
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
;  Creates an aperture to perform spectral extraction
;
; CALLING SEQUENCE:
;
;  ap = SPH_IRD_LSS_APERTURE(sig,pla_sep,lambda,APER_DIAM=aper,PIXEL=pixel)
;
; DESCRIPTION:
;
;  This function creates an aperture of fixed diameter in lambda/D. It
;  can then be used for the spectral extraction of the reference star,
;  the atmospheric calibrator or the planet.
;
;  Two apertures are created: one in the normal data space, and one in
;  the spatialy rescaled space. This second aperture can be used on
;  speckle-subtracted data to avoid having to perform another
;  transformation step to the data.
;
; INPUTS:
;
;  SIG - 2D spectrum used to determine the size of final array created
;        by the function
;
;  LAMBDA - wavelength vector; in nanometers
;
;  PLANET_SEP - separation of the planet; in arcseconds
;  
; KEYWORD PARAMETERS:
;
;  APER_DIAM - diameter of the aperture at each wavelength; in lambda/D
;
;  PIXEL - pixel scale; in arcsecond/pixel
;
;  LAMBDA_REF - reference wavelength for the spatial rescaling; in
;               nanometers
;               
;  APER_RESCALED - photometric aperture in the spatially rescaled
;                  space
;
;  UPSCALED - if set, the APER_RESCALED aperture is actually scaled
;             according to the original space
;
; OUTPUTS:
;
;  Returns an array of the same size as the input spectrum, containing
;  the aperture at the proper separation.
;
; MODIFICATION HISTORY:
;
;  16/10/2014 - arthur.vigan@lam.fr
;  add: UPSCALE keyword and corresponding calculations
;
;  23/01/2013 - arthur.vigan@lam.fr
;  add: added calculation of the photometric aperture in the spatially
;       rescaled space
;
;  22/01/2013 - arthur.vigan@lam.fr
;  imp: changed some keywords to improve the name
;  imp: added descriptive header with proper documentation
;
;  30/07/2012 - arthur.vigan@lam.fr
;  Created
;
;-

function sph_ird_lss_aperture,sig,lambda,planet_sep,APER_DIAM=aper,PIXEL=pixel, $
                              LAMBDA_REF=lambda_ref,APER_RESCALED=aper_rescaled, $
                              UPSCALE=upscale
  if ~keyword_set(aper)       then aper  = 1.0
  if ~keyword_set(pixel)      then pixel = 0.01225
  if ~keyword_set(lambda_ref) then lambda_ref = lambda[0]
  
  nlambda = n_elements(lambda)
  h = (size(sig,/dim))[1]

  ;; standard aperture
  sep = planet_sep/pixel+(h-1)/2.
  aper_lin = aper*lambda*1d-9/8D*180/!pi*3600/pixel
  aper_min = sep-aper_lin/2.0
  aper_min_int = ceil(aper_min)
  aper_min_rem = aper_min_int-aper_min
  aper_max = sep+aper_lin/2.0
  aper_max_int = floor(aper_max)
  aper_max_rem = aper_max-aper_max_int
  
  ;; rescaled aperture
  if ~keyword_set(upscale) then begin
     sep = planet_sep/pixel*lambda_ref/lambda+(h-1)/2.
     aper_r_lin = aper*lambda_ref*1d-9/8D*180/!pi*3600/pixel
  endif else begin
     sep = planet_sep/pixel*lambda/lambda_ref+(h-1)/2.
     aper_r_lin = aper*lambda*1d-9/8D*180/!pi*3600/pixel
  endelse
  aper_r_min = sep-aper_r_lin/2.0
  aper_r_min_int = ceil(aper_r_min)
  aper_r_min_rem = aper_r_min_int-aper_r_min
  aper_r_max = sep+aper_r_lin/2.0
  aper_r_max_int = floor(aper_r_max)
  aper_r_max_rem = aper_r_max-aper_r_max_int

  ;; loadct,39,/silent
  ;; plot,lambda,[0],ys=1,yr=[min(aper_r_min)-5,max(aper_max)+5]
  ;; oplot,lambda,floor(aper_min),color=150
  ;; oplot,lambda,aper_min_int,color=50
  ;; oplot,lambda,aper_min,color=250
  ;; oplot,lambda,aper_max,color=250
  ;; oplot,lambda,aper_max_int,color=50
  ;; oplot,lambda,ceil(aper_max),color=150
  ;; oplot,lambda,floor(aper_min),color=150
  ;; oplot,lambda,aper_r_min_int,color=50
  ;; oplot,lambda,aper_r_min,color=250
  ;; oplot,lambda,aper_r_max,color=250
  ;; oplot,lambda,aper_r_max_int,color=50
  ;; oplot,lambda,ceil(aper_r_max),color=150

  dist_lin    = findgen(h)
  aper_phot   = fltarr(nlambda,h)
  aper_r_phot = fltarr(nlambda,h)
  for l=0,nlambda-1 do begin
     ;; standard aperture
     aper_phot[l,*] = (aper_min_int[l] le dist_lin) and (dist_lin le aper_max_int[l])
     
     min = aper_min_int[l]-1
     aper_phot[l,min] = aper_min_rem[l]
     
     max = aper_max_int[l]+1
     aper_phot[l,max] = aper_max_rem[l]

     ;; rescaled aperture
     aper_r_phot[l,*] = (aper_r_min_int[l] le dist_lin) and (dist_lin le aper_r_max_int[l])
     
     min = aper_r_min_int[l]-1
     aper_r_phot[l,min] = aper_r_min_rem[l]
     
     max = aper_r_max_int[l]+1
     aper_r_phot[l,max] = aper_r_max_rem[l]     
  endfor

  aper_rescaled = aper_r_phot

  return,aper_phot  
end 
