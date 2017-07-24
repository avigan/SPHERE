;+
; NAME:
;
;  SPH_IRD_LSS_RESCALE
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
;  Performs spatial rescaling of IRDIS/LSS spectrum at all wavelengths
;
; CALLING SEQUENCE:
;
;  nsig = sph_ird_lss_rescale(sig,lambda,REVERSE=reverse,LAMBDA_REF=lambda_ref, $
;                             PIXEL=pixel,OVERSAMPLE=scale, $
;                             CORO_RAD=coro_rad,CORO_MASK=coro_mask, $
;                             PLANET_SEP=planet_sep,PLANET_SIZE=planet_size, $
;                             PLANET_MASK=pla_mask,USE_MASK=use_mask, $
;                             SILENT=silent)
;
; DESCRIPTION:
;
;  This function rescales an inpit 2D IRDIS LSS spectrum using the
;  wavelengt vector as reference for the scale. The rescaling can be
;  done either in a direct way (initial spectrum ==> spectrum with
;  straight speckles) or in a reverse way (spectrum with straight
;  speckles ==> initial scale). It assumes that the spectrum is
;  perfectly centered in the image along the spatial dimension.
;
;  The rescaling itself is done in an optimal way using the function
;  rescale_padding, which does zero-padding in the direct and Fourrier
;  space to obtain the best possible rescaling factor. The use of this
;  function implies that the spectrum must be an array of even size
;  along the spatial dimension.
;
;  The rescaling process is the first step to perform spectral
;  deconvolution to be used to estimate and subtract the speckles. See
;  the following papers for more details:
;
;    * Sparks & Ford, 2002, ApJ, 578, 543
;    * Vigan, Langlois, Moutou & Dohlen, 2008, A&A, 489, 1345
;
; INPUTS:
;
;  SIG - 2D spectrum to be rescaled
;
;  LAMBDA - wavelength vector; in nanometers
;
; KEYWORD PARAMETERS:
;
;  REVERSE - switch to indicate wether the direction of the rescaling
;
;  LAMBDA_REF - reference wavelength for the spatial rescaling; in
;               nanometers
;
;  PIXEL - pixel scale; in arcsecond/pixel
;
;  OVERSAMPLE - oversampling scale for the rescaling
;
;  CORO_RAD - radius of the coronagraph; in arcseconds
;
;  PLANET_SEP - separation of the planet; in arcseconds
;
;  PLANET_SIZE - size of the planet; in lambda/D. It is used to create
;                the planet mask
;
;  SILENT - suppress verbose outputs
;
; OUTPUTS:
;
;  Input signaled spatially rescaled at all wavelengths, with a size
;  determined by the oversampling parameter
;
; OPTIONAL KEYWORD PARAMETERS:
;
;  CORO_MASK - binary mask for the coronagraph. It has the same size
;              as the output rescaled signal
;
;  PLANET_MASK - binary mask for the planet. It has the same size as
;                the output rescaled signal 
;
;  USE_MASK - binary mask of the "usable" area, i.e. excluding the
;             coronagraph and the planet, in the rescaled space. It
;             has the same size as the output rescaled signal
;
; MODIFICATION HISTORY:
;
;  22/01/2013 - arthur.vigan@lam.fr
;  imp: changed some keywords to improve the name
;  fix: fixed a bug that created a non symmetric mask with respect to
;       the center
;  imp: added descriptive header with proper documentation
;
;  29/08/2013 - arthur.vigan@lam.fr
;  imp: improved to be able to handle the case of a null PSF size
;
;  08/11/2012 - arthur.vigan@lam.fr
;  add: SILENT keyword
;
;  24/07/2012 - arthur.vigan@lam.fr
;  Created
;
;-

function sph_ird_lss_rescale,sig,lambda,REVERSE=reverse,LAMBDA_REF=lambda_ref, $
                             PIXEL=pixel,OVERSAMPLE=scale, $
                             CORO_RAD=coro_rad,CORO_MASK=coro_mask, $
                             PLANET_SEP=planet_sep,PLANET_SIZE=planet_size, $
                             PLANET_MASK=pla_mask,USE_MASK=use_mask, $
                             SILENT=silent
  if ~keyword_set(scale)      then scale = 1
  if ~keyword_set(lambda_ref) then lambda_ref = lambda[0]
  
  nlambda = n_elements(lambda)
  scale = float(scale)

  if ~keyword_set(reverse) then begin
     ;; necessary keywords
     if ~keyword_set(pixel)       then pixel = 0.01225
     if ~keyword_set(coro_rad)    then coro_rad = 0.0
     if ~keyword_set(planet_size) then planet_size = 0.0

     ;; dimensions
     w = (size(sig,/dim))[0]
     h = (size(sig,/dim))[1]
     bigh  = scale*h
     psf_size = planet_size*lambda_ref*1d-9/8.0*206265/pixel

     ;; coronagraphic mask
     coro_mask = replicate(1.,w,h)
     coro_mask[*,(h-1)/2.-coro_rad/pixel:(h-1)/2.+coro_rad/pixel] = 0.0     
     
     ;; rescaling and oversampling --> direct
     sig_resc = fltarr(w,bigh)
     use_mask = fltarr(w,bigh)
     pla_mask = fltarr(w,bigh)  
     if ~keyword_set(silent) then print,'Spatial rescaling...'
     for l=0,nlambda-1 do begin
        if ~keyword_set(silent) then $
           print,cr()+' * lambda = '+numformat(lambda[l],deci=1), $
                 format='($,a)'

        ;; data
        col = fltarr(bigh)
        col[(bigh-h)/2:(bigh+h)/2-1] = sig[l,*]
        fact = lambda_ref / lambda[l] * scale
        ncol = rescale_padding_1d(tab=col,zoom_io=fact,/auto)
        sig_resc[l,*] = ncol

        ;; boundaries mask
        col = replicate(1,h)
        nh  = round(lambda_ref/lambda[l]*bigh,/L64)
        ncol = congrid(col,nh,cubic=-0.5)
        tmp = fltarr(bigh)
        tmp[(bigh-nh)/2.:(bigh+nh)/2.-1] = ncol
        
        ;; usable area mask
        col = replicate(1,bigh)
        extend = floor((bigh - bigh * lambda_ref / lambda[l]) / 2)
        if (extend gt 1) then begin
           col[0:extend-1] = 0.0
           col[bigh-extend:*] = 0.0
        endif
        if (coro_rad gt 0.0) then $
           col[bigh/2-coro_rad/pixel*fact:bigh/2+coro_rad/pixel*fact] = 0.0
        use_mask[l,*] = col

        ;; planet mask
        col = replicate(1,bigh)
        if (psf_size gt 0) then begin
           nsep = lambda_ref/lambda[l]*planet_sep/pixel*scale
           nc   = (bigh-1)/2.
           col[nc+nsep-psf_size*scale/2:nc+nsep+psf_size*scale/2] = 0.0
        endif
        pla_mask[l,*] = col
     endfor
     if ~keyword_set(silent) then begin
        print,cr()+'  --> done!',format='($,a)'
        print
     endif
  endif else begin
     w = (size(sig,/dim))[0]
     bigh = (size(sig,/dim))[1]
     h = bigh/scale

     ;; rescaling and oversampling --> reverse
     sig_resc = fltarr(w,h)
     if ~keyword_set(silent) then print,'Spatial rescaling...'
     for l=0,nlambda-1 do begin
        if ~keyword_set(silent) then $
           print,cr()+' * lambda = '+numformat(lambda[l],deci=2), $
                 format='($,a)'
        
        col = reform(sig[l,*],bigh)
        fact = lambda[l] / lambda_ref / scale
        ncol = rescale_padding_1d(tab=col,zoom_io=fact,/auto)
        sig_resc[l,*] = ncol[(bigh-h)/2:(bigh+h)/2-1]        
     endfor
     if ~keyword_set(silent) then begin
        print,cr()+'  --> done!',format='($,a)'
        print
     endif
  endelse

  return,sig_resc
end
