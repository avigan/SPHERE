;+
; NAME:
;
;  SYNTHETIC_VIGAN2008
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
;  Performs speckle estimation to create a "synthetic" spectrum using
;  the method described in Vigan et al. (2008).
;
; CALLING SEQUENCE:
;
;   synth = SYNTHETIC_VIGAN2008(sig_over,lambda,use_mask,pla_mask, $
;                               LIMIT=limit,TRANSMISSION=tr,TR_LIMIT=tr_lim, $
;                               SILENT=silent)
;
; DESCRIPTION:
;
;  This function is designed to estimate the speckles in the rescaled
;  spectrum. It will create a synthetic spectrum that represent the
;  speckles and stellar halo. It uses the method described in the
;  paper Vigan, Langlois, Moutou & Dohlen, 2008, A&A, 489, 1345.
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
;  LIMIT - limit the extent to which the signal is used to create the
;          "model" spectrum that will be fitted to the speckle; in
;          pixels
;
;  TRANSMISSION - estimation of the atmosphere transmission at all
;                 wavelengths
;
;  TR_LIMIT - limit transmission below which the points will not be
;             considered in the amplitude fit of the speckles
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
;  22/01/2013 - arthur.vigan@lam.fr
;  imp: added descriptive header with proper documentation
;  add: added the TR_LIMIT keyword
;  add: added the LIMIT keyword
;
;  08/11/2012 - arthur.vigan@lam.fr
;  add: SILENT keyword
;
;  24/07/2012 - arthur.vigan@lam.fr
;  Created
;
;-

function adjust_factor_fun,p,SPECTRUM=spectrum,MODEL=model,WEIGHT=weight,_EXTRA=extra
  return,weight*(spectrum-p[0]*model)^2.
end

function synthetic_vigan2008,sig_over,lambda,use_mask,pla_mask, $
                             LIMIT=limit,TRANSMISSION=tr,TR_LIMIT=tr_lim, $
                             SILENT=silent
  nlambda = n_elements(lambda)
  w = (size(sig_over,/dim))[0]
  bigh = (size(sig_over,/dim))[1]

  if ~keyword_set(limit)  then limit = bigh
  if ~keyword_set(tr)     then tr = replicate(1.0,nlambda)
  if ~keyword_set(tr_lim) then tr_lim = 0.0
  
  all_mask = use_mask*pla_mask
  
  ;; average spectrum
  linmask = uint(total(all_mask,1))
  index = findgen(bigh)-bigh/2
  good  = where((linmask eq max(linmask)) and (abs(index) le limit))
  model = total(sig_over[*,good],2) / max(linmask)

  ;; remove speckles
  if ~keyword_set(silent) then print,'Speckles attenuation (Vigan et al. 2008)...'
  synthetic_over = fltarr(w,bigh)
  for l=0,bigh-1 do begin
     if ~keyword_set(silent) then $
        print,cr()+' * separation '+numformat(l+1)+' / '+ $
              numformat(bigh),format='($,a)'
     spectrum = reform(sig_over[*,l],nlambda)
     mask     = reform(all_mask[*,l],nlambda)     
     
     ;; bad astmospheric transmission
     bad = where(tr lt tr_lim,nbad)
     if (nbad gt 0) then mask[bad] = 0

     weight  = mask
     par     = [max(spectrum)/max(model)]
     fargs   = {spectrum:spectrum,model:model,weight:weight}
     factor  = mpfit('adjust_factor_fun',par,functargs=fargs,quiet=1)
     
     synthetic_over[*,l] = model*factor[0]

     ;; loadct,39,/silent
     ;; plot,lambda,spectrum
     ;; oplot,lambda,model*factor[0],color=150
     ;; oplot,lambda,mask*max(spectrum),color=250
  endfor
  if ~keyword_set(silent) then begin
     print,cr()+'  --> done!',format='($,a)'
     print
  endif

  return,synthetic_over
end
