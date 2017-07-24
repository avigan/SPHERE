;+
; NAME:
;
;  SPH_IRD_LSS_FILL
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
;  Fills the area behind the coronagraph in a 2D IRDIS spectrum with
;  continuous data
;
; CALLING SEQUENCE:
;
;   nsig = SPH_IRD_LSS_FILL(sig_over,lambda,use_mask,LAMBDA_REF=lambda_ref, $
;                           PIXEL=pixel,OVERSAMPLE=scale,CORO_RAD=coro_rad,
;                           FILLER=filler)
;
; DESCRIPTION:
;
;  This function fills the area hidden behind the coronagraph with
;  continuous data at all wavelength to avoid having discontinuites in
;  the data that will create artifacts when spatially rescaled with
;  FFT. It uses a 3rd order polynomial to fill the gap, with two data
;  points as anchors for the fit on both sides of the coronagraph.
;
; INPUTS:
;
;  SIG_OVER - 2D rescaled (and sometimes oversampled) spectrum
;
;  LAMBDA - wavelength vector; in nanometers
;
;  USE_MASK - binary mask of the "usable" area, as created by the
;             SPH_IRD_LSS_RESCALE function
;
; KEYWORD PARAMETERS:
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
;  FILLER - original data that can be used to fill the gap with
;           continuous data
;
; OUTPUTS:
;
;  Same 2D spectrum as the input, except in the areas behind the
;  coronagraph where the spurious data has been replaced with the
;  polynomial fit performed in the function.
;
; MODIFICATION HISTORY:
;
;  09/04/2015 - arthur.vigan@lam.fr
;  add: possibility to use the original data to fill the gap
;
;  22/01/2013 - arthur.vigan@lam.fr
;  imp: changed some keywords to improve the name
;  imp: added descriptive header with proper documentation
;
;  24/07/2012 - arthur.vigan@lam.fr
;  Created
;
;-

function sph_ird_lss_fill,sig_over,lambda,use_mask,LAMBDA_REF=lambda_ref, $
                          PIXEL=pixel,OVERSAMPLE=scale,CORO_RAD=coro_rad, $
                          FILLER=filler
  if ~keyword_set(pixel)      then pixel = 0.01225
  if ~keyword_set(scale)      then scale = 1
  if ~keyword_set(coro_rad)   then coro_rad = 0.2
  if ~keyword_set(lambda_ref) then lambda_ref = lambda[0]

  nlambda = n_elements(lambda)
  w = (size(sig_over,/dim))[0]
  bigh = (size(sig_over,/dim))[1]

  ;; copy input signal
  synthetic_over = sig_over
  
  ;; fill unused area with "continuous" data to avoid FFT artifacts
  if keyword_set(filler) then begin
     ;; polynomial fill
     ext0 = coro_rad/pixel * scale
     ext1 = coro_rad/pixel * lambda_ref / lambda[nlambda-1] * scale + 1
     mask = fltarr(w,bigh)
     fill_mask = use_mask
     fill_mask[*,0:bigh/2-ext0] = 1.0
     fill_mask[*,bigh/2+ext0:*] = 1.0
     for l=0,nlambda-1 do begin
        line = reform(synthetic_over[l,*]*fill_mask[l,*])
        mask = reform(fill_mask[l,*])
        
        i0 = min(where(mask eq 0))-1
        i1 = fix(bigh)/2
        i2 = max(where(mask eq 0))+1
        
        x = [i0-1,i0,i2,i2+1]
        y = line[x]
        p = poly_fit(x,y,3)
        
        nx = findgen(bigh)
        ny = p[0]+p[1]*nx+p[2]*nx^2.+p[3]*nx^3.
        
        line[i0+1:i2-1] = ny[i0+1:i2-1]   
        synthetic_over[l,*] = line
        
        ;; loadct,39,/silent
        ;; plot,nx,line,xs=1,xr=[bigh/2-50,bigh/2+50]
        ;; oplot,nx[i0:i2],ny[i0:i2],color=150
        ;; plots,i0,!y.crange,linestyle=1
        ;; plots,i2,!y.crange,linestyle=1
        ;; if (l eq 50) then stop
     endfor
  endif else begin     
     ;; fill using original data
     ext0 = coro_rad/pixel * scale
     ext1 = coro_rad/pixel * lambda_ref / lambda[nlambda-1] * scale + 1
     mask = fltarr(w,bigh)
     fill_mask = use_mask
     fill_mask[*,0:bigh/2-ext0] = 1.0
     fill_mask[*,bigh/2+ext0:*] = 1.0

     nx = findgen(bigh)
     for l=0,nlambda-1 do begin
        line = reform(synthetic_over[l,*]*fill_mask[l,*])
        fill = filler[l,*]
        mask = reform(fill_mask[l,*])
        
        i0 = min(where(mask eq 0))-1
        i1 = fix(bigh)/2
        i2 = max(where(mask eq 0))+1
        
        line[i0:i2] = fill[i0:i2]
        synthetic_over[l,*] = line
        
        ;; loadct,39,/silent
        ;; plot,nx,line,xs=1,xr=[bigh/2-100,bigh/2+100]
        ;; oplot,nx,fill,color=150
        ;; plots,i0,!y.crange,linestyle=1
        ;; plots,i2,!y.crange,linestyle=1
        ;; ;; if (l eq 50) then stop
        ;; hak
     endfor
  endelse
  
  return,synthetic_over
end
