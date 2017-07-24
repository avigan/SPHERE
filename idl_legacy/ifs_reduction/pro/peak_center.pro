;+
; NAME:
;
;  PEAK_CENTER
;
; PURPOSE:
;
;  Finds the center of a peak inside the image
;
; CALLING SEQUENCE:
;
;   PEAK_CENTER, IMG [, EXTENSION=, HIDE_EDGES=, MASK=, SMOOTH=,
;                       PARAM=, FIT=, DISP=, CREF= ]
;
; DESCRIPTION:
; 
;  Find the center of peak inside the image using MPFIT2DPEAK(). The
;  basic use would be just to pass the image to the function, but it
;  can fail if there are several bright objects in the image. Some
;  parameters can be used to improve the fitting.
;
; INPUTS:
;
;  IMG - the image where we want to find the peak
;
; OPTIONAL INPUTS:
;
;  EXTENSION - half-size of the box where we are going to fit the
;              peak. Default in 10 (i.e. a 21x21 box for the fit);
;              in pixel
;
;  HIDE_EDGES - hide the edges of the images by a certain number of
;               pixels; in pixel
;
;  MASK - provide a mask that tells which pixels to consider for the
;         fit; array of the same size as the image
;
;  SMOOTH - smooth the image by a certain number of pixels; in pixel
;
;  DISP - display the fit
;
;  CREF - estimate of the approximate center to use as a starting
;         point; [x,y] values in pixel
;
; OUTPUTS:
;
;  The coordinates of the center within the image
;
; OPTIONAL OUTPUTS
;
;  FIT - fit output from MPFIT
;
;  PARAM - fit parameters from MPFIT
;
; MODIFICATION HISTORY:
;
;  arthur.vigan - 08/2015 - public version of personnal tool
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

function peak_center,img,EXTENSION=ext,HIDE_EDGES=hide,MASK=mask, $
                     SMOOTH=smooth,PARAM=param,FIT=fit,DISP=disp, $
                     CREF=cref
  if ~keyword_set(ext) then ext = 10

  wimg = img
  
  ;; hide edges
  if keyword_set(hide) then begin
     dim = (size(wimg,/dim))
     wimg[0:hide-1,*] = 0.0
     wimg[*,0:hide-1] = 0.0
     wimg[dim[0]-hide-1:*,*] = 0.0
     wimg[*,dim[1]-hide-1:*] = 0.0
  endif

  ;; mask NaN
  bad = where(finite(wimg) eq 0,nbad)
  if (nbad ne 0) then wimg[bad] = 0.0
  
  ;; apply mask
  if keyword_set(mask) then wimg = wimg*mask

  ;; find peak
  if keyword_set(smooth) then tmp = median(wimg,5) $
  else tmp = wimg

  if (keyword_set(cref) and (n_elements(cref) eq 2)) then begin
     cx_int = fix(cref[0])
     cy_int = fix(cref[1])
  endif else begin
     max = max(tmp,imax)
     pos = array_indices(wimg,imax)
     cx_int = pos[0]
     cy_int = pos[1]
  endelse
  
  sub = wimg[cx_int-ext:cx_int+ext,cy_int-ext:cy_int+ext,0]
  fit = mpfit2dpeak(sub,param,TILT=1)
  cx  = param[4] - ext + cx_int
  cy  = param[5] - ext + cy_int

  ;; display
  if keyword_set(disp) then begin
     ds9
     !v->im,wimg
     !v->box,cx_int,cy_int,2*ext+1,2*ext+1,0,/fixed,color='red',thick=2
     !v->point,cx,cy,/fixed,color='green',thick=2
  endif
  
  return,[cx,cy]
end
