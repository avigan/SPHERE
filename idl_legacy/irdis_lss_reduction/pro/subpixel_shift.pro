; $Id: subpixel_shift.pro,v 1.2 2009/03/30 10:17:28 lmugnier Exp $

;+
; NAME:
;	SUBPIXEL_SHIFT - shift a 2D image by any non-integer amount via a FFT.
;
; CATEGORY:
;	Array and Image Processing Routines
;
; CALLING SEQUENCE:
;	shifted_image = subpixel_shift(IMAGE, XSHIFT = xshift, YSHIFT = yshift 
;	                               [, /DOUBLE] [, /VERSION] [, /HELP])
;
; PURPOSE:
;   This function returns an image shifted by a non-integer amount via a
;   Fourier domain computation.
;
;   The IMAGE must be square and of even width.
;
; POSITIONAL PARAMETERS:
;	IMAGE: (input) image to be shifted, 2D real (float or double) array.
;
; KEYWORD PARAMETERS:
;	XSHIFT   : (input) amount of desired shift in X direction, float or double.
;
;	YSHIFT   : (input) amount of desired shift in Y direction, float or double.
;
;   /DOUBLE  : (optional input) if this keyword is set, or if the image is in
;              double precision, the computation is performed in double
;              precision, and the result is returned in double precision.
;   
;   /VERSION : (optional input) prints version number before execution.
;   
;   /HELP    : (optional input) prints the documentation and exits.
;
; AUTHORS:
; $Author: lmugnier $
;
; RESTRICTIONS:
;   This code is copyright (c) Laurent MUGNIER, ONERA, 2009. 
;   
; EXAMPLE:
;   image = fltarr(32,32) & image[16,16] = 1.0
;   shifted_image = subpixel_shift(image, xshift=2.0, yshift=1.0) 
;   tvwin, rebin(image+2.0*shifted_image, 256,256,/sample)
;   ; + visually check the [2.0, 1.0] shift in the sum
;   
; SEE ALSO:
;   IDL "shift" routine for integer shifts. 
;   
; HISTORY:
;   The initial revision is simply a robustified version of ONERA's "shift_2d"
;   routine, with English keywords and a routine name not already used by
;   Sphere's IDL software.
;   $Log: subpixel_shift.pro,v $
;   Revision 1.2  2009/03/30 10:17:28  lmugnier
;   HELP keyword now correctly honored.
;
;   Revision 1.1  2009/03/16 19:49:52  lmugnier
;   Initial revision
; 
; -

FUNCTION SUBPIXEL_SHIFT, image, XSHIFT = xshift, YSHIFT = yshift, DOUBLE = double, $
                         VERSION = version, HELP = help

  on_error,2
  IF keyword_set(version) THEN $
     printf, -2, '% SUBPIXEL_SHIFT: $Revision: 1.2 $, $Date: 2009/03/30 10:17:28 $'

  IF ((n_params() NE 1L) OR (n_elements(xshift) NE 1L) OR $
      (n_elements(yshift) NE 1L) OR keyword_set(help)) THEN BEGIN
     message, "Usage : shiftedimage = SUBPIXEL_SHIFT(image, XSHIFT=xshift, " + $
              "YSHIFT=yshift [, /DOUBLE] [, /VERSION] [, /HELP])"
     retall
  ENDIF


  sz = size(image)
  NP = sz[1]                    ; image assumed square NPxNP and NP additionally assumed even
  NL = sz[2]                    ; Number of Lines in the image
  type = sz[sz[0]+1]            ; data type: 4=float, 5=double, 6=complex, 9=double-complex, etc. 

  IF ((type EQ 6) OR (type EQ 9)) THEN $
     message, 'SUBPIXEL_SHIFT only works on real-valued images.'

  IF ((NP NE NL) OR ((NP MOD 2) NE 0)) THEN $ 
     message, 'SUBPIXEL_SHIFT only works on square images of even width.'

  IF NOT (keyword_set(double) OR (type EQ 5)) THEN BEGIN
     x_ramp = (fltarr(NP)+1.)##(findgen(NP)-NP/2) ; ramp along x (2D array)
     y_ramp = (findgen(NP)-NP/2)##(fltarr(NP)+1.) ; ramp along y (2D array)

     tilt = (-2.*!pi/float(NP)) * (float(xshift)*x_ramp+float(yshift)*y_ramp)

     shifted_image = real_part(fft(fft(image, -1) * $
                                   shift(complex(cos(tilt),sin(tilt)), -NP/2, -NP/2), 1))
  ENDIF ELSE BEGIN
     x_ramp = (dblarr(NP)+1D)##(dindgen(NP)-NP/2) ; ramp along x (2D array)
     y_ramp = (dindgen(NP)-NP/2)##(dblarr(NP)+1D) ; ramp along y (2D array)

     tilt = (-2D*!dpi/double(NP)) *(double(xshift)*x_ramp+double(yshift)*y_ramp)

     shifted_image = real_part(fft(fft(image, -1) * $
                                   shift(dcomplex(cos(tilt),sin(tilt)), -NP/2, -NP/2),1))
  ENDELSE

  return, shifted_image

END
