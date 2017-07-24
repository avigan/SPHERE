; $Id: subpixel_shift.pro,v 1.2 2009/03/30 10:17:28 lmugnier Exp $

;+
; NAME:
;	SUBPIXEL_SHIFT - shift a 1D vector by any non-integer amount via a FFT.
;
; CATEGORY:
;	Array and Image Processing Routines
;
; CALLING SEQUENCE:
;	shifted_image = subpixel_shift_1d(VECTOR, PSHIFT = pshift
;	                               [, /DOUBLE] [, /VERSION] [, /HELP])
;
; PURPOSE:
;   This function returns an image shifted by a non-integer amount via a
;   Fourier domain computation.
;
;   The IMAGE must be square and of even width.
;
; POSITIONAL PARAMETERS:
;	VECTOR: (input) vector to be shifted, 1D real (float or double) array.
;
; KEYWORD PARAMETERS:
;	PSHIFT   : (input) amount of desired shift, float or double.
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

FUNCTION SUBPIXEL_SHIFT_1D, vector, PSHIFT = pshift, DOUBLE = double, $
                   VERSION = version, HELP = help

on_error,2
IF keyword_set(version) THEN $
    printf, -2, '% SUBPIXEL_SHIFT: $Revision: 1.2 $, $Date: 2009/03/30 10:17:28 $'

IF ((n_params() NE 1L) OR (n_elements(pshift) NE 1L) OR $
    keyword_set(help)) THEN BEGIN
    message, "Usage : shiftedvector = SUBPIXEL_SHIFT(vector, PSHIFT=pshift, " + $
             "[, /DOUBLE] [, /VERSION] [, /HELP])"
    retall
ENDIF


sz = size(vector)
NP = sz[1]               ; vector of size NP and NP assumed even
type = sz[sz[0]+1]       ; data type: 4=float, 5=double, 6=complex, 9=double-complex, etc. 

IF ((type EQ 6) OR (type EQ 9)) THEN $
    message, 'SUBPIXEL_SHIFT only works on real-valued vectors.'

IF ((NP MOD 2) NE 0) THEN $ 
    message, 'SUBPIXEL_SHIFT only works on vectors of even size.'

IF NOT (keyword_set(double) OR (type EQ 5)) THEN BEGIN
   ramp = (findgen(NP)-NP/2)

   tilt = (-2.*!pi/float(NP)) * (float(pshift)*ramp)
   
   shifted_vector = real_part(fft(fft(vector, -1) * $
                                  shift(complex(cos(tilt),sin(tilt)), -NP/2), 1))
ENDIF ELSE BEGIN
   ramp = (dindgen(NP)-NP/2)
   
   tilt = (-2D*!dpi/double(NP)) *(double(pshift)*ramp)
   
   shifted_vector = real_part(fft(fft(vector, -1) * $
                                  shift(dcomplex(cos(tilt),sin(tilt)), -NP/2),1))
ENDELSE

return, shifted_vector

END
