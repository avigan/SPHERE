; $Id: rescale_padding_1d.pro,v 1.6 2009/05/13 12:19:08 lmugnier Exp $

FUNCTION RESCALE_PADDING_1D,TAB_INPUT = tab_input, $
                            AUTO = auto, ZOOM_IO=zoom_io, FECHS_IO = fechs_io, $
                            KD_IO = kd_io, KF_IO = kf_io, ALT = alt, $
                            VERSION = version, HELP = help
  
;+
; NAME:
;
;	RESCALE_PADDING - Scales an image by zero padding in direct and Fourier spaces.
;
; CATEGORY:
;	
;	Array and Image Processing Routines 
;
; CALLING SEQUENCE: 
;   resampled_image = RESCALE_PADDING(TAB_INPUT = tab_input,
;                                     ZOOM_IO=zoom_io XOR FECHS_IO = fechs_io,
;                                     /AUTO
;                                     [, KD_IO = kd_io, KF_IO = kf_io] [, /ALT]
;                                     [, /VERSION] [, /HELP])
;
; PURPOSE:
;
;   This function scales a square image of even width by adding (or
;   removing) pixels on its sides, in both direct and Fourier space. Both the
;   input image and the scaled one must be correctly (Shannon) sampled for
;   the procedure to work correctly.
;
;   In direct space, in order to avoid any truncation problem it is only
;   possible to add (zero-valued) pixels, not to remove any. In Fourier space
;   addition and removal of pixels are possible depending on the rescaling
;   factor.
;
;   The scaling factor is given by the ratio between the support size in
;   Fourier space (after pixel addition or removal) and the support size (in
;   pixels) in direct space (after pixel addition): if N is the image
;   original size, N' the support size in direct space (after pixels addition)
;   and N" the support size in Fourier space, the scaling factor is
;   ZOOM_IO = N"/N' = 1/FECHS_IO.
;
;   In order to avoid time consuming operations maximal size increase in
;   direct and Fourier space is limited to N'(max)=N"(max)=2N. Minimal support
;   size in Fourier space is limited to N"(min)=N/2. As a consequence scale
;   factor must be no less than 0.5 and no more than 4. After scaling, the
;   scaled image is truncated so that the support size of the returned image
;   is N.
;
;   Note: As null pixels might be added in direct space the image must have
;         zero values close to its edges for the procedure to work correctly.
;
;  KEYWORD PARAMETERS:
;  
;    TAB_INPUT : (input) double or float square image of even width to be
;                rescaled, the image is supposed to be centered between 4
;                pixels (center of the image in [(N-1)/2, (N-1)/2]).
;  
;    ZOOM_IO   : (input for normal usage), float or double, scaling factor of
;                the image physical pixel size. Setting this keyword and /AUTO
;                is the recommended way of calling this procedure, cf below.
;
;    FECHS_IO  : (input for normal usage) kept for compatibility with earlier
;                code. Alternative to setting ZOOM_IO: FECHS_IO = 1/ZOOM_IO
;  
;    /AUTO     : (input) if set, which is recommended, KD_IO and KF_IO are
;                automatically computed from FECHS_IO. Otherwise FECHS_IO and
;                ZOOM_IO are computed from KD_IO and KF_IO.
;  
;    KD_IO     : (optional output for normal usage) positive or null integer,
;                number of pixels to be added on each side of the image in
;                direct space: N' = N + 2*KD_IO. If AUTO is set then KD_IO is
;                computed from N and FECHS_IO and set to the computed value on
;                output.
;  
;    KF_IO     : (optional output for normal usage) signed integer, number of
;                pixels to be added on each side of the image in Fourier
;                space: N" = N + 2*KF_IO. If AUTO is set KF_IO is computed
;                from N and FECHS_IO and set to the computed value on output.
;  
;    /ALT      : (optional input) alternate way of optimizing the choice of
;                KD_IO and KF_IO. Experts only!
;
;    /VERSION  : (optional input) prints version number before execution.
;    
;    /HELP     : (optional input) prints the documentation and exits.
;  
;
; AUTHORS:
;
;   Laurent Mugnier and Nicolas Védrenne.
;
;
; RESTRICTIONS:
;
;   This code is copyright (c) Laurent Mugnier and Nicolas Védrenne, ONERA,
;   2008-2009. 
;
;   The copy of this software that you have received is for personal and for
;   local use only (local meaning by colleagues within your lab/office and
;   thus within your institution). You may use this software for your research
;   needs, but you may not make a commercial use of it. You may not
;   redistribute this software without permission.
;   
; EXAMPLE:
;   Typical call:
;   rescaled_image = RESCALE_PADDING(tab=image, /auto, ZOOM_IO=lambda1/lambda2)
;	For details, see the example in the code after the END of the function.
;
; SEE ALSO:
;   INTERPOLATE (IDL interpolation routines),
;   FFTSHIFT2 (used by this routine).
;   
; HISTORY:
;   $Log: rescale_padding.pro,v $
;   Revision 1.6  2009/05/13 12:19:08  lmugnier
;   Removed last remaining call to rgen and call to abs2 (local library).
;
;   Revision 1.5  2009/05/07 15:18:54  mugnier
;   Added and documented ZOOM_IO keyword to replace FECHS_IO
;   (FECHS_IO is the inverse zoom factor + has a difficult name),
;   Added VERSION and HELP keywords.
;
;   Revision 1.4  2009/05/04 16:40:53  mugnier
;   - Cleanup of code: removed some unused and some temp variables,
;     renamed some variables for clarity, added sanity tests on input,
;     simplified the code a bit, removed dependencies on some local
;     code (rgen, routine_courante),
;   - improved documentation, added quite a few comments in code,
;   - introduced a possible alternate optimization criterion.
;
;   Revision 1.3  2008/09/16 15:16:09  nvedrenn
;   variable names modifications in order to match SPHERE specifications in terms
;   of software furniture
;
;   Revision 1.2  2008/07/15 11:32:31  nvedrenn
;   Examples modifications with more approriate cases
;
;   Revision 1.1  2008/07/10 15:11:19  nvedrenn
;   Initial revision
;
;-

  on_error,2
  
  IF keyword_set(version) THEN $
     printf, -2, '% RESCALE_PADDING: $Revision: 1.6 $, $Date: 2009/05/13 12:19:08 $'
  
  IF (n_params() NE 0) OR keyword_set(help) THEN BEGIN 
     message, 'Help required or incorrect syntax. Documentation:', /INFO
     doc_library, 'rescale_padding'
     retall
  ENDIF

  ss = (size(tab_input))
  dim = ss[0]
  IF (dim NE 1L) THEN message, 'Only 1D arrays are supported'
  NP = ss[1]

  IF (Np MOD 2) THEN $
     message, 'Only arrays of even size are supported'
  
  IF NOT keyword_set(auto) THEN BEGIN
     IF (keyword_set(kf_io) AND keyword_set(kd_io)) THEN BEGIN
        zoom_inside = (double(np)+2D*kf_io)/(double(np)+2D*kd_io)
        fechs_io = 1D/zoom_inside
        zoom_io = zoom_inside
     ENDIF ELSE BEGIN
        message, 'KF_IO and/or KD_IO parameter(s) are missing' 
     ENDELSE
  ENDIF ELSE BEGIN              ; automatic way (recommended)
     IF NOT (keyword_set(fechs_io) XOR keyword_set(zoom_io)) THEN BEGIN 
        message, 'parameter (FECHS_IO XOR ZOOM_IO) is missing' 
     ENDIF ELSE BEGIN
;; We want Fechs *close* *to* N'/N" where N" = N + 2*KF, N' = N + 2*KD
;; => N" = 2*round(N'/(2*Fechs))
;; => KF = (N"-N)/2 = round(N'/(2*Fechs) - N/2) = round(N/2*(1/Fechs-1) +KD/2)
;; We call yy=N/2*(1/Fechs-1) +KD/2
;; We minimize this difference between the `ideal' N" and its closest integer value      
;; Compared to the ALTernate criterion below, this one favors small
;; values of N" ie little truncation in Fourier space.
        
        IF keyword_set(fechs_io) THEN $
           zoom_inside = 1D/fechs_io $
        ELSE $
           zoom_inside = zoom_io
        
        kd_array = lindgen(np/2L+1L) ;from 0 to np/2L=Max Number Of Added Pixels
        yy = double(np)/2D*(zoom_inside-1D)+double(kd_array)*double(zoom_inside) 
        kf_array = round(yy)
        minkf_array = min(abs(yy-kf_array), imin)
        
        kd_io = kd_array[imin]
        kf_io = kf_array[imin]
        
        IF keyword_set(alt) THEN BEGIN ; alternate criterion: minimize error |Fechs-N'/N"|
                                       ; or |zoom-N"/N'| ?
           error = abs(1D/zoom_inside - double(np + 2L*kd_array)/double(np + 2L*kf_array))
;         print, 'Error on N" for chosen kd = ', error[imin], format = '(A, F25.10)'
           minerror = min(error, iminerror)
;          print, 'Minimum error on Fechs of', minerror, ' for i=', $
;              iminerror, format = '(A, F25.10, A, F25.10)'
           kd_io = kd_array[iminerror]
           kf_io = kf_array[iminerror]
        ENDIF
        
     ENDELSE
  ENDELSE
  
  
;;; Extract a part of, or expand, tab_input to NPP pixels:
  npp = np + 2L*kd_io           ; NP Prime (long)
  
  IF npp GT np THEN BEGIN       ; if kd is positive:
     temp = dblarr(npp)
     temp[kd_io] = tab_input
  ENDIF ELSE BEGIN              ; if kd is negative:
     temp = tab_input[-kd_io:-kd_io+npp-1]
  ENDELSE
  
;;; Fourier-transform the result:
  tab_inputf = FFTSHIFT2_1D(temp, /ec)*double(npp)^2
  
;;; Extract a part of or expand the result to NPPP pixels:
  nppp = np + 2L*kf_io          ; NP Prime-Prime (long)
  
  IF nppp GT npp THEN BEGIN
     temp = dcomplexarr(nppp)
     temp[(nppp-npp)/2] = tab_inputf 
  ENDIF ELSE BEGIN              ; kf is smaller than kd, ie there are (kd-kf)
                                ; pixels to be removed on each side of image:
     temp = tab_inputf[kd_io-kf_io:kd_io-kf_io+nppp-1]
  ENDELSE
  
;;; Inverse Fourier-transform the result:
  tab_sc = real_part(FFTSHIFT2_1D(temp, /inv, /ce))/double(nppp)^2
  temp = 0B
  
;;; Extract a part of or expand the result back to NP pixels:
  IF nppp GT np THEN BEGIN      ; kf is positive:
     tab_sctrunc = tab_sc[kf_io:kf_io+np-1]
  ENDIF ELSE BEGIN               ; kf is negative:
     tab_sctrunc = tab_input*0B  ; zero-valued array of same type and size as tab_input
     tab_sctrunc[-kf_io] =tab_sc ; tab_sc is converted to type of tab_sct if needed
  ENDELSE
  
  return, tab_sctrunc
  
END

