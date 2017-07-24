; $Id: fftshift2_1d.pro sauvage Exp $ ;(Mettre en 1ere ligne)
FUNCTION fftshift2_1d, array, direct = direct, fourier = fourier, inverse = $
                       inverse, cc = cc, ee = ee, ce = ce, ec = ec
;+
;NOM :
;  FFTSHIFT2  - calcul de fft avec origine variable
;   
;CATEGORIE :
;   Signal Processing Routines
;
;SYNTAXE :
;   farray = fftshift2(array )
;
;DESCRIPTION : calcul de fft avec comme reference de position le centre du
;              tableau, entre 4 pixels ([(N-1)/2.,(N-1)/2.]) (dans espace direct et Fourier).
;   
;   ARGUMENTS :
;
;    array      : tableau dont on veut calculer la TF
;    direct   : position [x,y] en pixels de l'origine dans l'espace direct
;    fourier   : position [x,y] en pixels de l'origine dans l'espace de fourier
;
;AUTEUR :
;   $Author: sauvage $
;
;
;-

NP = double(n_elements(array))

IF NOT(keyword_set(inverse)) THEN sens = -1 ELSE sens = 1

IF NOT(keyword_set(direct)) THEN direct = (NP-1)/2D
IF NOT(keyword_set(fourier)) THEN fourier = (NP-1)/2D

IF keyword_set(cc) THEN BEGIN
   direct = (NP)/2D
   fourier = (NP)/2D
ENDIF
IF keyword_set(ce) THEN BEGIN
   direct = (NP)/2D
   fourier = (NP-1)/2D
ENDIF
IF keyword_set(ec) THEN BEGIN
   direct = (NP-1)/2D
   fourier = (NP)/2D
ENDIF
IF keyword_set(ee) THEN BEGIN
   direct = (NP-1)/2D
   fourier = (NP-1)/2D
ENDIF


;-------------------------------------------------------------------
; paramètres utiles pour le 
; centrage fin des fft
;-------------------------------------------------------------------
X = dindgen(NP) MOD NP
Y = transpose(X)
j = dcomplex(0, 1)
;-------------------------------------------------------------------

;-------------------------------------------------------------------
; décalage dans fourier par multiplication dans direct et calcul TF
;-------------------------------------------------------------------
farray = fft(array * $
             exp((-sens)*2D*!dpi*j*fourier*X/NP), $
             sens)
;-------------------------------------------------------------------



;-------------------------------------------------------------------
; décalage dans direct par multiplication dans fourier
;-------------------------------------------------------------------
farray *= exp((-sens)*2D*!dpi*j*direct*X/NP) 
;-------------------------------------------------------------------



;-------------------------------------------------------------------
; normalisation par le facteur complexe
;-------------------------------------------------------------------
farray *= exp(sens * (2D*j*!dpi/NP) * total(direct*fourier))
;-------------------------------------------------------------------




return, farray

END 
