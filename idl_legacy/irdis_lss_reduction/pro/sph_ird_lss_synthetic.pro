;+
; NAME:
;
;  SPH_IRD_LSS_SYNTHETIC
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
;  Performs speckle estimation to create a "synthetic" spectrum
;
; CALLING SEQUENCE:
;
;   synth = SPH_IRD_LSS_SYNTHETIC(func,sig_over,lambda,use_mask,pla_mask, $
;                                 FILL=fill,_EXTRA=ex)
;
; DESCRIPTION:
;
;  This function is a wrapper to call specific functions designed to
;  estimate the speckles in the rescaled spectrum. These functions
;  will create a synthetic spectrum that represent the speckles and
;  stellar halo. This synthetic spectrum can then be subtracted to the
;  spatially rescaled spectrum or be rescaled back and subtracted to
;  the original spectrum.
;
;  A useful option is the ability to fill the area below the
;  coronagraph with some continuous signal that will not create
;  artifacts when rescaled with FFT. It uses the binary masks created
;  by the SPH_IRD_LSS_RESCALE function.
;
; INPUTS:
;
;  FUNC - string containing the name of the function that will be
;         called to perform the speckle estimate. The predefined
;         parameters for such a function are the following:
;
;    SIG_OVER - 2D rescaled (and sometimes oversampled) spectrum
;
;    LAMBDA - wavelength vector; in nanometers
;
;    PLANET_MASK - binary mask for the planet, as created by the
;                  SPH_IRD_LSS_RESCALE function
;
;    USE_MASK - binary mask of the "usable" area, as created by the
;               SPH_IRD_LSS_RESCALE function
;
;    _EXTRA - extra keyword parameters to be passed to the function
;             
;  FILL - switch indicating wether the area below the coronagraph must
;         be filled with continuous data or not
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
;
;  24/07/2012 - arthur.vigan@lam.fr
;  Created
;
;-

function sph_ird_lss_synthetic,func, $
                               sig_over,lambda,use_mask,pla_mask,_EXTRA=ex, $
                               FILL=fill

  synthetic_over = call_function(func,sig_over,lambda,use_mask,pla_mask,_EXTRA=ex)
  
  if keyword_set(fill) then begin
     synthetic_over = sph_ird_lss_fill(synthetic_over,lambda,use_mask,_EXTRA=ex, $
                                       filler=sig_over)
  endif
     
  return,synthetic_over
end
