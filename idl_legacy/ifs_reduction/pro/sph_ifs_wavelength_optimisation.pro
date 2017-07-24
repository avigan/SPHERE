;+
; NAME:
;
;  SPH_IFS_WAVELENGTH_OPTIMISATION
;
; PURPOSE:
;
;  Wavlength optimisation routine for MPFIT() minimization.
;
; CALLING SEQUENCE:
;
;  Automatically called by MPFIT() during minimization
;
; DESCRIPTION:
; 
;  Interpolate the wavelength law at the position of the laser lines
;  and compares to the known wavelength of the lasers
;
; INPUTS:
;
;  PAR - fitted parameyers
;
;  LAM_LASER - wavelength of the lasers; in nm
;
;  POS_LASER - position of the laser peals in the wavelength vector;
;              in pixel
;
;  REF_IDX - index of the reference wavelength
;
;  NLAMBDA - number of wavelengths
;
;  SCALE_AVERAGE - scaling factor vector
;
; OUTPUTS:
;
;  Maximum of the difference between the know wavelength of the lasers
;  and the wavelength estimated from the wavelength calibration
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


function sph_ifs_wavelength_optimisation,par,lam_laser=lam_laser,pos_laser=pos_laser, $
                                         ref_idx=ref_idx,nlambda=nlambda, $
                                         scale_average=scale_average
  ref_lam = par[0]
  
  idx       = indgen(nlambda)
  lambda    = replicate(ref_lam,nlambda) * scale_average/scale_average[ref_idx]
  lam_peaks = interpol(lambda,idx,pos_laser)
  
  diff = lam_peaks-lam_laser
    
  return,max(abs(diff))
end
