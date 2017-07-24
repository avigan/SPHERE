;+
; NAME:
;
;  SPH_IRD_LSS_CLUSTER_BPM
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
;  Returns the bad pixel map for clusters of bad pixels. The file
;  should be stored in the same directory as the routine
;
; CALLING SEQUENCE:
;
;   map = SPH_IRD_LSS_CLUSTER_BPM()
;
; OUTPUTS:
;
;  Bad pixel map in the same format as the IRDIS detector.
;
; MODIFICATION HISTORY:
;
;  28/01/2016 - arthur.vigan@lam.fr
;  Created
;
;-

function sph_ird_lss_cluster_bpm
  path = routine_filepath('sph_ird_lss_cluster_bpm',/either)
  path = file_dirname(path)

  bpm = readfits(path+'/bpm_cluster.fits.gz')

  return, bpm
end
