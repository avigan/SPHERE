;+
; NAME:
;
;   SPH_IRD_TRANSMISSION
;
; PURPOSE:
;
;   Provides the neutral density value in a given CPI ND + IRDIS
;   filter setup
;
; CALLING SEQUENCE:
;
;   ND_value  = SPH_IRD_TRANSMISSION( BB_filter, DB_filter, ND_filter,
;                                     [ROOT=... ])
;
; DESCRIPTION:
;
;   The function provides the neutral density value as a function of
;   the CPI ND filter and the IRDIS filter combination. It works for
;   all IRDIS broad-band (BB), dual-band (DB) and narrow-band (ND)
;   filters. The list of available filters is below.
;
;   Broad-band filters:
;    - B_Y
;    - B_J
;    - B_H
;    - B_ND-H
;    - B_Ks
;
;   Dual-band filters:
;    - D_Y23
;    - D_J23
;    - D_H23
;    - D_H32
;    - D_H34
;    - D_K12
;
;   Narrow-band filters:
;    - N_BrG
;    - N_CO
;    - N_CntH
;    - N_CntJ
;    - N_CntK1
;    - N_CntK2
;    - N_FeII
;    - N_H2
;    - N_Hel
;    - N_PaB
;
;   Neutral density filters:
;    - OPEN
;    - ND_1.0
;    - ND_2.0
;    - ND_3.5
;
;   The ND value is calculated as the ratio between the transmission
;   of the filters with and without multiplication by the transmission
;   curve of the CPI ND filter. For the B_ND-H filter, the IRDIS ND
;   filter is also taken into account. The transmission of each filter
;   is read from text files stored on disk.
;
;   The input values for the function can be obtained directly from
;   the raw files FITS headers using the following keywords:
;
;    HIERARCH ESO INS1 FILT NAME  ==> BB_filter
;    HIERARCH ESO INS1 OPTI2 NAME ==> DB_filter
;    HIERARCH ESO INS4 FILT2 NAME ==> ND_filter
;
;   The returned value is in units of neutral density. To convert it
;   to the proper multiplicative factor, one has to apply 10^ND.
;
; INPUTS:
;
;   BB_FILTER - broad-band or narrow-band filter name from the list
;               above, or extracted from the header of raw IRDIS data;
;               format: text
;
;   DB_FILTER - dual-band filter name from the list above, or
;               extracted from the header of raw IRDIS data. When no
;               DB filter is used, the value 'CLEAR' must be passed;
;               format: text
;
;   ND_FILTER - CPI neutral density filter name from the list above,
;               or extracted from the header of raw IRDIS data;
;               format: text
;
; KEYWORD PARAMETERS:
;
;   ROOT - root directory where the filter transmission files are
;          stored. The default value can be modified at the
;          beginning of the programe; format: text
;
; OUTPUTS:
;
;   Neutral density value in the filter(s). In any case, two values
;   are returned, corresponding to the IRDIS left and right fields. In
;   the case of BB or NB filters the two values are identical, while
;   for DB filters the values correspond to each of the DB filters.
;
; EXAMPLE:
;
;   Example use for dual-band imaging (DBI):
;
;   IDL> PRINT, SPH_IRD_TRANSMISSION('B_H','D_H23','ND_3.5')
;         3.28147      2.97996
;
;   In this case, the correction factors to apply to IRDIS left and
;   right images would:
;
;   Image_left_corrected  = Image_left_raw  * 10^(3.28147)
;   Image_right_corrected = Image_right_raw * 10^(2.97996)
;
;   Example use for broad-band imaging (also known as "classical
;   imaging", CI, in the documentation):
;
;   IDL> PRINT, SPH_IRD_TRANSMISSION('B_H','CLEAR','ND_3.5')
;         3.07009      3.07009
;
; MODIFICATION HISTORY:
;
;   arthur.vigan - 07/2015 - additonnal documentation for distribution
;                            ND filter specified by its name rather
;                            than value
;                            added some safeguards
;
;   arthur.vigan - 03/2014 - added basic support for DPI data (i.e. no
;                            crash when polarizers are used in IRDIS)
;
;   arthur.vigan - 01/2013 - added support for B_ND-H filter
;
;   arthur.vigan - 11/2012 - original version received from A.
;                            Boccaletti. Improved to be more robust
;                            and support all filters
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
;   because this code is dedicated to the SPHERE/IRDIS subsystem,
;   please cite the papers relevant to your observations from the
;   following list:
;
;    * IRDIS general descripton: Dohlen et al., 2008, SPIE, 2014
;    * Dual-Band Imaging mode: Vigan et al., 2010, MNRAS, 407, 71
;    * Dual-Polarization Imaging mode: Langlois et al., 2014, SPIE, 9147
;
;   We are grateful for your effort, and hope that this tool will
;   contribute to your scientific work and discoveries. Please feel
;   free to report any bug or possible improvement to the author(s)
;
;-

function sph_ird_transmission, BB_filter, DB_filter, ND_filter, ROOT=root
  on_error, 2

  if ~keyword_set(root) then root = '~/IDL/work/Data/Filters/'

  ;; removing spaces and using uppercase
  BB_filter = strupcase(strtrim(BB_filter,2))
  DB_filter = strupcase(strtrim(DB_filter,2))
  ND_filter = strupcase(strtrim(ND_filter,2))
  
  ;;
  ;; BB filter
  ;;
  test = file_test(root+'SPHERE_IRDIS_'+BB_filter+'.txt')
  if (test eq 0) then message,'File not available. Either the ROOT variable is ' + $
                              'not set properly or '+BB_filter+' is not a valid filter'
  readcol,root+'SPHERE_IRDIS_'+BB_filter+'.txt',w_bb,t_bb,/silent

  ;;
  ;; DBI filter
  ;;
  if ((DB_filter eq 'CLEAR') or $
      (DB_filter eq 'P0-90') or (DB_filter eq 'P45-135')) then begin
     w_db = w_BB
     t_db_0 = replicate(1.,n_elements(w_db))
     t_db_1 = replicate(1.,n_elements(w_db))
  endif else begin
     test = file_test(root+'SPHERE_IRDIS_'+DB_filter+'.txt')
     if (test eq 0) then message,'File not available. Either the ROOT variable is ' + $
                                 'not set properly or '+DB_filter+' is not a valid filter'
     
     readcol,root+'SPHERE_IRDIS_'+DB_filter+'.txt',w_db,t_db_0,t_db_1,/silent
  endelse

  ;;
  ;; ND CPI
  ;;
  readcol,root+'SPHERE_CPI_ND.txt',w_cpi_nd,t_cpi_nd_0,t_cpi_nd_1, $
          t_cpi_nd_2,t_cpi_nd_35,/silent
  case ND_filter of
     'OPEN':   t_cpi_nd = t_cpi_nd_0
     'ND_1.0': t_cpi_nd = t_cpi_nd_1
     'ND_2.0': t_cpi_nd = t_cpi_nd_2
     'ND_3.5': t_cpi_nd = t_cpi_nd_35
     else: message,ND_filter+' is not available.'
  endcase

  ;;
  ;; ND IRDIS
  ;;
  if (BB_filter eq 'B_ND-H') then begin
     readcol,root+'SPHERE_IRDIS_ND.txt',w_ird_nd,t_ird_nd,/silent     
  endif else begin
     w_ird_nd = w_bb
     t_ird_nd = replicate(1.0,n_elements(w_ird_nd))
  endelse

  ;;
  ;; interpolation
  ;;
  lambdainterp    = findgen(2500.-900.+1)+900.
  t_bb_interp     = interpol(t_bb >0 ,w_bb,lambdainterp)
  t_db_0_interp   = interpol(t_db_0 >0 ,w_db,lambdainterp)
  t_db_1_interp   = interpol(t_db_1 >0 ,w_db,lambdainterp)
  t_cpi_nd_interp = interpol(t_cpi_nd >0 ,w_cpi_nd,lambdainterp)
  t_ird_nd_interp = interpol(t_ird_nd >0 ,w_ird_nd,lambdainterp)

  ;;
  ;; integrated ND value
  ;;
  t_final_0 = total(t_bb_interp * t_db_0_interp)
  t_final_1 = total(t_bb_interp * t_db_1_interp)
  t_nd_final_0 = total(t_bb_interp * t_db_0_interp * t_cpi_nd_interp * t_ird_nd_interp)
  t_nd_final_1 = total(t_bb_interp * t_db_1_interp * t_cpi_nd_interp * t_ird_nd_interp)

  NDvalue = -alog10([t_nd_final_0 / t_final_0,t_nd_final_1 / t_final_1])
  
  return,NDvalue
end
