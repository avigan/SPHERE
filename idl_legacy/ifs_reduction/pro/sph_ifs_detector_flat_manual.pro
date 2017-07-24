;+
; NAME:
;
;  SPH_IFS_DETECTOR_FLAT_MANUAL
;
; PURPOSE:
;
;  Creates a SPHERE/IFS detector flat field.
;
; CALLING SEQUENCE:
;
;  Intended to be called from command line:
;
;    idl -e sph_ifs_preprocess -args ${SOF} ${FFNAME} ${BPNAME}
;
; DESCRIPTION:
;
;  This routine creates the master detector flats for the IFS. It is
;  used as a replacement of the sph_ifs_master_detector_flat recipe
;  from the ESO pipeline.
;
;  It takes as input 2 and only 2 raw flat field frames, and some
;  optional bad pixel map. The raw flat field frames are simply
;  subtracted. The subtracted frame is then normalized to its median
;  value, and bad pixels are corrected using the MASKINTERP() routine
;  and sigma-clipping.
;
;  The routine outputs a flat field frame and a bad pixel map
;  calculated as the pixels departing from more than 10% for the
;  median flat field value (normally equal to 1.0)
;
; INPUTS:
;
;  SOF - path to the SOF file
;
; OPTIONAL INPUTS:
;
;  FFNAME - file name of the output flat field
;
;  BPNAME - file name of the output bad pixel map
;
; OUTPUTS:
;
;  Detector flat field and bad pixel map
;
; EXAMPLE:
;
;  The SOF file should always look like this:
;  
;    /path/to/raw/data/SPHER.2015-01-01T00:00:00.000.fits    IFS_DETECTOR_FLAT_FIELD_RAW
;    /path/to/raw/data/SPHER.2015-01-01T00:01:00.000.fits    IFS_DETECTOR_FLAT_FIELD_RAW
;
;    /path/to/cal/data/badpixel_map.fits                     IFS_STATIC_BADPIXELMAP
;
;  The call sequence would be:
;
;   idl -e sph_ifs_detector_flat_manual -args flat.sof flatfield_filename.fits badpix_filename.fits
;
; MODIFICATION HISTORY:
;
;  arthur.vigan - 08/2015 - commented for distribution
;
;  arthur.vigan - 07/2015 - original version with basic features
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
;   because this code is dedicated to the SPHERE/IFS subsystem, please
;   cite the papers relevant to your observations from the following
;   list:
;
;    * IFS general descripton: Claudi et al., 2008, SPIE, 7014
;    * performance: Mesa et al., 2015, A&A, 576, 121
;
;   And in particular, if you are using this routine:
;
;    * reduction pipeline: Vigan et al., 2015, MNRAS, 454, 129
;
;   We are grateful for your effort, and hope that this tool will
;   contribute to your scientific work and discoveries. Please feel
;   free to report any bug or possible improvement to the author(s)
;   
;-

pro sph_ifs_detector_flat_manual
  cl_args = command_line_args(count=count)
  print,'sph_ifs_detector_flat_manual: ',cl_args
  print
  
  if (count ne 1) and (count ne 3) then begin
     message,string(count,format='("Unexpected number of arguments: ' + $
                    '",I0," instead of 1 or 3")'),/inform
     exit,status=1
  endif

  ;; root = '/Users/avigan/data/SPHERE-GTO/IRDIFS/HR4796_IFS-J/'
  ;; line = root+'sof/flat_white.sof '+root+'ff.fits '+root+'bp.fits'
  ;; cl_args = strsplit(line,' ',/extract)
  ;; count = n_elements(cl_args)
  
  ;;
  ;; arguments
  ;;
  sof = cl_args[0]
  if (count eq 3) then begin
     ffname = cl_args[1]
     bpname = cl_args[2]

     ffname = strmid(ffname,0,strlen(ffname)-5)
     bpname = strmid(bpname,0,strlen(bpname)-5)
  endif else begin
     ffname = 'master_detector_flat'
     bpname = 'dff_badpixels'
  endelse
    
  print,'Producing IFS detector flat field:'
  print,' SOF FILE  = '+sof
  print,' FF FILE   = '+ffname
  print,' BP FILE   = '+bpname  
  print

  ;;
  ;; SOF file
  ;;
  readcol,sof,file,type,format='A,A',count=nfiles,/silent
  if (nfiles eq 0) then begin
     message,'SOF file empty or not formated correctly',/inform
     exit,status=1
  endif

  ff = where(type eq 'IFS_DETECTOR_FLAT_FIELD_RAW',nff)
  bb = where(type eq 'IFS_STATIC_BADPIXELMAP',nbb)

  print,'SOF file content:'
  print,' * ',nff,' flat field file(s)',format='(A0,I3,A0)'
  print,' * ',nbb,' badpix file(s)',format='(A0,I3,A0)'
  print

  ;;
  ;; flat field files
  ;;
  if (nff eq 2) then begin
     ff_file = file[ff]
  endif else begin
     message,'Their should be 2 raw flat field files in the SOF file',/inform
     exit,status=1
  endelse

  ;;
  ;; bad pixel map files
  ;;
  if (nbb ne 0) then begin
     bpm_file = file[bb]
     
     bpm = readfits(bpm_file[0])     
     bpm = fltarr(size(bpm,/dim))
     for f=0,n_elements(bpm_file)-1 do begin
        tmp = readfits(bpm_file[f])
        bpm = bpm or tmp
     endfor
     bpm  = byte(abs(1-bpm))
  endif else bpm = replicate(0,2048,2048)
  
  ;;
  ;; lamps
  ;;
  ff0 = readfits(ff_file[0],hdr0)
  ff1 = readfits(ff_file[1],hdr1)
  
  lamp_switches0 = [sxpar_eso(hdr0,'HIERARCH ESO INS2 LAMP1 ST'), $
                    sxpar_eso(hdr0,'HIERARCH ESO INS2 LAMP2 ST'), $
                    sxpar_eso(hdr0,'HIERARCH ESO INS2 LAMP3 ST'), $
                    sxpar_eso(hdr0,'HIERARCH ESO INS2 LAMP4 ST'), $
                    sxpar_eso(hdr0,'HIERARCH ESO INS2 LAMP5 ST')]
  ilamp0 = where(lamp_switches0 ne 0)+1

  lamp_switches1 = [sxpar_eso(hdr1,'HIERARCH ESO INS2 LAMP1 ST'), $
                    sxpar_eso(hdr1,'HIERARCH ESO INS2 LAMP2 ST'), $
                    sxpar_eso(hdr1,'HIERARCH ESO INS2 LAMP3 ST'), $
                    sxpar_eso(hdr1,'HIERARCH ESO INS2 LAMP4 ST'), $
                    sxpar_eso(hdr1,'HIERARCH ESO INS2 LAMP5 ST')]
  ilamp1 = where(lamp_switches1 ne 0)+1

  if (ilamp0 ne ilamp1) then begin
     message,'Flat fields were not taken with the same lamp',/inform
     exit,status=1
  endif     

  ;;
  ;; create master flat
  ;;
  dit0 = sxpar_eso(hdr0,'HIERARCH ESO DET SEQ1 DIT')
  dit1 = sxpar_eso(hdr1,'HIERARCH ESO DET SEQ1 DIT')
  
  ;; order by increasing DIT
  if (dit1 lt dit0) then begin
     tmp  = dit1
     dit1 = dit0
     dit0 = tmp

     tmp  = ff1
     ff1  = ff0
     ff0  = tmp
  endif
  
  if ((size(ff0))[0] eq 3) then ff0 = median(ff0,dim=3)
  if ((size(ff1))[0] eq 3) then ff1 = median(ff1,dim=3)

  flat_sub = ff1 - ff0

  fake_flat = replicate(1.0,2048,2048)
  fake_dark = replicate(0.0,2048,2048)
  
  nflat = maskinterp(flat_sub,bpm,1,6,'csplinterp',gpix=10,gpoints=5,cdis=2)
  nflat = sigma_filter(nflat,5.,n_sigma=5,n_change=nchange)
  nflat = sigma_filter(nflat,5.,N_sigma=3,radius=3,/iterate,n_change=n_change)

  nflat = nflat / median(nflat)

  fbpm  = (nflat lt 0.9) or (nflat gt 1.1)
  fbpm  = byte(abs(1-fbpm))
  nflat = maskinterp(nflat,fbpm,1,6,'csplinterp',gpix=10,gpoints=5,cdis=2)
  nflat = sigma_filter(nflat,5.,n_sigma=5,n_change=nchange)

  final_flat = nflat / median(nflat)
  final_bpm  = (final_flat lt 0.95) or (final_flat gt 1.1)

  ;;
  ;; field edges
  ;;
  array = shift(dist(10),5,5)
  kern0 = array le 5  
  mask  = dilate(dilate(erode(final_bpm,kern0),kern0),kern0)

  final_flat[where(mask eq 1)] = 1
  
  ;;
  ;; save result
  ;;
  lamp_suffix = '_l'+string(ilamp0,format='(I1)')
  
  writefits,ffname+lamp_suffix+'.fits',final_flat,hdr0
  writefits,bpname+lamp_suffix+'.fits',final_bpm,hdr0
     
  fin:
end
