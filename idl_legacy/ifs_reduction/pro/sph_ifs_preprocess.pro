;+
; NAME:
;
;  SPH_IFS_PREPROCESS
;
; PURPOSE:
;
;  Performs the pre-processing of SPHERE/IFS data.
;
; CALLING SEQUENCE:
;
;  Intended to be called from command line:
;
;    idl -e sph_ifs_preprocess -args ${SOF} ${DO_COLLAPSE}
;                                    ${DO_BKG_SUB} ${DO_BADPIXEL}
;                                    ${DO_CROSSTALK} ${COLLAPSE_TYPE}
;                                    ${COLLAPSE_VAL} ${COLLAPSE_TOL}
;
; DESCRIPTION:
;
;  This routine performs the pre-processing steps for the SPHERE/IFS
;  raw data that are necessary to be done before the creation of the
;  (x,y,lambda) data cubes by the ESO pipeline. The purpose is to
;  improve the results in terms of bad pixels and cross-talk
;  correction.
;
;  The routine is written to be called from the command line, with
;  some mandatory arguments: sof files and control keywords.
;
;  The routine offers several possibilities:
;   - collapse/bin the frames according to different binning schemes
;   - subtract background directly in the raw frames
;   - correct bad pixels directly in the raw frames
;   - correct the spectral crosstalk
;
;  The SOF file is formated in the same way as the ones used by the
;  ESO pipeline. The possible data types are:
;   - IFS_RAW: raw IFS science files (minimum: 1)
;   - IFS_MASTER_DARK: IFS master dark/background frame (maximum: 1)
;   - IFS_STATIC_BADPIXELMAP: bad pixel map(s)
;
;  Then the pre-processing is controled using the input keywords. A
;  value of "0" for the keyword does not execute the correspondong
;  step; any other value executes it, provided that the required files
;  are present in the SOF file (at least one bad pixel map for BP
;  correction, one background for background subtraction).
;
;  The COLLAPSE step (binning) offers a large degree of control
;  through the use of three keywords. There are three types of
;  collapse/binning:
;
;   - MEAN: the whole IFS cube is collapsed as a single frame using a
;           "mean" operation. Useful for a quicklook analysis of the
;           data or when the field rotation per cube is very small;
;
;   - ANGLE: the frames are binned in a way so that a constant field
;            rotation is obtained in the resulting binned
;            frames. Useful to make the companions smearing uniform
;            within a sequence and reduce the biases of the ADI;
;            
;   - COADD: always bins together a fixed number of frames. Useful to
;            decrease the number of frames, while still keeping a finer
;            sampling of the field rotation than the MEAN collapse
;            
;  For the ANGLE collapse, there are two additional required keywords:
;  
;   - COLLAPSE_VAL: the field rotation covered in the binned frame
;                   (e.g. 0.5 deg)
;                   
;   - COLLAPSE_TOL: tolerance on the previous value, because it is
;                   impossible to have exactly the field rotation
;                   given by the COLLAPSE_VAL keyword
;                   
;  Binning/collapse is ALWAYS performed within a data cube file,
;  i.e. frames from different data cube files are NEVER mixed together
;  because of the unknown overheads between two consecutive
;  cubes. Frames that do do satisfy the constraints given by the
;  combination COLLAPSE_VAL/COLLAPSE_TOL are dropped. For example, in
;  a cube of 50 DITs, if 24 frames are required to satisfy
;  COLLAPSE_VAL, the routine will create 2 binned frames out of 48 in
;  the original cube, and it will drop the last 2 frames.
;
;  The bad pixels correction uses the MASKINTERP function with a
;  bicubic spline interpolation. This was found to provide better
;  results than the classical FIXPIX routine, which replaces the bad
;  pixel by the median of surrounding good pixels. A sigma-clipping
;  procedure is used to provide an additional level of correction.
;
;  Spectral crosstalk correction is performed using the
;  SPH_IFS_CROSSTALK routine, which a faster update to the original
;  routine written by Dino Mesa. It gives virtually the same output as
;  the original routine. See documentation for details. By default,
;  only the small-scale crosstalk is corrected, and one would have to
;  manually set remove_large_scale=1 in the call to
;  SPH_IFS_CROSSTALK() to also correct the large-scale crosstalk.
;
;  The final result is saved in two files:
;  
;   - a science file for which the name is identical to the original
;     science file, but with additional suffixes appended that
;     describe how it was processed:
;      * _col+collapse_type when data was collapsed
;      * _bkg when dark/background was subtracted
;      * _bp when bad pixels were corrected
;      * _ct when crosstalk was corrected
;    The original header is written into the file
;    
;  - a "sidecar" FITS table that contains a structure with all relevant
;    informations for each of the frames in the saved science
;    files. These informations are:
;      * file:       file name
;      * img:        DIT number
;      * dit:        DIT value
;      * time:       mean LST time
;      * time_start: start LST time
;      * time_end:   end LST time
;      * ha:         mean hour angle
;      * ha_start:   start hour angle
;      * ha_end:     end hour angle
;      * pa:         mean parallactic angle
;      * pa_start:   start parallactic angle
;      * pa_end:     end parallactic angle
;      * seeing:     seeing (read from the header)
;      * pupoff:     pupil offset (unused at the moment)
;      * filtname:   name of the IFS mode used for the observation
;                    (read from the header)
;      * filtswap:   filter swaping (IRDIS only)
;      * nwave:      number of spectral channels (always 39)
;      * nd:         CPI neutral density (read from the header)
;
; INPUTS:
;
;  SOF - path to the SOF file
;
; OUTPUTS:
;
;  Corrected science file + side car file (see description above)
;
; EXAMPLE:
;
;  If for a set of file we want to do:
;   - collapse using an ANGLE binning of 0.5 deg, with a tolerance of
;     0.05 degree
;   - no dark/background subtraction
;   - bad pixels correction
;   - crosstalk correction
;
;  Then the SOF file would look like this:
;  
;    /path/to/raw/data/SPHER.2015-01-01T00:00:00.000.fits       IFS_RAW
;    /path/to/raw/data/SPHER.2015-01-01T00:01:00.000.fits       IFS_RAW
;    /path/to/raw/data/SPHER.2015-01-01T00:02:00.000.fits       IFS_RAW
;
;    /path/to/cal/data/bad_pixel_map_1.fits      IFS_STATIC_BADPIXELMAP
;    /path/to/cal/data/bad_pixel_map_2.fits      IFS_STATIC_BADPIXELMAP
;
;  And the call sequence would be:
;
;   idl -e sph_ifs_preprocess -args preproc.sof 1 0 1 1 MEAN 0.5 0.05
;
; MODIFICATION HISTORY:
;
;  arthur.vigan - 02/2016 - updated the pupil offset value
;
;  arthur.vigan - 08/2015 - commented for distribution
;
;  arthur.vigan - 07/2015 - added some sigma-clipping to improve bad
;                 pixels correction
;
;  arthur.vigan - 12/2014 - replaced FIXPIX with MASKINTERP for bad
;                 pixels correction. Fixed a bug in the determination
;                 of start/end parallactic angles. Added COADD binning
;                 support
;
;  arthur.vigan - 10/2014 - major update to all components, with
;                 additonal features and bug fixes
;
;  arthur.vigan - 08/2014 - possibility to do background subtraction
;
;  arthur.vigan - 06/2014 - original version with basic features
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

pro sph_ifs_preprocess
  cl_args = command_line_args(count=count)
  print,'sph_ifs_preprocess: ',cl_args
  print
  
  expect = 8
  if (count ne expect) then begin
     message,string(count,expect,format='("Unexpected number of arguments: ' + $
                    '",I0," instead of ",I0)'),/inform
     exit,status=1
  endif
 
  ;; line = '/Users/avigan/Work/ZELDA/validation/2015-12-19/fourier_ifs/sof/preproc.sof ' + $
  ;;        '0 1 1 0 MEAN 0.0 0.0'
  ;; cl_args = strsplit(line,' ',/extract)
  
  ;;
  ;; arguments
  ;;
  sof     = cl_args[0]
  coll    = uint(cl_args[1])
  bkgsub  = uint(cl_args[2])
  bpcor   = uint(cl_args[3])
  xtalk   = uint(cl_args[4])
  colltyp = strlowcase(cl_args[5])
  collval = float(cl_args[6])
  colltol = float(cl_args[7])
  
  print,'Preprocessing of IFS data cubes:'
  print,' SOF FILE  = '+sof
  print,' COLLAPSE  = ',keyword_set(coll),format='(A0,I1)'
  if keyword_set(coll) then begin
     print,' + TYPE      = '+colltyp
     print,' + VALUE     = ',collval,format='(A0,D0.2)'
     print,' + TOLERANCE = ',colltol,format='(A0,D0.2)'
  endif
  print,' BKG SUB   = ',keyword_set(bkgsub),format='(A0,I1)'
  print,' BAD PIX   = ',keyword_set(bpcor),format='(A0,I1)'
  print,' CROSSTALK = ',keyword_set(xtalk),format='(A0,I1)'
  print

  if ~keyword_set(coll) and ~keyword_set(bkgsub) and ~keyword_set(bpcor) and $
     ~keyword_set(xtalk) then begin
     print,' ==> no preprocessing required.'
  endif

  ;;
  ;; SOF file
  ;;
  readcol,sof,file,type,format='A,A',count=nfiles,/silent
  if (nfiles eq 0) then begin
     message,'SOF file empty or not formated correctly',/inform
     exit,status=1
  endif

  ff = where(type eq 'IFS_RAW',nff)
  dd = where(type eq 'IFS_MASTER_DARK',ndd)
  bb = where(type eq 'IFS_STATIC_BADPIXELMAP',nbb)

  print,'SOF file content:'
  print,' * ',nff,' raw file(s)',format='(A0,I3,A0)'
  print,' * ',ndd,' background file(s)',format='(A0,I3,A0)'
  print,' * ',nbb,' badpix file(s)',format='(A0,I3,A0)'
  print

  ;; raw files
  if (nff ne 0) then begin
     raw_file = file[ff]
  endif else begin
     message,'No raw files in SOF file',/inform
     exit,status=0
  endelse

  ;; background
  if (ndd eq 1) then begin
     bkg = readfits(file[dd])
  endif else begin
     if keyword_set(bkgsub) then begin
        message,'More than one background or no background in SOF file. ' + $
                'Skipping background subtraction',/inform
        bkgsub = 0
     endif
  endelse
  
  ;; bad pixel map files
  if (nbb ne 0) then begin
     bpm_file = file[bb]
     
     bpm = readfits(bpm_file[0])     
     bpm = fltarr(size(bpm,/dim))
     for f=0,n_elements(bpm_file)-1 do begin
        tmp = readfits(bpm_file[f])
        bpm = bpm or tmp
     endfor
     bpm  = byte(abs(1-bpm))
  endif else begin
     if keyword_set(bpcor) then begin
        message,'No bad pixel map in SOF file. Skipping BP correction',/inform
        bpcor = 0
     endif
  endelse
  
  ;;
  ;; frames information
  ;;
  print,'Extracting frame information...'
  template_frame = {template_frame_ifs,     $
                    file:'',                $
                    img:-1L,                $
                    dit:-1D,                $
                    time:-1D,               $
                    time_start:-1D,         $
                    time_end:-1D,           $
                    ha:!values.d_nan,       $
                    ha_start:!values.d_nan, $
                    ha_end:!values.d_nan,   $
                    pa:!values.d_nan,       $
                    pa_start:!values.d_nan, $
                    pa_end:!values.d_nan,   $
                    seeing:!values.d_nan,   $
                    pupoff:!values.d_nan,   $
                    filtname:'',            $
                    filtswap:0,             $
                    nwave:39,               $
                    nd:!values.f_nan}
  
  nframes = 0
  for f=0,nff-1 do begin
     hdr = headfits(raw_file[f],/silent)
     ndit = sxpar(hdr,'NAXIS3')
     ndit = (ndit gt 0) ? ndit : 1     
     nframes += ndit
  endfor

  frames = replicate(template_frame,nframes)
  idx = 0L
  for f=0,nff-1 do begin
     hdr = headfits(raw_file[f],/silent)
     
     ;; obs
     DIT  = double(sxpar_eso(hdr,'HIERARCH ESO DET SEQ1 DIT'))
     ndit = sxpar(hdr,'NAXIS3')
     ndit = (ndit eq 0) ? 1 : ndit
     mode = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS2 COMB IFS'),2)

     ;; neutral density
     CPI_ND_hdr  = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS4 FILT2 NAME'),2)
     if (CPI_ND_hdr ne 'OPEN') then CPI_ND = float(strmid(CPI_ND_hdr,3,3)) else CPI_ND = 0.0     
     
     ;; RA,DEC
     ra_drot  = double(sxpar_eso(hdr,'HIERARCH ESO INS4 DROT1 RA'))
     dec_drot = double(sxpar_eso(hdr,'HIERARCH ESO INS4 DROT1 DEC'))
     
     ra0 = fix(ra_drot/10000.)
     ra1 = fix((ra_drot-ra0*10000.)/100)
     ra2 = ra_drot-ra0*10000.-ra1*100.
     ra  = ten([ra0,ra1,ra2])
     
     dec0 = fix(dec_drot/10000.)
     dec1 = fix((dec_drot-dec0*10000.)/100)
     dec2 = dec_drot-dec0*10000.-dec1*100.
     dec  = ten([dec0,dec1,dec2])

     ;; derotator offset
     DROT_mode = strtrim(sxpar_eso(hdr,'HIERARCH ESO INS4 DROT2 MODE'),2)
     if (DROT_mode eq 'ELEV') then pupoff = 135.87D - 100.46D else pupoff = 0D
     
     ;; observatory for PA calculation
     geolat = sxpar_eso(hdr,'HIERARCH ESO TEL GEOLAT')
     geolon = sxpar_eso(hdr,'HIERARCH ESO TEL GEOLON')

     ;; julian date
     obs_in  = sxpar_eso(hdr,'DATE-OBS')
     obs_out = sxpar_eso(hdr,'HIERARCH ESO DET FRAM UTC')
     
     jul_in  = double(date_conv(obs_in, 'J'))
     jul_out = double(date_conv(obs_out,'J'))
     
     delta = (jul_out - jul_in) / ndit

     ;; seeing
     seeing_beg = sxpar_eso(hdr,'HIERARCH ESO TEL AMBI FWHM START')
     seeing_end = sxpar_eso(hdr,'HIERARCH ESO TEL AMBI FWHM END')
     
     for d=0,ndit-1 do begin
        ;; time for each DIT
        time_beg = jul_in + delta * d
        time_mid = jul_in + delta * d + (DIT / 24D / 3600D / 2D)
        time_end = jul_in + delta * d + (DIT / 24D / 3600D)

        ;; seeing
        seeing = interpol([seeing_beg,seeing_end],[time_beg,time_end],time_mid)
        
        ;; parallactic angle     
        ct2lst,lst_beg,geolon,dummy,time_beg
        ha_beg = lst_beg-ra
        pa_beg = parangle(ha_beg,dec,geolat)

        ct2lst,lst_mid,geolon,dummy,time_mid
        ha_mid = lst_mid-ra
        pa_mid = parangle(ha_mid,dec,geolat)     

        ct2lst,lst_end,geolon,dummy,time_end
        ha_end = lst_end-ra
        pa_end = parangle(ha_end,dec,geolat)

        ;; final values
        frames[idx].file       = file_basename(raw_file[f])
        frames[idx].img        = d
        frames[idx].dit        = DIT      ;; DIT value
        frames[idx].pupoff     = pupoff   ;; pupil position correction        
        frames[idx].time       = time_mid ;; local sideral time
        frames[idx].time_start = time_beg ;; local sideral time - start
        frames[idx].time_end   = time_end ;; local sideral time - end
        frames[idx].ha         = ha_mid   ;; hour angle
        frames[idx].ha_start   = ha_beg   ;; hour angle - start
        frames[idx].ha_end     = ha_end   ;; hour angle - end
        frames[idx].pa         = pa_mid   ;; parallactic angle
        frames[idx].pa_start   = pa_beg   ;; parallactic angle - start
        frames[idx].pa_end     = pa_end   ;; parallactic angle - end
        frames[idx].seeing     = seeing   ;; seeing
        frames[idx].filtname   = mode     ;; IFS obs mode
        frames[idx].filtswap   = 0        ;; no swap for IFS
        frames[idx].nd         = CPI_ND   ;; neutral density
        
        idx++
     endfor
  endfor
  print
  
  ;;
  ;; processing
  ;;
  print,'Processing...'
  for f=0,nff-1 do begin
     print,raw_file[f],f+1,nff,format='(" ==> ",A0," (",I0,"/",I0,")")'
     
     ;;
     ;; read data
     ;;
     img = readfits(raw_file[f],hdr)

     ;; frame info
     subframes  = frames[where(frames.file eq file_basename(raw_file[f]))]

     ;;
     ;; background subtraction
     ;;
     if keyword_set(bkgsub) then begin
        print,'   * background subtraction'
        
        dim = size(img)
        if (dim[0] eq 3) then ndit = (size(img))[3] $
        else ndit = 1
        for d=0,ndit-1 do begin
           frame = img[*,*,d]
           frame = frame - bkg
           img[*,*,d] = frame
        endfor
     endif

     ;;
     ;; collapse or binning
     ;;
     if keyword_set(coll) then begin
        if (colltyp eq 'mean') then begin
           print,'   * collapse: mean'

           ;; data
           dim = size(img)
           if (dim[0] eq 3) then begin
              img = total(img,3) / dim[3]
              idx = dim[3]-1
           endif else begin
              idx = 0
           endelse

           ;; ensures that pa values have the same sign
           if ((subframes[0].pa_start/subframes[idx].pa_end) lt 0) then begin
              if (subframes[0].pa_start gt subframes[idx].pa_end) then begin
                 subframes[idx].pa_end += 360D
              endif else begin
                 subframes[0].pa_start += 360D
              endelse
           endif
           
           ;; info
           nsubframes = replicate(template_frame,1)
           nsubframes[0].file       = subframes[0].file
           nsubframes[0].img        = 0
           nsubframes[0].dit        = subframes[0].dit
           nsubframes[0].time       = (subframes[0].time_start + subframes[idx].time_end) / 2D
           nsubframes[0].time_start = subframes[0].time_start
           nsubframes[0].time_end   = subframes[idx].time_end
           nsubframes[0].ha         = (subframes[0].ha_start + subframes[idx].ha_end) / 2D
           nsubframes[0].ha_start   = subframes[0].ha_start
           nsubframes[0].ha_end     = subframes[idx].ha_end
           nsubframes[0].pa         = (subframes[0].pa_start + subframes[idx].pa_end) / 2D
           nsubframes[0].pa_start   = subframes[0].pa_start
           nsubframes[0].pa_end     = subframes[idx].pa_end
           nsubframes[0].pupoff     = subframes[0].pupoff
           nsubframes[0].seeing     = mean(subframes[0:idx].seeing)
           nsubframes[0].filtname   = subframes[0].filtname
           nsubframes[0].filtswap   = subframes[0].filtswap
           nsubframes[0].nd         = subframes[0].nd
           subframes = nsubframes
        endif else if (colltyp eq 'angle') then begin
           print,'   * collapse: angle (',collval,' deg)',format='(A0,D0.2,A0)'

           ;; find frames to be binned
           imin = 0L
           imax = 0L
           nbin = 0L
           val_min = [!values.f_nan]
           idx_min = [-1]
           idx_max = [-1]
           for d=0L,ndit-1 do begin
              if (d lt imin) then continue
              
              delta = subframes[imin].pa_start - subframes[imin:*].pa_end
              min   = min(abs(abs(delta)-collval),imax)
              
              if (min gt colltol) then continue

              val_min = [val_min,min]
              idx_min = [idx_min,imin]
              idx_max = [idx_max,imin+imax]
              
              imin = imin + imax + 1
              nbin++
           endfor

           ;; error if no binning was possible
           if (nbin eq 0) then begin
              message,'Could not find enough frames to bin!',/inform
              exit,status=1              
           endif

           val_min = val_min[1:*]
           idx_min = idx_min[1:*]
           idx_max = idx_max[1:*]

           rbin = max(idx_max)+1
           drop = ndit - max(idx_max)-1

           print,'     + tolerance: ',colltol,' deg, max difference: ',max(val_min),' deg', $
                 format='(A0,D0.2,A0,D0.2,A0)'
           print,'     + ',ndit,' DITs in original cube',format='(A0,I0,A0)'
           print,'     + ',rbin,' DITs binned into ',nbin,format='(A0,I0,A0,I0)'
           print,'     + ',drop,' DITs dropped from original cube',format='(A0,I0,A0)'
           
           ;; bin data
           nimg       = fltarr(2048,2048,nbin)
           nsubframes = replicate(template_frame,nbin)
           for d=0,nbin-1 do begin
              imin = idx_min[d]
              imax = idx_max[d]

              ;; data
              if (imin eq imax) then begin
                 nimg[*,*,d] = img[*,*,imin]
              endif else begin
                 nimg[*,*,d] = total(img[*,*,imin:imax],3) / (imax-imin+1)
              endelse

              ;; ensures that pa values have the same sign
              if ((subframes[imin].time_start/subframes[imax].time_end) lt 0) then begin
                 if (subframes[imin].time_start gt subframes[imax].time_end) then begin
                    subframes[imax].time_end += 360D
                 endif else begin
                    subframes[imin].time_start += 360D
                 endelse
              endif
              
              ;; info
              nsubframes[d].file       = subframes[imin].file
              nsubframes[d].img        = d
              nsubframes[d].dit        = subframes[imin].dit
              nsubframes[d].time       = (subframes[imin].time_start + subframes[imax].time_end) / 2D
              nsubframes[d].time_start = subframes[imin].time_start
              nsubframes[d].time_end   = subframes[imax].time_end
              nsubframes[d].ha         = (subframes[imin].ha_start + subframes[imax].ha_end) / 2D
              nsubframes[d].ha_start   = subframes[imin].ha_start
              nsubframes[d].ha_end     = subframes[imax].ha_end
              nsubframes[d].pa         = (subframes[imin].pa_start + subframes[imax].pa_end) / 2D
              nsubframes[d].pa_start   = subframes[imin].pa_start
              nsubframes[d].pa_end     = subframes[imax].pa_end
              nsubframes[d].pupoff     = subframes[imin].pupoff              
              nsubframes[d].seeing     = mean(subframes[imin:imax].seeing)
              nsubframes[d].filtname   = subframes[imin].filtname
              nsubframes[d].filtswap   = subframes[imin].filtswap
              nsubframes[d].nd         = subframes[imin].nd
           endfor
           img       = nimg
           subframes = nsubframes
        endif else if (colltyp eq 'coadd') then begin
           collval = round(collval)
           print,'   * binning by coadd of ',collval,' consecutive DITs', $
                 format='(A0,I0,A0,D0.2)'
           
           nbin = 0
           idx_min = [-1]
           idx_max = [-1]              
           for d=0,ndit-1,collval do begin
              if (d+collval gt ndit) then continue

              idx_min = [idx_min,d]
              idx_max = [idx_max,d+collval-1]
              
              nbin++
           endfor

           ;; error if no binning was possible
           if (nbin eq 0) then begin
              message,'Could not find enough frames to bin!',/inform
              exit,status=1              
           endif

           idx_min = idx_min[1:*]
           idx_max = idx_max[1:*]

           drop = ndit - nbin*collval
           print,'     + ',ndit,' DITs in original cube',format='(A0,I0,A0)'
           print,'     + ',nbin*collval,' DITs binned into ',nbin,format='(A0,I0,A0,I0)'
           print,'     + ',drop,' DITs dropped from original cube',format='(A0,I0,A0)'

           ;; bin data
           nimg       = fltarr(2048,2048,nbin)
           nsubframes = replicate(template_frame,nbin)
           for d=0,nbin-1 do begin
              imin = idx_min[d]
              imax = idx_max[d]

              ;; data
              if (imin eq imax) then begin
                 nimg[*,*,d] = img[*,*,imin]
              endif else begin
                 nimg[*,*,d] = total(img[*,*,imin:imax],3) / (imax-imin+1)
              endelse

              ;; ensures that pa values have the same sign
              if ((subframes[imin].time_start/subframes[imax].time_end) lt 0) then begin
                 if (subframes[imin].time_start gt subframes[imax].time_end) then begin
                    subframes[imax].time_end += 360D
                 endif else begin
                    subframes[imin].time_start += 360D
                 endelse
              endif
              
              ;; info
              nsubframes[d].file       = subframes[imin].file
              nsubframes[d].img        = d
              nsubframes[d].dit        = subframes[imin].dit
              nsubframes[d].time       = (subframes[imin].time_start + subframes[imax].time_end) / 2D
              nsubframes[d].time_start = subframes[imin].time_start
              nsubframes[d].time_end   = subframes[imax].time_end
              nsubframes[d].ha         = (subframes[imin].ha_start + subframes[imax].ha_end) / 2D
              nsubframes[d].ha_start   = subframes[imin].ha_start
              nsubframes[d].ha_end     = subframes[imax].ha_end
              nsubframes[d].pa         = (subframes[imin].pa_start + subframes[imax].pa_end) / 2D
              nsubframes[d].pa_start   = subframes[imin].pa_start
              nsubframes[d].pa_end     = subframes[imax].pa_end
              nsubframes[d].pupoff     = subframes[imin].pupoff              
              nsubframes[d].seeing     = mean(subframes[imin:imax].seeing)
              nsubframes[d].filtname   = subframes[imin].filtname
              nsubframes[d].filtswap   = subframes[imin].filtswap
              nsubframes[d].nd         = subframes[imin].nd
           endfor
           img       = nimg
           subframes = nsubframes           
        endif else begin
           message,'Unknown collapse type "'+colltyp+'"',/inform
           message,' ==> available types are MEAN, ANGLE and COADD'
        endelse
     endif else begin
        dim = size(img)
        if (dim[0] eq 3) then ndit = (size(img))[3] $
        else ndit = 1
        
        print,'   * no collapse: ',ndit,' DIT(s)',format='(A0,I0,A0)'
     endelse

     ;;
     ;; bad pixel correction
     ;;
     if keyword_set(bpcor) then begin
        print,'   * bad pixels correction'
        
        dim = size(img)
        if (dim[0] eq 3) then ndit = (size(img))[3] $
        else ndit = 1
        for d=0,ndit-1 do begin
           frame = img[*,*,d]
           
           frame = maskinterp(frame,bpm,1,6,'csplinterp',gpix=10,gpoints=4,cdis=2)
           frame = sigma_filter(frame,5.,n_sigma=5,n_change=nchange)
           frame = sigma_filter(frame,5.,N_sigma=3,radius=2,/iterate,n_change=n_change)
           
           img[*,*,d] = frame
        endfor
     endif

     ;;
     ;; crosstalk
     ;;
     if keyword_set(xtalk) then begin
        print,'   * crosstalk correction'

        dim = size(img)
        if (dim[0] eq 3) then ndit = (size(img))[3] $
        else ndit = 1
        for d=0,ndit-1 do begin
           frame = img[*,*,d]

           frame = sph_ifs_crosstalk(frame,remove_large_scale=0)

           img[*,*,d] = frame
        endfor        
     endif

     ;;
     ;; save result
     ;;
     suffix = '_preproc'
     if keyword_set(coll)   then suffix = suffix+'_col_'+colltyp
     if keyword_set(bkgsub) then suffix = suffix+'_bkg'
     if keyword_set(bpcor)  then suffix = suffix+'_bp'
     if keyword_set(xtalk)  then suffix = suffix+'_ct'

     fname = raw_file[f]
     fname = strmid(fname,0,strlen(fname)-5)
     print,'   * saving into '+fname+suffix+'.fits'

     ;; write data and info in the same file
     writefits,fname+suffix+'.fits',img,hdr
     mwrfits,subframes,fname+suffix+'.fits',/silent

     ;; write info in a sidecar file for safety
     mwrfits,subframes,fname+suffix+'_info.fits',/silent,/create
     print
     
  endfor
  
  fin:
end
