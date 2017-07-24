;;
;; mplot
;;
  
pro mplot, EPS=eps, OPEN=open, CLOSE=close, XSIZE=xsize, YSIZE=ysize, ZOOM=zoom, $
           THICK=thick, CHARSIZE=charsize, WINDOW=wnum, FILENAME=filename, DISP=disp, $
           PDF=pdf

  if ~keyword_set(eps)   then eps   = (!d.name eq 'PS')
  if ~keyword_set(filename) then filename = 'idl.eps'
  if ~keyword_set(zoom)  then zoom  = 1.0
  if ~keyword_set(xsize) then xsize = 6.4
  if ~keyword_set(ysize) then ysize = 4.8
  if ~keyword_set(bits)  then bits  = 8
  if ~keyword_set(wnum)  then wnum  = 0
  
  if not keyword_set(open) and not keyword_set(close) then $
     message,"You should either open or close the output device"
  
  if keyword_set(open) and keyword_set(close) then $
     message,"You can't both open and close the output device"  

  if keyword_set(eps) then begin
     ;; open EPS
     if keyword_set(open) then begin
        set_plot,'PS'
        device,FILENAME=filename,XSIZE=xsize*zoom,YSIZE=ysize*zoom, $
               /SCALE,/COLOR,/INCHES,/ENCAP,BITS=bits

        path = filename

        if ~keyword_set(thick)    then t = 2.0 else t = thick
        if ~keyword_set(charsize) then c = 1.0 else c = charsize
        !p.thick     = t
        !p.charthick = t
        !p.charsize  = c
        !p.font      = 0
        !x.thick     = t
        !x.charsize  = c
        !y.thick     = t
        !y.charsize  = c
     endif

     ;; close EPS
     if keyword_set(close) then begin
        if keyword_set(disp) or keyword_set(pdf) then file = FSTAT(!D.UNIT)
        
        device,/CLOSE
        set_plot,'X'
        
        t = 1.0
        c = 1.0
        !p.thick     = t
        !p.charthick = t
        !p.charsize  = c
        !p.font      = -1
        !x.thick     = t
        !x.charsize  = c
        !y.thick     = t
        !y.charsize  = c

        if keyword_set(disp) then begin
           if (!version.os_name eq 'Mac OS X') then begin
              cmd = 'open '+file.name
           endif else if (!version.os_name eq 'linux') then begin
              cmd = 'gv '+file.name
           endif
           spawn,cmd
        endif

        if keyword_set(pdf) then begin
           outfile = strmid(file.name,0,strlen(file.name)-3)+'pdf'
           
           spawn,'ps2pdf -dEPSCrop -dAutoFilterColorImages=false '+ $
                 '-dAutoFilterGrayImages=false -dEncodeColorImages=true '+ $
                 '-dEncodeGrayImages=true -dEncodeMonoImages=true '+ $
                 '-dColorImageFilter=/FlateEncode -dGrayImageFilter=/FlateEncode '+ $
                 '-dMonoImageFilter=/FlateEncode -dOptimize=true '+file.name+' '+outfile
        endif
     endif
  endif else begin
     ;; open window
     if keyword_set(open) then begin
        set_plot,'X'
        wopen,wnum,XSIZE=xsize*zoom*100,YSIZE=ysize*zoom*100
        
        t = 1.0
        if ~keyword_set(charsize) then c = 1.0 else c = charsize
        !p.thick     = t
        !p.charthick = t
        !p.charsize  = c
        !p.font      = -1
        !x.thick     = t
        !x.charsize  = c
        !y.thick     = t
        !y.charsize  = c
     endif

     ;; close window
     if keyword_set(close) then begin
        set_plot,'X'

        t = 1.0
        c = 1.0
        !p.thick     = t
        !p.charthick = t
        !p.charsize  = c
        !p.font      = -1
        !x.thick     = t
        !x.charsize  = c
        !y.thick     = t
        !y.charsize  = c
     endif     
  endelse
end
