;;
;; wopen
;; open a specific window if it is
;; closed, otherwise select it
;;

pro wopen,win,XSIZE=xs,YSIZE=ys,STRICT=strict,_EXTRA=extra
  on_error,2

  ;; work only with X display
  if (!d.name eq 'X') then begin
     ;; get windows status
     device,window_state=window_state
     
     if (n_elements(win) gt 1) then $
        message,'Argument must be a scalar.'

     ;; window number too large
     if (win gt 31) then $
        message,"Can't open more than 32 windows manually."

     ;; open window if needed
     if (window_state[win] eq 0) then begin
        window,win,xsize=xs,ysize=ys,_extra=extra
     endif else begin
        wset,win
        if keyword_set(strict) and (keyword_set(xs) or keyword_set(ys)) then begin
           if (xs ne !d.x_size) or (ys ne !d.y_size) $
           then window,win,xsize=xs,ysize=ys,_EXTRA=extra
        endif
     endelse
  endif
end
