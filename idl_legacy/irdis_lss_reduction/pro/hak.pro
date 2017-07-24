pro hak,NoPrompt=noprompt,Main=main
  if keyword_set(main) then on_error,1 else on_error,2
  
  ;; Clear typeahead buffer before input
  while get_kbrd(0) do junk = 1
  
  ;; Get input
  if Keyword_set(noprompt) then begin 
     junk = get_kbrd(1)
  endif else begin
     print, 'Hit any key to continue...'
     junk = get_kbrd(1)
  endelse
  
  ;; Clear typeahead buffer after input
  while get_kbrd(0) do tmp = 1
  
  ;; ESCape character
  if (byte(junk) eq 27) $
  then message,' --> escape key pressed! Stopping...', $
               /noname,/noprefix,/reset
end
 
