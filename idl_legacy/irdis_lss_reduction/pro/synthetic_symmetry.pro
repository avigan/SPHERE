;;
;; synthetic_symmetry
;;

function synthetic_symmetry,sig_over,lambda,use_mask,pla_mask, $
                            SILENT=silent
  nlambda = n_elements(lambda)
  w = (size(sig_over,/dim))[0]
  bigh = (size(sig_over,/dim))[1]

  ;; remove speckles
  if ~keyword_set(silent) then print,'Speckles attenuation with symmetry...'
  synthetic_over = fltarr(w,bigh)

  synthetic_over[*,0:bigh/2-1] = reverse(sig_over[*,bigh/2:*],2)
  synthetic_over[*,bigh/2:*]   = sig_over[*,bigh/2:*]

  if ~keyword_set(silent) then begin
     print,cr()+'  --> done!',format='($,a)'
     print
  endif

  return,synthetic_over
end
