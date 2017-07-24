;;
;; lss_noise
;;

function lss_noise,tab,lambda,pixel,NLSD=nlsd
  width  = (size(tab,/dim))[0]
  height = (size(tab,/dim))[1]
  noise  = fltarr(width,height)
  l = 7

  if ~keyword_set(nlsd) then nlsd = 2
  
  for x=0,width-1 do begin
     hlsd = ceil(nlsd*lambda[x]*1d-9/8.*180/!pi*2600/pixel/2)
     
     for y=hlsd,height-1-hlsd do begin
        sub = tab[x,y-hlsd:y+hlsd]
        noise[x,y] = stdev(sub)
     endfor
  endfor
  
  return,noise
end
