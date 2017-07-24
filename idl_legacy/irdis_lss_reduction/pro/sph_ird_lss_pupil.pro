;;
;; sph_ird_lss_pupil
;;

function sph_ird_lss_pupil,taille,tailleM1,M1_ONLY=m1_only,ORIENT=orient

  if ~keyword_set(orient) then orient = 0.0
  
  tailleM2 = tailleM1*0.14
  taillesp = tailleM1*0.005 > 1.
  
  cc = (taille-1)/2.
  xx = (findgen(taille) # replicate(1,taille)) - cc
  yy = transpose(xx)
  rho = sqrt(xx^2+yy^2)
  
  pupil = rho le (tailleM1/2.)
  
  if ~keyword_set(M1_only) then begin
     ;; secondary
     pupil -= rho le (tailleM2/2.)

     ;; spiders
     largeur = taillesp

     x1 = (taille-1)/2.
     y1 = (taille-1)/2.+tailleM1/2

     ref = fltarr(taille,taille)
     ref[x1-largeur/2.+1:x1+largeur/2.-1+1, $
         taille/2+tailleM2/(2*sqrt(2)):taille/2+tailleM1/2<taille-1] = 1.
     
     croix1 = rot(ref,-5.5,1,x1,y1,/pivot,/cubic)     
     croix2 = rot(rotate(ref,3),5.5,1,y1,x1,/pivot,/cubic)
     
     croix  = croix1+croix2
     croix  = rot(shift(croix,-1,-1)+rot(croix,180),45+orient,cubic=-0.5,missing=0.0)
     
     pupil = pupil*(1-shift(croix,1,1)) >0 <1
  endif
  
  return,pupil
end
