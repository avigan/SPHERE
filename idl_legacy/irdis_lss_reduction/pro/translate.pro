function translate,im,dx,dy,missing=missing

  if (dx eq 0 and dy eq 0) then return,im
  if n_params() lt 3 then dy=0
  
  sz   = size(im)
  dimx = sz[1]
  dimy = sz[2]
  
  xc = (dimx-1)/2.0
  yc = (dimy-1)/2.0
  
  ;;imt=rot(im,0,1.0,xc-dx,yc-dy,missing=missing,cubic=-0.5)
  ;;return,imt
  
  ;;les lignes suivantes font le shift sans interpoler la ou il y
  ;;a des NAN, beaucoup plus rapide
  imt   = make_array(size=sz,value=missing)
  igood = where(finite(im) eq 1,cgood,complement=ibad)

  if cgood lt n_elements(im) then imt[ibad]=!values.f_nan
  if cgood eq 0 then return,imt

  xgood = (igood mod dimx)+round(dx)
  ygood = (igood/dimx)+round(dy)

  imt[xgood,ygood] = interpolate(im,xgood-dx,ygood-dy,cubic=-0.5,missing=missing)
  
  return,imt
end
