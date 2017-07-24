;;
;; sph_ird_lss_sky_wave_calib
;;

pro sph_ird_lss_sky_wave_calib,lambda_std,std,skyfile,kernel_fwhm=fwhm, $
                               wave_sky_init=wave_sky_init, $
                               wave_std_init=wave_std_init, $
                               smooth=smo,lambda_final=lambda_final,lambda_par=lambda_par
  on_error,1
  
  if ~keyword_set(xrange) then xrange = [900,1800]  
    
  ;;
  ;; sky transmisison file from ESO sky model tool
  ;;    http://www.eso.org/observing/etc/bin/gen/form?INS.MODE=swspectr+INS.NAME=SKYCALC
  ;; must be calculated with dlambda = 0.1 nm
  ;;
  d = mrdfits(skyfile,1)

  ;; convolution kernel
  nlam_k = 5*fwhm/0.1
  lam_k  = (findgen(2*nlam_k+1) - nlam_k)*0.1
  sigma  = fwhm / (2*sqrt(2*alog10(2)))
  val_k  = exp(-0.5*(lam_k/sigma)^2)

  ;; atsmospheric transmission
  lambda_sky = d.lam*1000
  transm_sky = convol(d.trans,val_k,total(val_k))

  ;; scale atmospheric spectrum
  std = std / max(std)
  if keyword_set(smo) then std = smooth(std,smo)
    
  ii = where((xrange[0] le lambda_std) and (lambda_std le xrange[1]))
  jj = where((xrange[0] le lambda_sky) and (lambda_sky le xrange[1]))

  transm_sky += mean(std[ii])-mean(transm_sky[jj])

  ;; ---------------------------------
  ;; review selected lines
  ;;
  wopen,0,xs=1200,ys=800
  loadct,39,/silent
  !p.multi = [0,1,2]
  
  nwave = n_elements(wave_sky_init)
  for w=0,nwave-1 do begin
     ;; global view
     plot,lambda_std,std,xs=1,xr=[900,1850],ys=1,yr=[0,1],/nodata,xticks=19,xminor=5
     plots,wave_sky_init[w],!y.crange,linestyle=1
     oplot,lambda_std,std
     oplot,lambda_sky,transm_sky,color=250

     for i=0,nwave-1 do begin
        istd = min(where(lambda_std ge wave_std_init[i]))
        plots,lambda_std[istd],std[istd],psym=1
        
        isky  = min(where(lambda_sky ge wave_sky_init[i]))
        plots,lambda_sky[isky],transm_sky[isky],psym=1,color=250
     endfor     
     
     ;; line-by-line view
     xrange = wave_sky_init[w] + [-50,50]     
     plot,lambda_std,std,xs=1,xr=xrange,ys=1,yr=[0,1]
     oplot,lambda_sky,transm_sky,color=250
     
     for i=0,nwave-1 do begin
        istd = min(where(lambda_std ge wave_std_init[i]))
        plots,lambda_std[istd],std[istd],psym=1
        
        isky  = min(where(lambda_sky ge wave_sky_init[i]))
        plots,lambda_sky[isky],transm_sky[isky],psym=1,color=250
     endfor
     hak,/main
  endfor

  ;; ---------------------------------
  ;; fit lines
  ;;
  !p.multi = [0,2,1]
  pix = findgen(n_elements(lambda_std))
  pix_std_final  = fltarr(nwave)
  wave_sky_final = fltarr(nwave)
  for i=0,nwave-1 do begin
     ext = 4. ;; nm
     
     istd = where(((wave_std_init[i]-ext) le lambda_std) and $
                  (lambda_std le (wave_std_init[i]+ext)))
     sub_l = pix[istd]
     sub_t = std[istd]
     
     fit = mpfitpeak(sub_l,sub_t,par,nterms=5)
     
     pix_std_final[i] = par[1]
     
     plot,sub_l,sub_t,/yno,xtitle='Position [pix]',title='Std',/nodata
     oplot,sub_l,sub_t,thick=!p.thick*2
     oplot,sub_l,fit,color=250,thick=!p.thick*2,linestyle=2
     
     isky  = where(((wave_sky_init[i]-ext) le lambda_sky) and $
                   (lambda_sky le (wave_sky_init[i]+ext)))
     sub_l = lambda_sky[isky]
     sub_t = transm_sky[isky]

     fit = mpfitpeak(sub_l,sub_t,par,nterms=5)

     wave_sky_final[i] = par[1]
     
     plot,sub_l,sub_t,/yno,xtitle='Wavelength [nm]',title='Sky simulation',/nodata
     oplot,sub_l,sub_t,thick=!p.thick*2
     oplot,sub_l,fit,color=250,thick=!p.thick*2,linestyle=2

     al_legend,['data','fit'],linestyle=[0,2],color=[255-!p.background,250],/top,/left,box=0
     
     hak,/main
  endfor
  !p.multi = 0

  fit = poly_fit(pix_std_final,wave_sky_final,1)
  lambda_final = fit[0]+pix*fit[1]
  lambda_par   = fit
  
  print,'residuals:',wave_sky_final-(fit[0]+pix_std_final*fit[1])   
  
  plot,pix_std_final,wave_sky_final,psym=1,/yno
  oplot,pix,lambda_final,linestyle=1,color=250

end
