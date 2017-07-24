;+
; NAME:
;
;  SPH_IFS_CROSSTALK
;
; PURPOSE:
;
;  Corrects the spectral crosstalk in SPHERE/IFS images.
;
; CALLING SEQUENCE:
;
;  NIMG = SPH_IFS_CROSSTALK(IMG)
;
; DESCRIPTION:
;
;  This routines corrects for the SPHERE/IFS spectral crosstalk at
;  small scales and (optionally) at large scales. This correction is
;  necessary to correct the signal that is "leaking" between
;  lenslets. See Antichi et al. (2009ApJ...695.1042A) for a
;  theoretical description of the IFS crosstalk. Some informations
;  regarding its correction are provided in Vigan et al. (2015), but
;  this procedure still lacks a rigorous description and performance
;  analysis.
;  
; INPUTS:
;
;  IMG - input raw SPHERE/IFS image
;
; KEYWORDS:
;
;  REMOVE_LARGE_SCALE - remove large scale crosstalk
;
; OUTPUTS:
;
;  Input image corrected from the spectral crosstalk
;
; SIDE EFFECTS:
;
;  Since the correction of the crosstalk involves a convolution by a
;  kernel of size 41x41, the values at the edges of the frame depend
;  on how you choose to apply the convolution. Current implementation
;  is EDGE_TRUNCATE. In other parts of the image (i.e. far from the
;  edges), the result is identical to original routine by Dino
;  Mesa. Note that in the original routine, the convolution that was
;  coded did not treat the edges in a clean way defined
;  mathematically. The IDL CONVOL() routine offers different
;  possibilities for the edges that are all documented.
;
; MODIFICATION HISTORY:
;
;   arthur.vigan - 08/2015 - commented for distribution
;
;   arthur.vigan - 07/2015 - option to enable large scale crosstalk
;                            subtraction 
;
;   arthur.vigan - 10/2014 - bug fix for the histogram part
;
;   arthur.vigan - 10/2014 - major optimization and simplification.
;                            Speed improvement byfactor >10
;
;   dino.mesa    - 00/0000 - original version sent to arthur.vigan for
;                            analysis of IFS data
;
; AUTHORS:
;
;   Original version:
;    Dino Mesa
;    INAF/Osservatorio Astronomico di Padova
;    dino.mesa@oapd.inaf.it
;
;   Current version:
;    Arthur Vigan
;    Laboratoire d'Astrophysique de Marseille
;    arthur.vigan@lam.fr
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


function sph_ifs_crosstalk,img,remove_large_scale=remove_large_scale
  dimimg = 2048
  sepmax = 20   ;; definition of the dimension of the matrix

  dimmat = sepmax*2+1
  matsub = dblarr(dimmat,dimmat)
  pc = sepmax
  bfac = 0.727986/1.8

  ;; all calculations done in double
  img = double(img)
  
  ;; defines a matrix to be used around each pixel
  ;; (the value of the matrix is lower for greater
  ;; distances form the center.
  for k=0,dimmat-1 do begin
     for j=0,dimmat-1 do begin
        if abs(pc-k) gt 1 or abs(pc-j) gt 1 then begin
           rdist = sqrt((float(k)-float(pc))^2+(float(j)-float(pc))^2)
           matsub[k,j] = 1./(1+rdist^3/bfac^3)
        endif
     endfor
  endfor

  ;; convolution and subtraction
  conv = convol(img,matsub,/edge_truncate)
  imgsub = img - conv

  ;; remove large scale crosstalk
  if keyword_set(remove_large_scale) then begin
     ;; copy
     imgfin = imgsub
     
     ;; Step 2 ==> calculation of the
     ;; cross-talk on large scale
     ;; on sub-images of 64x64 pixels  
     img = imgfin
     dimimg = 2048
     dimsub = 64
     dimimgct = dimimg/dimsub
     valinix = intarr(dimimgct*dimimgct)
     valfinx = intarr(dimimgct*dimimgct)
     valiniy = intarr(dimimgct*dimimgct)
     valfiny = intarr(dimimgct*dimimgct)
     konta = 0
     
     ;; defines the positions of the subimages
     for j=0,dimimgct-1 do begin
        for k=0,dimimgct-1 do begin
           valinix[konta] = k*dimsub
           valfinx[konta] = valinix[konta]+dimsub-1
           valiniy[konta] = j*dimsub
           valfiny[konta] = valiniy[konta]+dimsub-1
           konta = konta+1
        endfor
     endfor

     mdnimg = median(img)
     for k=0,dimimg-1 do begin
        for j=0,dimimg-1 do begin
           if abs(img[k,j]) gt 30000. then begin
              img[k,j] = mdnimg  ;; still bad pixels. The bad pixel procedure did not work perfectly.
              ;; print, 'Problemi immagine '+strtrim(numf,1)
           endif
        endfor
     endfor
     
     ;; For each subimage it creates an histogram and defines
     ;; the value of the maximum of the pixel counts distribution
     stephist = 10D
     imgct = dblarr(dimimgct,dimimgct)
     for k=0,n_elements(valinix)-1 do begin
        imgsub  = img[valinix[k]:valfinx[k],valiniy[k]:valfiny[k]]

        ncount = histogram(imgsub,binsize=stephist,loc=vcount,/nan)
        vcount += stephist/2D
        
        rs = where(ncount eq max(ncount))
        valct = vcount(rs(0))
        cy = fix(float(k)/float(dimimgct))
        cx = k-cy*dimimgct
        imgct[cx,cy] = valct
     endfor

     ;; ###########################################################################################
     ;; Step 3 ==> subtraction of the large
     ;; scale cross-talk
     img = imgfin
     dimimg = 2048
     dimsub = 64
     dimimgct = dimimg/dimsub
     valinix = intarr(dimimgct*dimimgct)
     valfinx = intarr(dimimgct*dimimgct)
     valiniy = intarr(dimimgct*dimimgct)
     valfiny = intarr(dimimgct*dimimgct)
     konta = 0

     for j=0,dimimgct-1 do begin
        for k=0, dimimgct-1 do begin
           valinix[konta] = k*dimsub
           valfinx[konta] = valinix[konta]+dimsub-1
           valiniy[konta] = j*dimsub
           valfiny[konta] = valiniy[konta]+dimsub-1
           konta = konta+1
        endfor
     endfor
     
     imgsub = dblarr(dimimg,dimimg)
     for k=0,n_elements(valinix)-1 do begin
        cy = fix(float(k)/float(dimimgct))
        cx = k-cy*dimimgct
        imgsub[valinix[k]:valfinx[k],valiniy[k]:valfiny[k]] = $
           img[valinix[k]:valfinx[k],valiniy[k]:valfiny[k]] - imgct[cx,cy]
     endfor

     imgsub = float(imgsub)
  endif else message,'Skipping large scale crosstalk subtraction',/inform
  
  ;; final output
  return,imgsub
end
