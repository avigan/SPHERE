pro fixpix,imgs,badpix,outimgs, $
           npix=npix,weight=weight,search_box=search_box, $
           noise=noise,sigma=sigma,dc=dc, $
           silent=silent,badvalmask=badvalmask,NaN=NaN,quick=quick
  
;+
; NAME:
; 		fixpix
; PURPOSE:
; given a image or stack of images and a bad pixel
; mask, will fill in bad pixels by finding the NPIX nearest
; good pixels, toss the highest and lowest ones of the bunch,
; and then arithmatically average. 
;
; NOTES:
;
;		pro fixpix,imgs,badpix,outimgs, $
;		   npix=npix,weight=weight, $
;		   noise=noise,sigma=sigma,dc=dc, $
;		   silent=silent,badvalmask=badvalmask,NaN=NaN
;		   
; bad pixel list is processed in order array is stored in
; memory, i.e. row by row, and is defined by:
;	 0 = bad pixel
;    not 0 = good pixel
;
; If /badval is set, badpix is ignored and the procedure determines
; 	a new bad pixel list by examining the image for pixels set to BADVAL.
; Alternatively, if /NaN is set, any IEEE Not-a-Number values are 
;   considered to be bad pixels and fixed.
;
; NPIX = # of adjacent pixels used for correction (default = 8)
; /weight: weight adjacent pixels by inverse of their distances
;	   in averaging process
;
; checks to see if some pixels are equidistant; if so,
; use all of them, even if the total # of pixels then exceeds NPIX
;
; WARNINGS:
;  - will work for entire bad columns/rows, but 
;    may not be the most sensible thing to use
;  - do not pass it the same array for input & output as this
;    might corrupt the correction algorithm
;
; 7/3/95 MCL
;
; added /noise: replace bad pixels with random noise with
; user defined sigma and dc level                9/24/95 MCL
; badpix can now be a 2d or 3d array             4/18/96 MCL
; uses MESSAGE now instead of regular PRINT     04/25/01 MCL
; Added /badval keyword                       08/27/2001 MDP
; Added /NaN keyword						  2003-04-12 MDP
; Added /Quick keyword						  2004-07-04 MDP
; Made /quick work with arbitrary size imgs	  2004-09-15 MDP
; Made /quick sense left and right sides independently
; 											  2005-03-16 MDP
;
;-
;###########################################################################
;
; LICENSE
;
; This software is OSI Certified Open Source Software.
; OSI Certified is a certification mark of the Open Source Initiative.
;
; Copyright © 1995-2003 by Michael Liu & Marshall Perrin
;   
; This software is provided "as-is", without any express or
; implied warranty. In no event will the authors be held liable
; for any damages arising from the use of this software.
;   
; Permission is granted to anyone to use this software for any
; purpose, including commercial applications, and to alter it and
; redistribute it freely, subject to the following restrictions:
;   
; 1. The origin of this software must not be misrepresented; you must
;    not claim you wrote the original software. If you use this software
;    in a product, an acknowledgment in the product documentation
;    would be appreciated, but is not required.
;   
; 2. Altered source versions must be plainly marked as such, and must
;    not be misrepresented as being the original software.
;   
; 3. This notice may not be removed or altered from any source distribution.
;   
; For more information on Open Source Software, visit the Open Source
; web site: http://www.opensource.org.
;   
;###########################################################################

  BADVAL = -1.e6
  
  if keyword_set(npix) eq 0 then npix = 8
  if n_params() lt 2 then begin
     print,'pro fixpix,imgs,badpix,outimgs,'
     print,'           [npix=',numformat(npix),'],[weight],'
     print,'           [noise],[sigma=],[dc=],[silent]'
     retall
  endif
  
  ;;if n_elements(outimgs) ne 0 then $
  ;;  if (total(outimgs-imgs) eq 0) then $
  ;;	message,'original and output images are the same!'
  outimgs = imgs  ;; don't change originals!
  
  if keyword_set(noise) and not(keyword_set(sigma)) then begin
     read,'enter sigma for random noise: ',sigma  
     read,'enter dc level for noise: ',dc
  endif
  if keyword_set(noise) and not(keyword_set(dc)) then dc=0
  
  sz = size(imgs)
  if sz(0) eq 2 then begin
     nimg = 1 
     img = imgs
  endif else if sz(0) eq 3 then begin
     nimg = sz(3) 
     img = imgs(*,*,0)
  endif else begin
     print,'** input image is not 2 or 3-d! **'
     return
  endelse
  
  if keyword_set(NaN) then begin
     wnotfinite = where(finite(imgs) eq 0,notfinitecount)
     if notfinitecount eq 0 then return ; bail if no bad pix
     badpix=fltarr(sz[1],sz[2],nimg)+1B
     badpix[wnotfinite]=0
  endif
  if keyword_set(badvalmask) then begin
     badpix=fltarr(sz[1],sz[2],nimg)+1B
     wbad = where(imgs eq BADVAL,badcount)
     if badcount gt 0 then badpix[wbad]=0
  endif 
  
  ;; check there are any bad pixels at all
  wbad = where(badpix eq 0.0,badcount)
  if badcount le 0 then begin
     ;; already set outimgs = imgs above, so just return
     return
  endif
  
  szb = size(badpix)
  if szb(0) eq 3 and szb(3) ne sz(3) then begin
     message, '3d badpix file different size than image file',/info 
     return
  endif
  
  ;; loop through images
  for j=0,nimg-1 do begin
     if szb(0) eq 3 then bp = badpix(*, *, j) $
     else bp = badpix
     
     wbp = where(bp eq 0.0, nbad)
     if not keyword_set(silent) then begin
        if (nimg gt 1) then print, format = '("image ",A,"/",A,": ")', numformat(j),numformat(nimg-1)
        message, numformat(nbad)+' bad pixels to fix', /info
     endif
     if sz(0) eq 3 then img = imgs(*,*,j)
     imoffset = j * sz(1) * sz(2)
     
     ;; quickly deal with large contiguous regions of bad pixels at the left
     ;; and/or right side of an image.
     if keyword_set(quick) and (bp[0,0] eq 0) then begin
        ;; quick fixing of the left edge
        badleft = search2d(bp,0,0,0,0)
        wgood = where(bp eq 1)
        med = median(img[wgood])
        outimgs[imoffset+badleft]=med		
        bp[badleft]=1
        wbp = where(bp eq 0.0, nbad)
     endif
     if keyword_set(quick) and (bp[sz[1]-1,0] eq 0) then begin
        ;; quick fixing of the right edges
        badright = search2d(bp,sz[1]-1,0,0,0)
        wgood = where(bp eq 1)
        med = median(img[wgood])
        outimgs[imoffset+badright]=med		
        bp[badright]=1
        wbp = where(bp eq 0.0, nbad)
     endif
     
     ;; loop through pixels
     if keyword_set(noise) then begin
        ;; replace with noise
        
        for i=0L,nbad-1 do begin
           newval = randomn(seed) * sigma + dc
           outimgs(wbp(i)+imoffset) = newval
        endfor                
     endif else begin
        ;; interpolate

        if ~keyword_set(search_box) then ddmax = 1d4 else ddmax = search_box
        
        for i=0L,nbad-1 do begin         
           ;; default search box is 2*dd+1 on a side
           dd    = 2
           found = -1

           ;; check if there is any good pixel within the search_box
           if keyword_set(search_box) then begin
              y = floor(wbp(i)/sz(1))
              x = wbp(i) - y*sz(1)
              x0 = x-dd > 0
              x1 = x+dd < (sz(1)-1)
              y0 = y-dd > 0
              y1 = y+dd < (sz(2)-1)
              bpcut = bp(x0:x1,y0:y1)              
              wgood = where(bpcut ne 0.0,ngood)
              if (ngood lt npix) then continue
           endif
           
           ;; determine search region
           repeat begin
              y = floor(wbp(i)/sz(1))
              x = wbp(i) - y*sz(1)
              x0 = x-dd > 0
              x1 = x+dd < (sz(1)-1)
              y0 = y-dd > 0
              y1 = y+dd < (sz(2)-1)
              bpcut = bp(x0:x1,y0:y1)
              cut = img(x0:x1,y0:y1)
              wgood = where(bpcut ne 0.0,ngood)
              if ngood lt npix then begin
                 dd = dd + 2
              endif else found = 1
              if (dd gt ddmax) then break
           endrep until (found ne -1)
           if (dd gt ddmax) then continue
           
           ;; calculate distances to adjacent good pixels 
           dist_circle,distarr,2*dd+1,dd,dd
           gdist = distarr(wgood)
           gpix = cut(wgood)
           
           ;; sort good pixels by distance
           ss = sort(gdist)
           gdist = gdist(ss)
           gpix = gpix(ss)
           
           ;; accounting for pixels with the same distance at the edge
           mm = where(gdist(npix-1:*) eq gdist(npix-1),nn)
           nn = nn - 1
           ;; if nn gt 0 then print,'  multiplicty = ',numformat(nn+1)
           
           ;; error checking
;;            if i eq 0 then begin 
;;               for k=0,npix+nn-1 do print,gdist(k),gpix(k)
;;               print,x,y
;;               print,cut
;;            endif
           
           ;; get the values of the relevant pixels
           gpix = gpix(0:npix+nn-1)
           ss2 = sort(gpix)
           gpix = gpix(sort(gpix))
           gdist = gdist(ss2)
           
           ;; calculate new pixel value, tossing the highest
           ;; and lowest pixels of the bunch, then weighting
           ;; by the inverse of the distances if desired
           if keyword_set(weight) then begin
              newval = total(gpix(1:npix+nn-2)/gdist(1:npix+nn-2))
              newval = newval / total(1./gdist(1:npix+nn-2))
           endif else begin
              newval = total(gpix(1:npix+nn-2)) / (npix+nn-2.0)
           endelse
           
           ;; more error checking
           if i eq 0 and keyword_set(silent) eq 0 then $
              print, '  oldval=', numformat(outimgs(wbp(i)+imoffset)), ' ' + $
                     'newval = ', numformat(newval) 
           
           outimgs(wbp(i)+imoffset) = newval        
        endfor      
     endelse
  endfor  
end
