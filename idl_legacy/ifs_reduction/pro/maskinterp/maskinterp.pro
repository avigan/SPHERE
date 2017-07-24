function maskinterp, rawdata, m, apstep, maxap, func, gpix=gpix, $
	gpoints=gpoints, cdis=cdis, debug=debug
;
; NAME:
;	maskinterp
;
; PURPOSE:
;	Given an image data array(2D) and a mask(of the same size)
;	which identifies a bad pixel on the image with a 0 and a 
;	good pixel with a 1, this function attempts to correct the
;       bad pixels with a user-specified interpolating function.  
; EXPLANATION:
;	Each bad pixel is fixed using only data from neighboring 
;	good pixels.  A circular aperture centered on the bad 
;	pixel determines the neighboring pixels to be used in the 
;	interpolation.  The aperture radius starts at zero and 
;	increases on successive iterations by an aperture step up 
;	to a maximum radius.  Since the number and the locations 
;	of the interpolating pixels vary, we set certain conditions 
;	for fixing the bad pixel.  The pixel is fixed only if (1) 
;	there is a sufficient number of good pixels, (2) the good 
;	pixels are not weighted to one side of the aperture, and 
;	(3) the conditions in the interpolating function are 
;	satisfied.  The method of determining the second condition 
;	is the calculation of the center of mass.  A good pixel is 
;	weighted by 1 and the bad pixel is weighted by 0.  The user 
;	specifies the maximum allowed distance between the center 
;	of mass and the center of the aperture, the minimum number 
;	of good pixels, and the minimum percentage of good pixels 
;	within the aperture.
;
; CALLING SEQUENCE:
;	Result = maskinterp(rawdata,m,apstep,maxap,"func")
;
; INPUTS: rawdata	A 2-dimensional data array. 
;         mask    	A 2-dimensional integral array of the same size as
;			rawdata. Allowed values are 1 and 0. 1 specifies
;			that the corresponding pixel in raw data is good.
;			0 specifies that it is bad.
;	  apstep	Aperture step. 
;	  maxap		Maximum radius of the aperture. 
;	  func		Interpolating function.
; KEYWORDS:
;	  gpix		Minimum percentage of good pixels in the aperture.
;	  gpoints    	Minimum number of good pixels in the aperture.
;	  cdis		Maximum distance from the center of aperture
;			to the "center of mass" of good pixels in the 
;			aperture.
;         debug		Set to 1 for debugging. If set, display a matrix
;			which identifies corrected pixels with a 1,
;			and the uncorrected pixels with a 0.
; OUTPUT: im		Interpolated data
; EXAMPLE:
;	raw = indgen(51,55)		;create image data array
;       mask = intarr(51,55)+1		;create a mask
;	mask(25,27) = 0
;	mask(23:27,4) = 0
;	z =  maskinterp(raw,mask,1,6,"gausfit",gpix=10,gpoints=5,cdis=2)
;	Interpolate with gausfit function, with at least 10% of good pixels
;       in an image aperture, and at least 5 good pixels.  Maximum distance
;	from the center of aperture to the center of mask is 2 pixels. 
; FUNCTIONS USED: keyword_defined,disc, dist_circle
; MODIFCATION HISTORY:
;	written, Siree Vatanavigkit, June 25, 1999
;       modified, Alex Ruane, April 7, 2000
;       modfied, W. Landsman, October 2001, use dist_circle for aperture mask
;determine values of keyword variables
on_error, 1
if not keyword_defined(gpix) then gpix = 20
if not keyword_defined(gpoints) then gpoints = 4
if not keyword_defined(cdis) then cdis = 3
if not keyword_defined(debug) then debug=0
if not isarray(m) then message, "m must be an array"  ;check mask data type
if not isarray(rawdata) then message, "rawdata must be an array"  ;check data's data type
if not isarray(where(m eq 0)) then return, rawdata     ;checks for lack of bad pixels

info=size(rawdata)
xdim=info(1)		; num of cols of rawdata
ydim=info(2)		; num of rows of rawdata
ap = apstep		; initalize aperture radius

;create a new mask which is the old mask padded with 0`s on the sides
minfo = size(m)           ; size of orginal mask
mxdim = minfo(1)+2*maxap  ; num of cols of new mask
mydim = minfo(2)+2*maxap  ; num of rows of new mask
mask = intarr(mxdim,mydim)
mask(maxap,maxap) = m

;create a new data matrix which is rawdata padded with 0's on the sides
im0 = dblarr(mxdim,mydim)
im0(maxap,maxap) = rawdata

;final mask describes if a bad pixel is fixed
fmask=intarr(xdim,ydim)

repeat begin		
   ;create aperture mask
   rad = ceil(ap)

   dist_circle, apmask, 1+2*rad      
   apmask = apmask LE ap

   ; number of good points in an aperture mask
   appts = n_elements(where(apmask eq 1))
   if(debug) then print, "ap rad ", rad
   for j = 0, ydim - 1 do begin		; traverse rows of rawdata
     for i = 0, xdim - 1 do begin	; traverse cols of rawdata
	   ;find a bad pixel	
           if (mask(i+maxap,j+maxap) eq 0) then begin
		umask = mask(i+maxap-rad:i+maxap+rad,j+maxap-rad:j+maxap+rad)*apmask		
        	;determine if conditions for fixing are met.
		if disc_interp(umask,appts,gpix,gpoints,cdis) then begin
                  s = call_function(func, $
                  	im0(i-rad+maxap:i+rad+maxap,j-rad+maxap:j+rad+maxap),$
			umask,debug=debug)
                  im0(i+maxap,j+maxap) = s.image
                  fmask(i,j) = s.fixed
		endif
           endif      
      endfor
   endfor
   ap = ap + apstep
   done = ( total(fmask(where (m eq 0))) eq n_elements(where (m eq 0)) )
endrep until ( (ap gt maxap) or done)
if(debug) then print, fmask+m
im = im0(maxap:mxdim-maxap-1,maxap:mydim-maxap-1)
return, im 
end
