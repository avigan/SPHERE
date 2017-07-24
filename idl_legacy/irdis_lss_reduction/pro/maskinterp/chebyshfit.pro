function chebyshfit, raw,mask,n=n,maxn=maxn,debug=debug
;+
;NAME: 	chebyshfit
;PURPOSE:  
;	An interpolating function implementing a 2D chebyshef
;	polynomial interpolation.
;	Given an image data array(2D) and a mask(of the same size)
;	which identifies a bad pixel on the image with a 0 and a 
;	good pixel with a 1, this function attempts to correct the
;       bad pixel at the center of the grid.  
;EXPLANATION:
;	Use curve fitting to find a function that best defines the 
;	values of the good pixels in the centered row and another 
;	function that best defines the values in the centered column. 
;	Plug in the centered coordinates to these two functions to
;	obtain two interpolated values. The result is the average of 
;	those numbers.
;	In general, all functions can be written in the following terms.
;		f(x) = sum(Ti*Ci)  - .5*Co,
;	where i = 0,1,...infinity, Ti`s are Chebyshef polynomials,
;	Ci's are constants.  
;	Chebyshef Polynomials:
;	To(x) = 1, T1(x) = x, Tn+1(x) = 2*x*Tn - Tn-1 
;	Let (xo,yo) be the coordinates of the centered pixel.
;	Let    x = a vector of the x coordinates of the good pixels,
;	       z = a vector of the image data corresponding to (x,yo), 
;	       coefs = a vector containing Ci's.
;       We obtain the following linear systems.
;	    	z = [To(x)-.5, T1(x), T2(x), ..., Tn(x)]*coefs
;	Call matrix [To(x)-.5, T1(x), T2(x), ..., Tn(x)], matrix A 
;	with m rows and n colums.  m must be greater than n+1 in order 
;	for the system to be solvable.  Since m corresponds to the number 
;	of good pixels in the centered row, m may be greater than n,
;	indicating that the matrix may be overconstrained. We use 
;	the method of least squares to find the best fit coefs.
;               coefs = inverse(transpose(A)*A)*z
; 	        result1 = sum(Ti(xo)*coefs(i)) - .5*coefs(0),
;	where i = 0,1,...n
;	Result2 in the y direction is obtained in the similar fashion with
;	the change from x to y and row to column.
;       	result = (result1+result2)/2
;CALLING SEQUENCE:
;	Result = chebyshfit(raw,mask)
;INPUTS:raw	A 2-dimensional data array.
;	mask	A 2-dimensional integral array of the same size as 
;		raw data.  Allowed values are 1 and 0. 1 specifies
;		that the corresponding pixel in rawdata is good.
;		0 specifies that it is bad.
;KEYWORDS:maxn	maximum degree of Chebyshef polynomial.  If not set, maxn
;		is default to 10. 
;	debug	Set to 1 to display the coeficients of the interpolated
;		functions. 
;OUTPUT:str     A structure with the following fields:
;		image 	Value of the corrected bad pixel.
;		fixed	Set to 1 if the bad pixel is corrected, 0 otherwise.
;EXAMPLE:
;	raw = indgen(5,7)	;create image data array
;	raw(2,3) = -99
;	mask = intarr(5,7)+1	;create mask
;	mask(2,3) = 0
;	ch = chebyshfit(raw,mask)
;FUNCTIONS USED: ludc,lusol,chvarchg
;MODIFICATION HISTORY:
;	WRITTEN: Siree Vatanavigkit. June 26,1999
;	MODIFIED: Alex Ruane. November 20, 2000
;		bug fixes, chvarchg application	
;       MODIFIED: W. Landsman October 2001, Use TOTAL() instead of ROWSUM	
;-
if not(keyword_defined(debug)) then debug = 0
if not(keyword_defined(maxn)) then maxn = 10

;determine position of badpixel to be interpolated
info = size(raw)
rangex = [0,(info[1]-1)]		
rangey = [0,(info[2]-1)]		
xo = info[1]/2
yo = info[2]/2

;determine points for interpolation
row = total(mask,1)    ;Row sum
col = total(mask,2)    ;Column sum
hn   = row(yo)
vn   = col(xo)

;restrictions for fit
fit = 1
;Need at least four points to create a good polynomial(at least degree3)
if (min([hn,vn]) lt 4) then fit =0

if fit then begin
	;row
	if (hn gt maxn) then hn = maxn
	pixel = where(mask(*,yo) eq 1)*1d
	z = raw(pixel,yo)
	for k = 0, n_elements(pixel)-1 do begin
	  pixel[k] = chvarchg(pixel[k], rangex)
	endfor
	A = dblarr(hn,n_elements(pixel))
	for j = 0,n_elements(pixel)-1 do begin
		for i = 0,hn-1 do begin
	    		case i of
			0: A[i,j] = 1d
			1: A[i,j] = pixel(j)
		     else: A[i,j] = 2*pixel(j)*A(i-1,j)-A(i-2,j)
		    endcase
		endfor
		A[0,j] = 0.5d
	endfor
	AA = transpose(A)##A  ;AA is sym pos definite
	b  = transpose(A)##z
	b  = transpose(b)
	;solve for least square ie. find x that minimizes (A^TAx=A^Tz)
	ludc,AA,in
	coefs = lusol(AA,in,b)
	xcoefs = coefs		;acr addition
	if (debug) then print,"Xcoefs",coefs
	;calculate val1
	xof = chvarchg(xo, rangex)
	T=dblarr(hn)
	T[0] =	1d
	T[1] =  xof
	for i = 2,hn-1 do begin
	     T[i] = 2*xof*T(i-1)-T(i-2)
	endfor
	T[0] = 0.5d
	T = transpose(T)
	val1 = coefs##T

	;col
	if (vn gt maxn) then vn = maxn		;hn =maxn
	pixel = where(mask(xo,*) eq 1)*1d
	z = raw(xo,pixel)
	for k = 0, n_elements(pixel)-1 do begin
	  pixel[k] = chvarchg(pixel[k], rangey)
	endfor
	A = dblarr(vn,n_elements(pixel))
	for j = 0,n_elements(pixel)-1 do begin
		for i = 0,vn-1 do begin
	    		case i of
			0: A[i,j] = 1d
			1: A[i,j] = pixel(j)
		     else: A[i,j] = 2*pixel(j)*A(i-1,j)-A(i-2,j)
		    endcase
		endfor
		A[0,j] = 0.5d
	endfor
	AA = transpose(A)##A  ;AA is sym pos definite
	b  = transpose(A)##z
	b  = transpose(b)
	;solve for least square ie. find x that minimizes (A^TAx=A^Tz)
	ludc,AA,in
	coefs = lusol(AA,in,b)
	ycoefs = coefs
	If (debug) then print,"Ycoefs",coefs
	;calculate val2
	yof = chvarchg(yo, rangey)
	T=dblarr(vn)
	T[0] =	1.0d
	T[1] =  yof
	for i = 2,vn-1 do begin
	     T[i] = 2*yof*T(i-1)-T(i-2)
	endfor
	T[0] = 0.5d
	T = transpose(T)
	val2 = coefs##T
	str = {image:(val1+val2)/2.,fixed:1}
endif else str = {image:raw(xo,yo),fixed:0}
return, str
end
