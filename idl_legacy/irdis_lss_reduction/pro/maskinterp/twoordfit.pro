function twoordfit, raw, mask, debug=debug
;NAME: 	twoordfit
;PURPOSE:  
;	An interpolating function implementing second-order surface fit.
;	Given an image data array(2D) and a mask(of the same size)
;	which identifies a bad pixel on the image with a 0 and a 
;	good pixel with a 1, this function attempts to correct the
;       bad pixel at the center of the grid.  
;EXPLANATIONS:
;	Fit a second-order suface to the good data. By plugging 
;	in the coordinates of the centered pixel as variables to the 
;	fitted 2D gaussian equation, the interpolated value of the 
;	centered pixel is obtained.
;	   2D gaussian equation: 
;	z = a*xo^2 + b*xo*yo + c*yo^2 + d*xo + eyo + f
;	where xo, yo are the coordinates of the centered pixel.
;	   Let x = a vector of the x coordinates of the good pixels,
;	       y = a vector of the y coordinates of the good pixels,
;	       z = a vector of the image data corresponding to (x,y), 
;	and coefs = a vector containing the coeficients [a b c d e f] that 
;	satisfy the system  
;	       z = [x^2 xy y^2 x y 1]*coefs
;	Call matrix [x^2 xy y^2 x y 1], matrix A with m rows and n colums.  
;	m must be greater than 6 in order for the system to be
;       solvable.  Since m corresponds to the number of good
;	pixels in the image, m may be greater then n, indicating that 
;       the matrix may have too many constraints. We use the method of 
;	least square to find the best fit coefs.
;               coefs = inverse(transpose(A)*A)*(ln(z+c))
;	   Let (xo,yo) be the coordinates of the centered pixel.
; 	        result = axo^2 + bxo + cyo^2 + dyo + exoyo + f
;CALLING SEQUENCE:
;	Result = twoordfit(raw,mask)
;INPUTS:raw	A 2-dimensional data array.
;	mask	A 2-dimensional integral array of the same size as 
;		raw data.  Allowed values are 1 and 0. 1 specifies
;		that the corresponding pixel in rawdata is good.
;		0 specifies that it is bad.
;KEYWORDS:
;	debug	Set to 1 to show the coeficients of the 2D gaussian equation.
;OUTPUT:str     A structure with the following fields:
;		image 	Value of the corrected bad pixel.
;		fixed	Set to 1 if the bad pixel is corrected, 0 otherwise.
;EXAMPLE:
;	raw = indgen(5,7)	;create image data array
;	raw(2,3) = -99
;	mask = intarr(5,7)+1	;create mask
;	mask(2,3) = 0
;	g = twoordfit(raw,mask)
;FUNCTIONS USED:     coord.pro
;MODIFCATION HISTORY:
;	WRITTEN: Alex Ruane, May, 18, 2000
if not(keyword_defined(debug)) then debug = 0
;determine position of the centered bad pixel
info = size(raw)
xo = info[1]/2
yo = info[2]/2
;determine points for interpolation
points = coord(mask)
sp = size(points)
np = sp[2]	; num of fitting points

;restrictions of gaussian surface fit
if (np ge 8) then fit = 1 else fit = 0

if (fit) then begin
	z = raw(points(0,0:np-1),points(1,0:np-1))  ;values in z-direction
	A = dblarr(6,np)+1
	A[0,0:np-1] = (points[0,0:np-1])^2
	A[1,0:np-1] = (points[0,0:np-1])*(points[1,0:np-1])
	A[2,0:np-1] = (points[1,0:np-1])^2
	A[3,0:np-1] =  points[0,0:np-1]
	A[4,0:np-1] =  points[1,0:np-1]
	AA = transpose(A)##A  ;AA is sym pos definite
	b  = transpose(A)##z
	b  = transpose(b)
	;solve for least square ie. find x that minimizes (A^TAx=A^Tz)
	ludc,AA,in
	coefs = lusol(AA,in,b)
	if(debug) then print, "coefs",coefs
	val = (coefs[0]*(xo^2)+coefs[1]*xo*yo+coefs[2]*(yo^2)+coefs[3]*xo+coefs[4]*yo+coefs[5])
 	str = {image:val,fixed:1}
endif else str = {image:raw[xo,yo],fixed:0}
return, str
end
