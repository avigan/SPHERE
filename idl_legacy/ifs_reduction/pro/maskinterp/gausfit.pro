function gausfit, raw,mask,debug=debug
;NAME: 	gausfit
;PURPOSE:  
;	An interpolating function implementing a gaussian surface fit.
;	Given an image data array(2D) and a mask(of the same size)
;	which identifies a bad pixel on the image with a 0 and a 
;	good pixel with a 1, this function attempts to correct the
;       bad pixel at the center of the grid.  
;EXPLANATIONS:
;	Fit a 2D gaussain function to the good data. By pluging 
;	in the coordinates of the centered pixel as variables to the 
;	fitted 2D gaussian equation, the interpolated value of the 
;	centered pixel is obtained.
;	   2D gaussian equation: 
;	z = (1/2*pi*thx)*exp(-.5*( ((x-mx)/thx)^2+((y-my)/thy)^2 ))- c
;	where mx,my are the means, thx,thy are the widths 
;       of the gaussian, and c is a positive number.
;       Rearranging the equation, we obtain the following.
;	 ln(z+c) = (-.5/thx^2)x^2 + (-mx/thx^2)x +
;		   (-.5/thy^2)y^2 + (-my/thy^2)y +
;		   (-.5((mx/thx)^2 +(my/thy)^2)+ln(1/2*pi*thx))
;	   Let x = a vector of the x coordinates of the good pixels,
;	       y = a vector of the y coordinates of the good pixels,
;	       z = a vector of the image data corresponding to (x,y), 
;	and coefs = a vector containing the coeficients [q r s t u] where
;	q = (-.5/thx^2) r = (-mx/thx^2) s = (-.5/thy^2) t = (-my/thy^2)
;       v =  (-.5((mx/thx)^2 +(my/thy)^2)+ln(1/2*pi*thx)).  		
;       We obtain the following linear system.
;	     ln(z+c) = [x^2 x y^2 y 1]*coefs
;	Call matrix [x^2 x y^2 y 1], matrix A with m rows and n colums.  
;	m must be greater than 5 in order for the system to be
;       solvable.  Since m corresponds to the number of good
;	pixels in the image, m may be greater then n, indicating that 
;       the matrix may have too many constraints. We use the method of 
;	least square to find the best fit coefs.
;               coefs = inverse(transpose(A)*A)*(ln(z+c))
;	   Let (xo,yo) be the coordinates of the centered pixel.
; 	        result = exp(qxo^2 + rxo + syo^2 + tyo + u) - c
;CALLING SEQUENCE:
;	Result = gausfit(raw,mask)
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
;	g = gausfit(raw,mask)
;FUNCTIONS USED:     coord.pro
;MODIFCATION HISTORY:
;	WRITTEN: Siree Vatanavigkit. June 26,1999
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
if (np ge 5) then fit = 1 else fit = 0

if (fit) then begin
	z = raw(points(0,0:np-1),points(1,0:np-1))  ;values in z-direction
	;shift z to make the lowest point in ln(z) > zero
	if ( min(z) le 0 ) then c = abs(min(z))+1 else c = 0
        z = alog(z + c)
	A = dblarr(5,np)+1
	A[0,0:np-1] = (points[0,0:np-1])^2
	A[1,0:np-1] =  points[0,0:np-1]
	A[2,0:np-1] = (points[1,0:np-1])^2
	A[3,0:np-1] =  points[1,0:np-1]
	AA = transpose(A)##A  ;AA is sym pos definite
	b  = transpose(A)##z
	b  = transpose(b)
	;solve for least square ie. find x that minimizes (A^TAx=A^Tz)
	ludc,AA,in
	coefs = lusol(AA,in,b)
	if(debug) then print, "gauscoefs",coefs
	val = coefs[0]*(xo^2)+coefs[1]*xo+coefs[2]*(yo^2)+coefs[3]*yo+coefs[4] 
        val = exp(val)-c
 	str = {image:val,fixed:1}
endif else str = {image:raw[xo,yo],fixed:0}
return, str
end
