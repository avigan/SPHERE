function plsfit, raw,mask,debug=debug
;NAME: 	plsfit
;PURPOSE:  
;	An interpolating function implementing a plane surface fit.
;	Given an image data array(2D) and a mask(of the same size)
;	which identifies a bad pixel on the image with a 0 and a 
;	good pixel with a 1, this function attempts to correct the
;       bad pixel at the center of the grid.  
;EXPLANATION:
;	Fit a plane to the good data.  By plugging in the coordinates 
;	of the centered pixel as variables to the fitted planar 
;	equation, the interpolated value of the centered pixel is 
;	obtained.
;	   Planar equation: z = ax + by + c
;	   Let x = a vector of the x coordinates of the good pixels,
;	       y = a vector of the y coordinates of the good pixels,
;	       z = a vector of the image data corresponding to (x,y), 
;	 and coefs = a vector containing the coefficients: [c a b].
;       We obtain the following linear system.
;	    	z = [1 x y]*coefs
;	Call matrix [1 x y], matrix A with m rows and n colums.  
;	m must be at least 3 in order for the system to be
;       solvable.  Since m corresponds to the number of good
;	pixels in the image, m may be greater than n, indicating that 
;       the matrix may have too many constraints. We use the method of 
;	least squares to find the best fit coefficients.
;               coefs = inverse(transpose(A)*A)*z
;	   Let (xo,yo) be the coordinates of the centered pixel.
; 	        result = c + a*xo + b*yo
;CALLING SEQUENCE:
;	Result = plsfit(raw,mask)
;INPUTS:raw	A 2-dimensional data array.
;	mask	A 2-dimensional integral array of the same size as 
;		raw data.  Allowed values are 1 and 0. 1 specifies
;		that the corresponding pixel in rawdata is good.
;		0 specifies that it is bad.
;KEYWORDS:
;	debug	Set to 1 to show the coeficients of the planar equation.
;OUTPUT:str     A structure with the following fields:
;		image 	Value of the corrected bad pixel.
;		fixed	Set to 1 if the bad pixel is corrected, 0 otherwise.
;EXAMPLE:
;	raw = indgen(5,7)	;create image data array
;	raw(2,3) = -99
;	mask = intarr(5,7)+1	;create mask
;	mask(2,3) = 0
;	p = plsfit(raw,mask)
;FUNCTIONS USED:  keyword_defined,coord,ludc,lusol
;MODFICATION HISTORY:
;	WRITTEN: Siree Vatanavigkit, June 26,1999

if not(keyword_defined(debug)) then debug = 0
;determine position of the centered bad pixel
info = size(raw)
xo = info[1]/2
yo = info[2]/2
;determine points for interpolation
points = coord(mask)
sp = size(points)
np = sp[2]	; num of fitting points

;restrictions of plane surface fit
;At least three good pixels that are not lined up in one direction
fit = 0
if(np ge 3) then begin
	;check if x's are different
        c = 0
        difx = 0
	while( (c lt np-1) and ( difx eq 0 )) do begin
	    if (points[0,c] ne points[0,c+1] ) then $
                difx = 1
	    c = c + 1
 	endwhile
        if (difx eq 1) then begin
		m = ( points[1,c]-points[1,c-1] ) / ( points[0,c] - points[0,c-1] )
		c = 2
		while ( (c lt np) and (fit eq 0) ) do begin
		    left = points[1,c] - points[1,0]
		    right = m*(points[0,c] - points[0,0])
       		     if ( left ne right ) then fit = 1
		     c = c+1
		endwhile
	endif 
endif 
if (fit) then begin
	z = raw(points[0,0:np-1],points[1,0:np-1])  ;values in z-direction
	A = dblarr(3,np)+1
	A[1,0:np-1] = points[0,0:np-1]
	A[2,0:np-1] = points[1,0:np-1]
	AA = transpose(A)##A  ;AA is sym pos definite
	b  = transpose(A)##z
	b  = transpose(b)
	;solve for least square ie. find x that minimizes (A^TAx=A^Tz)
	ludc,AA,in
	coefs = lusol(AA,in,b)
	if (debug) then print,"plsfit coefs:",coefs
	val = coefs[0]+coefs[1]*xo+coefs[2]*yo 
	str = {image:val,fixed:1}
endif else str = {image:raw[xo,yo],fixed:0}
return, str
end
