function splinterp,raw,mask,debug=debug
;+
;NAME: 	splinterp
;PURPOSE:  
;	An interpolating function implementing bicubic spline
;	interpolation. Given an image data array(2D) and a mask
;	(of the same size) which identifies a bad pixel on the 
;	image with a 0 and a good pixel with a 1, this function 
;	attempts to correct the bad pixel at the center of the grid.
;EXPLANATION:
;	By applying a spline to the centered column of the data 
;	array and another spline to the centered row, we obtain 
;	two cubic splines.  Plugging in the x coordinate of the 
;	centered pixel to the row spline, and the y coordinate to 
;	the column spline, we obtain two interpolated values.  
;	The result is the average of these two values. 
;	See function csplinterp which implements a different
;	idl function(spline instead of spl_init and spl_interp) 
;	for finding splines.
;CALLING SEQUENCE:
;	Result = splinterp(raw,mask)
;INPUTS:raw	A 2-dimensional data array.
;	mask	A 2-dimensional integral array of the same size as 
;		raw data.  Allowed values are 1 and 0. 1 specifies
;		that the corresponding pixel in rawdata is good.
;		0 specifies that it is bad.
;KEYWORDS:debug No use for this keyword. Placed here to make function
;	        call compatible with maskinterp.pro
;OUTPUT:str     A structure with the following fields:
;		image 	Value of the corrected bad pixel.
;		fixed	Set to 1 if the bad pixel is corrected, 0 otherwise.
;EXAMPLE:
;	raw = indgen(5,7)	;create image data array
;	raw(2,3) = -99
;	mask = intarr(5,7)+1	;create mask
;	mask(2,3) = 0
;	s = splinterp(raw,mask)
;FUNCTIONS USED: spl_init,spl_interp,ludc,lusol
;MODIFICATION HISTORY:
;	WRITTEN: Siree Vatanavigkit. June 26,1999
;       MODIFIED: W. Landsman October 2001, Use TOTAL instead of ROWSUM()
;-
if not(keyword_defined(debug)) then debug = 0
;determine position of the centered bad pixel
info = size(raw)
xdim = info[1]
ydim = info[2]
xo = xdim/2
yo = ydim/2

;restriction for cubic spline interpolation
; need at least two good points in row xo and col yo
row = total(mask,1)       ;Row sum
col = total(mask,2)       ;Column sum
if ( row(xo) ge 2 and col(yo) ge 2) then fit = 1 $
else fit = 0

if fit then begin
	;interpolate from row, then column
	x =  where(mask(*,yo) eq 1)
	z =  raw(x,yo)
	d2 = spl_init(x,z,/double)
	val1 = spl_interp(x,z,d2,xo)

	y = where(mask(xo,*) eq 1)
	z = raw(xo,y)
	d2 = spl_init(y,z,/double)
	val2 = spl_interp(y,z,d2,yo)

	str = {image:(val1+val2)/2.,fixed:1}
endif else str = {image:raw[xo,yo],fixed:0}	
return, str
end
