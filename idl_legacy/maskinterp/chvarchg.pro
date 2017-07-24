function chvarchg, x, range
;NAME: chvarchg
;PURPOSE:
;	A change of variable function restraining a value to a range of -1, 1. 
;	Given a value in an original range, this function will identify the 
;	value's equivilent location in a range of (-1, 1).
;EXPLANATION:
;	To restrain a value to a range of (-1, 1) its location relative to the
;	center of the original range (a, b) must be determined using the 
;	function:
;		xrc = x - 0.5(b+a)
;	To then restrict the range to (-1, 1) from (a, b) you must divide by
;	one-half the distance covered by the range using the function:
;		xr0 = xrc / (0.5(b-a))
;	The resulting position is the original value's equivalent location on 
;	the (-1, 1) range.
;CALLING SEQUENCE:
;	result = chvarchg (x, range)
;INPUTS:
;	x	value to be redefined on the -1, 1 range.
;	range	range vector [beginning of range, end of range] 
;OUTPUT: 
;	xf	value of equivalent location on the (-1,1) range
;EXAMPLE:
;	range = [80, 100]
;	x = 94
;	varch = chvarchg(x, range)
;MODIFICAION HISTORY:
;	WRITTEN:	Alex Ruane. November, 16, 2000
x = x*1.0
xrc = x-0.5*(range[1]+range[0])
xr0 = xrc/(0.5*(range[1]-range[0]))
return, xr0

end