function disc_interp,umask,appts,goodpix,goodpoints,cdis
;NAME: disc_interp
;PURPOSE:
;	This function determines whether a region of the
;	rawdata defined by umask is a good candidate
;	for interpolation.
;CALLING SEQUENCE:
;	f = disc(umask,appts,goodpix,goodpoints,cdis)
;INPUT:	umask		A 2-D matrix with 1 defining a good pixel 
;			and 0 defining a bad pixel.
;	appts		Number of good points in the aperture mask
;			Good points are those that are within
;			the radius of the aperture.
;	goodpix		Minimum percentage of good pixels
;	goodpoints  	Minimum number of good pixels
;	cdis		Minimum distance from the center of mass
;			to the center of aperture
;output: d		Set to 1 if aperture mask is good
;			for interpolation, 0 otherwise.
; Modified W. Landsman   October 2001   Remove loops
pos = where(umask eq 1)
points = n_elements(pos)
percent = 100*points/appts
info = size(umask)
xdim = info(1)	      ;always odd	
ydim = info(2)        ;always odd
xo = floor(xdim/2)
yo = floor(ydim/2)
;center of mass
xcm = total( pos mod xdim) / double(points) - xo
ycm = total( pos  /  xdim) / double(points) - yo

dist = sqrt( xcm^2 + ycm^2 )
d = ( points ge goodpoints ) and ( percent ge goodpix ) and (dist le cdis)
return, d
end
