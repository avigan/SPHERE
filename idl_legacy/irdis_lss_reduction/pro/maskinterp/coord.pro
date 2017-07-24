function coord, mask
;Given a 2D matrix, this function returns the coordinates of
;the points whose value is equal to 1.
; October 2001   W. Landsman     Remove loops
info = size(mask)
xdim = info(1)
ydim = info(2)
pos  = where(mask eq 1)
x    = pos mod xdim
y    = pos / xdim 

points = transpose([[x],[y]])
return, points
end
