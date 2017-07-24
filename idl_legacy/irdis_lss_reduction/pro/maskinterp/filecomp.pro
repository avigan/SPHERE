pro filecomp,file1,file2,eps
;compare file 1 with file2
;file1: given .fits file
;file2: generated .fits file
;eps  : the tolerance of error
a = readfits(file1,ah)
x = sxpar(ah,'naxis1')
y = sxpar(ah,'naxis2')
nopix = x*y
error = total( abs( readfits(file1) - readfits(file2) ) ) /nopix
if ( error le eps ) then print, file2,":good" else $
print, file2, " has an average error per pixel equal to ",error
end
