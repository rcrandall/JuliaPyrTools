# MTX = steer2HarmMtx(HARMONICS, ANGLES, REL_PHASES)
#
# Compute a steering matrix (maps a directional basis set onto the
# angular Fourier harmonics).  HARMONICS is a vector specifying the
# angular harmonics contained in the steerable basis/filters.  ANGLES 
# (optional) is a vector specifying the angular position of each filter.  
# REL_PHASES (optional, default = 'even') specifies whether the harmonics 
# are cosine or sine phase aligned about those positions.
# The result matrix is suitable for passing to the function STEER.

# Eero Simoncelli, 7/96.

# R. Crandall 8/13 ported to Julia from Matlab

function steer2HarmMtx(harmonics, angles, evenoroddStr)

##=================================================================
### Optional Parameters:

# Make HARMONICS a row vector
#harmonics = harmonics[:];

numh = 2*size(harmonics,2) - any(harmonics == 0)

#angles = angles[:]

##=================================================================


  if isequal(evenoroddStr,"even")
    evenorodd = 0;
  elseif strcmp(evenoroddStr,"odd")
    evenorodd = 1;
  else
    error("EVEN_OR_ODD should be the string  even or odd")
  end


## Compute inverse matrix, which maps Fourier components onto 
## steerable basis.
imtx = zeros(length(angles),numh)
col = 1
for h = harmonics
  args = h*angles
  if h == 0
    imtx[:,col] = ones(length(angles))
    col += 1
  elseif evenorodd
    imtx[:,col] = sin(args)
    imtx[:,col+1] = -cos(args)
    col += 2
  else
    imtx[:,col] = cos(args)
    imtx[:,col+1] = sin(args)
    col += 2;
  end
end
  
r = rank(imtx)
if  r != numh  &&  r != length(angles)
  println("WARNING: matrix is not full rank (imtx in steer2HarmMtx)")
end  

mtx = pinv(imtx)

end