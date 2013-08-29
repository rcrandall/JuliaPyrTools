# RES = pointOp(IM, LUT, ORIGIN, INCREMENT, WARNINGS)
#
# Apply a point operation, specified by lookup table LUT, to image IM.
# LUT must be a row or column vector, and is assumed to contain
# (equi-spaced) samples of the function.  ORIGIN specifies the
# abscissa associated with the first sample, and INCREMENT specifies the
# spacing between samples.  Between-sample values are estimated via
# linear interpolation.  If WARNINGS is non-zero, the function prints
# a warning whenever the lookup table is extrapolated.
#
# This function is much faster than MatLab's interp1, and allows
# extrapolation beyond the lookup table domain.  The drawbacks are
# that the lookup table must be equi-spaced, and the interpolation is
# linear.

# Eero Simoncelli, 8/96.

# R. Crandall 8/13 ported to Julia from Matlab.  Original code below

import Grid

function pointOp(im,lut,origin,increment)

yi = Grid.InterpGrid(lut, Grid.BCnearest, Grid.InterpLinear)

res = map(x -> yi[x],(im-origin)/increment + 1.0)
#res = zeros(size(im))

#for i = 1:length(res)

#	res[i] = yi[(im[i]-origin)/increment+1]
#end

end

# Original Matlab code: 

#function res = pointOp(im, lut, origin, increment, warnings)

#fprintf(1,'WARNING: You should compile the MEX code for "pointOp", found in the MEX subdirectory.  It is MUCH faster.\n');

#X = origin + increment*[0:(length(lut)-1)];
#Y = lut;

#res = reshape(interp1(X, Y, im(:), 'linear'),size(im));

#end