# [PYR, INDICES, STEERMTX, HARMONICS] = buildSFpyr(IM, HEIGHT, ORDER, TWIDTH)
#
# Construct a steerable pyramid on matrix IM, in the Fourier domain.
# This is similar to buildSpyr, except that:
#
#    + Reconstruction is exact (within floating point errors)
#    + It can produce any number of orientation bands.
#    - Typically slower, especially for non-power-of-two sizes.
#    - Boundary-handling is circular.
#
# HEIGHT (optional) specifies the number of pyramid levels to build. Default
# is maxPyrHt(size(IM),size(FILT));
#
# The squared radial functions tile the Fourier plane, with a raised-cosine
# falloff.  Angular functions are cos(theta-k\pi/(K+1))^K, where K is
# the ORDER (one less than the number of orientation bands, default= 3).
#
# TWIDTH is the width of the transition region of the radial lowpass
# function, in octaves (default = 1, which gives a raised cosine for
# the bandpass filters).
#
# PYR is a vector containing the N pyramid subbands, ordered from fine
# to coarse.  INDICES is an Nx2 matrix containing the sizes of
# each subband.  This is compatible with the MatLab Wavelet toolbox.
# See the function STEER for a description of STEERMTX and HARMONICS.

# Eero Simoncelli, 5/97.
# See http://www.cis.upenn.edu/~eero/steerpyr.html for more
# information about the Steerable Pyramid image decomposition.

# R. Crandall 8/13 ported to Julia from Matlab

require("rcosFn.jl")
require("pointOp.jl")
require("steer2HarmMtx.jl")
require("buildSFpyrLevs.jl")

function buildSFpyr(im, ht, orderIn)

#-----------------------------------------------------------------
## DEFAULTS:

max_ht = floor(log2(min(size(im)))+2)

#if (exist('ht') ~= 1)
#  ht = max_ht;
#else
  if ht > max_ht
    error(string("Cannot build pyramid, too many levels; max is ", max_ht))
  end
#end

#if (exist('order') ~= 1)
#  order = 3;
if (orderIn > 15)  || (orderIn < 0)
  println("Warning: ORDER must be an integer in the range [0,15]. Truncating")
  order = clamp(orderIn,0,15)
else
  order = round(orderIn)
end

nbands = order+1

#if (exist('twidth') ~= 1)
  twidth = 1;
#elseif (twidth <= 0)
#  fprintf(1,'Warning: TWIDTH must be positive.  Setting to 1.\n');
#  twidth = 1;
#end

#-----------------------------------------------------------------
## Steering stuff:

if iseven(nbands)
  harmonics = [0:(nbands/2)-1]*2 + 1
else
  harmonics = [0:(nbands-1)/2]*2
end

steermtx = steer2HarmMtx(harmonics, pi*[0:(nbands-1)]/nbands, "even")

#-----------------------------------------------------------------

dims = size(im)
ctr = ceil(([dims...]+0.5)/2)

xtmp = ([1:dims[2]]-ctr[2])./(dims[2]/2)
ytmp = ([1:dims[1]]-ctr[1])./(dims[1]/2)

angle = [atan2(x,y) for x = xtmp, y = ytmp]
log_rad = [sqrt(x.^2 + y.^2) for x = xtmp, y = ytmp]

#angle = atan2(yramp,xramp)
#log_rad = sqrt(xramp.^2 + yramp.^2)


log_rad[ctr[1],ctr[2]] =  log_rad[ctr[2],ctr[2]-1]
log_rad  = log2(log_rad)

## Radial transition function (a raised cosine in log-frequency):
(Xrcos,Yrcos) = rcosFn(twidth,-twidth/2,[0 1])
Yrcos = sqrt(Yrcos)

YIrcos = sqrt(1.0 - Yrcos.^2)

lo0mask = pointOp(log_rad, YIrcos, Xrcos[1], Xrcos[2]-Xrcos[1])
imdft = fftshift(fft(im))
lo0dft =  imdft .* lo0mask

(pyr,pind) = buildSFpyrLevs(lo0dft, log_rad, Xrcos, Yrcos, angle, ht, nbands)

hi0mask = pointOp(log_rad, Yrcos, Xrcos[1], Xrcos[2]-Xrcos[1])
hi0dft =  imdft .* hi0mask
hi0 = ifft(ifftshift(hi0dft))

pyr = [real(hi0[:]); pyr]
pind = [[size(hi0)...]'; pind]

return pyr,pind,steermtx,harmonics

end