# RES = reconSFpyr(PYR, INDICES, LEVS, BANDS, TWIDTH)
#
# Reconstruct image from its steerable pyramid representation, in the Fourier
# domain, as created by buildSFpyr.
#
# PYR is a vector containing the N pyramid subbands, ordered from fine
# to coarse.  INDICES is an Nx2 matrix containing the sizes of
# each subband.  This is compatible with the MatLab Wavelet toolbox.
#
# LEVS (optional) should be a list of levels to include, or the string
# 'all' (default).  0 corresonds to the residual highpass subband.  
# 1 corresponds to the finest oriented scale.  The lowpass band
# corresponds to number spyrHt(INDICES)+1.
#
# BANDS (optional) should be a list of bands to include, or the string
# 'all' (default).  1 = vertical, rest proceeding anti-clockwise.
#
# TWIDTH is the width of the transition region of the radial lowpass
# function, in octaves (default = 1, which gives a raised cosine for
# the bandpass filters).

# Eero Simoncelli, 5/97.

require("spyrNumBands")
require("spyrHt")
require("pointOp")
require("reconSFpyrLevs")

function reconSFpyr(pyr, pind)

##------------------------------------------------------------
## DEFAULTS:

#if (exist('levs') ~= 1)
  levs = "all"
#end

#if (exist('bands') ~= 1)
  bands = "all"
#end

#if (exist('twidth') ~= 1)
  twidth = 1;
#elseif (twidth <= 0)
#  fprintf(1,'Warning: TWIDTH must be positive.  Setting to 1.\n');
#  twidth = 1;
#end

##------------------------------------------------------------

nbands = spyrNumBands(pind)

maxLev =  1+spyrHt(pind)
if isequal(levs,"all")
  levs = 0:maxLev
#else
#  if (any(levs > maxLev) | any(levs < 0))
#    error(sprintf('Level numbers must be in the range [0, #d].', maxLev));
#  end
#  levs = levs(:);
end

if isequal(bands,"all")
  bands = 1:nbands
#else
#  if (any(bands < 1) | any(bands > nbands))
#    error(sprintf('Band numbers must be in the range [1,3].', nbands));
#  end
#  bands = bands(:);
end

#----------------------------------------------------------------------

dims = pind[1,:]
ctr = ceil((dims+0.5)/2)

xtmp = [(1:dims[2])-ctr[2]]./(dims[2]/2)
ytmp = [(1:dims[1])-ctr[1]]./(dims[1]/2) 

angle = [atan2(y,x) for y = ytmp, x = xtmp]
log_rad = [sqrt(x.^2 + y.^2) for x = xtmp, y = ytmp]

log_rad[ctr[1],ctr[2]] =  log_rad[ctr[1],ctr[2]-1]
log_rad = log2(log_rad)

## Radial transition function (a raised cosine in log-frequency):
(Xrcos,Yrcos) = rcosFn(twidth,(-twidth/2),[0 1])
Yrcos = sqrt(Yrcos)
YIrcos = sqrt(abs(1.0 - Yrcos.^2))

if (size(pind,1) == 2)
  if (any(levs.==1))
    resdft = fftshift(fft(pyrBand(pyr,pind,2)))
  else
    resdft = zeros(pind[2,:])
  end
else
  resdft = reconSFpyrLevs(pyr[1+prod(pind[1,:]):size(pyr,1)], 
                          pind[2:size(pind,1),:],
      log_rad, Xrcos, Yrcos, angle, nbands, levs, bands)
end

lo0mask = pointOp(log_rad, YIrcos, Xrcos[1], Xrcos[2]-Xrcos[1])
resdft = resdft .* lo0mask




## residual highpass subband
if any(levs .== 0)
  hi0mask = pointOp(log_rad, Yrcos, Xrcos[1], Xrcos[2]-Xrcos[1])

  hidft = fftshift(fft(reshape(pyr[1:prod(pind[1,:])],convert(Int32,pind[1,1]),convert(Int32,pind[1,2]))));	
  #hidft = fftshift(fft(subMtx(pyr, pind[1,:])))

  resdft = resdft + hidft .* complex(hi0mask)
end
 
res = real(ifft(ifftshift(resdft)))

end