# [PYR, INDICES] = buildSFpyrLevs(LODFT, LOGRAD, XRCOS, YRCOS, ANGLE, HEIGHT, NBANDS)
#
# Recursive function for constructing levels of a steerable pyramid.  This
# is called by buildSFpyr, and is not usually called directly.

# Eero Simoncelli, 5/97.

# R. Crandall 8/13 ported to Julia from Matlab

require("pointOp.jl")

function buildSFpyrLevs(lodft,log_radIn,Xrcos,Yrcos,angle,ht,nbandsIn)

nbands = convert(Int64,nbandsIn)
if ht <= 0

  lo0 = ifft(ifftshift(lodft))
  pyr = real(lo0[:])
  pind = [size(lo0)...]'

else

  bands = zeros(length(lodft), nbands)
  bind = zeros(nbands,2)

  log_rad = log_radIn + 1

  lutsize = 1024
  Xcosn = pi*[-(2*lutsize+1):(lutsize+1)]/lutsize  # [-2*pi:pi]
  order = convert(Int64,nbands-1)
  ## divide by sqrt(sum_(n=0)^(N-1)  cos(pi*n/N)^(2(N-1)) )
  ## Thanks to Patrick Teo for writing this out :)
  aconst = (exp2(2*order))*(factorial(order)^2)/(nbands*factorial(2*order))
  Ycosn = sqrt(aconst) * (cos(Xcosn)).^order
  himask = pointOp(log_rad, Yrcos, Xrcos[1], Xrcos[2]-Xrcos[1])

  for b = 1:nbands
    anglemask = pointOp(angle, Ycosn, Xcosn[1]+pi*(b-1)/nbands, Xcosn[2]-Xcosn[1])
    banddft = ((-1im)^(nbands-1)) .* lodft .* anglemask .* himask
    band = ifft(ifftshift(banddft))

    bands[:,b] = real(band[:])
    bind[b,1:2]  = [size(band)...]
  end

  dims = [size(lodft)...]
  ctr = ceil((dims+0.5)/2)
  lodims = ceil((dims-0.5)/2)
  loctr = ceil((lodims+0.5)/2)
  lostart = ctr-loctr+1
  loend = lostart+lodims-1

  log_rad = log_rad[lostart[1]:loend[1],lostart[2]:loend[2]]
  angleNew = angle[lostart[1]:loend[1],lostart[2]:loend[2]]
  lodftNew = lodft[lostart[1]:loend[1],lostart[2]:loend[2]]
  YIrcos = abs(sqrt(1.0 - Yrcos.^2));
  lomask = pointOp(log_rad, YIrcos, Xrcos[1], Xrcos[2]-Xrcos[1])

  lodftNew = lomask .* lodftNew;

  (npyr,nind) = buildSFpyrLevs(lodftNew, log_rad, Xrcos, Yrcos, angleNew, ht-1, nbands)

  pyr = [bands[:]; npyr]

  pind = [bind; nind]

end

  return pyr, pind

end