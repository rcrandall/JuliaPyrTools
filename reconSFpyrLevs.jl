# RESDFT = reconSFpyrLevs(PYR,INDICES,LOGRAD,XRCOS,YRCOS,ANGLE,NBANDS,LEVS,BANDS)
#
# Recursive function for reconstructing levels of a steerable pyramid
# representation.  This is called by reconSFpyr, and is not usually
# called directly.

# Eero Simoncelli, 5/97.

require("pointOp")
require("pyrBand") 
require("pyrBandIndices")

function reconSFpyrLevs(pyr,pind,log_radIn,Xrcos,Yrcos,angle,nbands,levs,bands)

lo_ind = nbands+1

dims = pind[1,:]
ctr = ceil((dims+0.5)/2)

#log_rad = log_radIn + 1
log_rad = log_radIn
Xrcos = Xrcos - log2(2)

if any(levs .> 1)

  lodims = ceil((dims-0.5)/2)
  loctr = ceil((lodims+0.5)/2)
  lostart = ctr-loctr+1
  loend = lostart+lodims-1
  nlog_rad = log_rad[lostart[1]:loend[1],lostart[2]:loend[2]]
  nangle = angle[lostart[1]:loend[1],lostart[2]:loend[2]]

  if  size(pind,1) > lo_ind
    nresdft = reconSFpyrLevs( pyr[1+sum(prod(pind[1:lo_ind-1,:],2)):size(pyr,1)],
    	                        pind[lo_ind:size(pind,1),:], 
    	                        nlog_rad, Xrcos, Yrcos, nangle, nbands,levs-1, bands )
  else
    nresdft = fftshift(fft(pyrBand(pyr,pind,lo_ind)))
  end

  YIrcos = sqrt(abs(1.0 - Yrcos.^2))
  lomask = pointOp(nlog_rad, YIrcos, Xrcos[1], Xrcos[2]-Xrcos[1])
  dtmp = (convert(Int32,round(dims[1])),convert(Int32,round(dims[2])))
  resdft = zeros(Complex{Float64},dtmp)
  resdft[lostart[1]:loend[1],lostart[2]:loend[2]] = nresdft.*complex(lomask)

else

  resdft = zeros(dims[1],dims[2])

end
	
if any(levs .== 1)

  lutsize = 1024
  Xcosn = pi*[-(2*lutsize+1):(lutsize+1)]/lutsize  # [-2*pi:pi]
  order = nbands-1
  ## divide by sqrt(sum_(n=0)^(N-1)  cos(pi*n/N)^(2(N-1)) )
  aconst = (exp2(2*order))*(factorial(BigInt(order))^2)/(nbands*factorial(BigInt(2*order)))
  Ycosn = sqrt(aconst) * (cos(Xcosn)).^order
  himask = pointOp(log_rad, Yrcos, Xrcos[1], Xrcos[2]-Xrcos[1])

  ind = 1;
  for b = 1:nbands
    if any(bands .== b)
      anglemask = pointOp(angle,Ycosn,Xcosn[1]+pi*(b-1)/nbands,Xcosn[2]-Xcosn[1])
      band = reshape(pyr[ind:ind+prod(dims)-1], convert(Int32,dims[1]),convert(Int32,dims[2]))
      banddft = fftshift(fft(band))
      resdft = resdft + (1im)^(nbands-1) * banddft.*anglemask.*himask
    end
    ind = ind + prod(dims)
  end
end
return resdft

end