# HEIGHT = maxPyrHt(IMSIZE, FILTSIZE)
#
# Compute maximum pyramid height for given image and filter sizes.
# Specifically: the number of corrDn operations that can be sequentially
# performed when subsampling by a factor of 2.

# Eero Simoncelli, 6/96.
#
# R. Crandall 8/13 ported to Julia from Matlab

function maxPyrHt(imszIn, filtszIn)

if any(imszIn .== 1) # 1D image
  imsz = prod(imszIn)
  filtsz = prod(filtszIn)
elseif any(filtszIn .== 1)              # 2D image, 1D filter
  filtsz = [filtszIn[1]; filtszIn[1]]
  imsz = imszIn
else
  filtsz = filtszIn
  imsz = imszIn
end

if any(imsz .< filtsz)
  height = 0
else
  height = 1 + maxPyrHt( floor(imsz/2), filtsz )
end

end