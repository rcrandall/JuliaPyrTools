# RES = pyrBand(PYR, INDICES, BAND_NUM)
#
# Access a subband from a pyramid (gaussian, laplacian, QMF/wavelet, 
# or steerable).  Subbands are numbered consecutively, from finest
# (highest spatial frequency) to coarsest (lowest spatial frequency).

# Eero Simoncelli, 6/96.

# R. Crandall 8/13 ported to Julia from Matlab

function pyrBand(pyr, pind, band)

res = reshape( pyr[pyrBandIndices[pind,band]], pind[band,1], pind[band,2] )

end
