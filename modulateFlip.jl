# [HFILT] = modulateFlipShift(LFILT)
#
# QMF/Wavelet highpass filter construction: modulate by (-1)^n,
# reverse order (and shift by one, which is handled by the convolution
# routines).  This is an extension of the original definition of QMF's
# (e.g., see Simoncelli90).

# Eero Simoncelli, 7/96.
#
# R. Crandall 8/13 ported to Julia

function modulateFlip(lfilt)

sz = length(lfilt)
sz2 = ceil(sz/2)

ind = sz:-1:1

hfilt = lfilt[ind].*(-1).^(ind-sz2)

end