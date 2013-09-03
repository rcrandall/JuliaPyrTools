
require("buildSFpyr")
require("reconSFpyr")


#A = 255.0*ones(Float64,256,256)
#A[10:30,40:80] = 124.0

A = abs(120.0*randn(512,512))

@time (pyr,pind,steer,harm) = buildSFpyr(A,3,7)

pind = convert(Array{Int64,2},pind)
@time res = reconSFpyr(pyr,pind)

println(string("Error: ", norm(A-res)/norm(A)))