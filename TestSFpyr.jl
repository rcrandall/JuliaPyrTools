
require("buildSFpyr")
require("reconSFpyr")


A = 255.0*ones(Float64,256,256)
A[10:30,40:80] = 124.0


println("Not working with NOR > 0!")
(pyr,pind,steer,harm) = buildSFpyr(A,3,7)

pind = convert(Array{Int32,2},pind)
res = reconSFpyr(pyr,pind)

println(string("Error: ", norm(A-res)/norm(A)))