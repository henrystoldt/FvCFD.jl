using Profile
using BenchmarkTools

include("finiteDifference.jl")
include("finiteVolume.jl")
include("mesh.jl")
include("output.jl")

################## Output ##################
println("Reading mesh")

### UnstructuredFVM from OpenFOAM Meshes ###
mesh = OpenFOAMMesh("OFwedgeMesh")
P = 100000
T = 300
U = [ 694.26, 0, 0 ]
cellPrimitives = initializeUniformSolution3D(mesh, P, T, U...)
P, U, T, rho = unstructured3DFVM(mesh, cellPrimitives, ShuOsher, initDt=0.0000001, endTime=0.0000045, targetCFL=0.95, silent=false)

println("Plotting results")
plot2DResult(mesh, P)
