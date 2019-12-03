using Profile
using BenchmarkTools

include("finiteVolume.jl")
include("mesh.jl")
include("output.jl")

################## Output ##################
println("Reading mesh")

### UnstructuredFVM from OpenFOAM Meshes ###
mesh = OpenFOAMMesh("OFcoarseWedgeMesh")
P = 100000
T = 300
U = [ 694.26, 0, 0 ]
cellPrimitives = initializeUniformSolution3D(mesh, P, T, U...)
boundaryConditions = [ supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [], symmetryBoundary, [], zeroGradientBoundary, [], wallBoundary, [], emptyBoundary, [] ]
P, U, T, rho = unstructured3DFVM(mesh, cellPrimitives, boundaryConditions, ShuOsher, initDt=0.0000001, endTime=0.002, targetCFL=0.95, silent=false)

println("Plotting results")
plot2DResult(mesh, P)
