using Profile
using BenchmarkTools

include("finiteVolume.jl")
include("mesh.jl")
include("output.jl")
include("shockTube.jl")

################## Output ##################
println("Reading mesh")

### UnstructuredFVM from OpenFOAM Meshes ###

### Supersonic Wedge ###
# Freestream Conditions
P = 100000
T = 300
U = [ 694.26, 0, 0 ]

# Choose a mesh
# meshPath = "OFcoarseWedgeMesh"
# boundaryConditions = [ emptyBoundary, [], supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [], symmetryBoundary, [], zeroGradientBoundary, [], wallBoundary, [] ]
meshPath = "OFmemesWedgeMesh"
boundaryConditions = [ symmetryBoundary, [], emptyBoundary, [], supersonicInletBoundary, [P, T, U...], wallBoundary, [], zeroGradientBoundary, [], zeroGradientBoundary, [] ]

# Solve
mesh = OpenFOAMMesh(meshPath)
cellPrimitives = unstructured3DFVM(mesh, initializeUniformSolution3D(mesh, P, T, U...), boundaryConditions, ShuOsher, initDt=0.0000001, endTime=0.002, targetCFL=1.5, silent=false, restart=false)

outputVTK(meshPath, cellPrimitives)
P = cellPrimitives[:,1]
println("Plotting results")
plot2DResult(mesh, P)

#### Simple convection ####
# mesh, cellPrimitives = initializeShockTube3DFVM(10)
# P = 100000
# T = 300
# U = [ 694.26, 0, 0 ]
# boundaryConditions = [ supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [], emptyBoundary, [] ]
# P, U, T, rho = unstructured3DFVM(mesh, initializeUniformSolution3D(mesh, P, T, U...), boundaryConditions, initDt=0.0000001, endTime=0.00002, targetCFL=0.1, silent=false)
# plotShockTubeResults_PyPlot(P, U, T, rho)
