using Profile
using BenchmarkTools
using ProfileView

include("finiteVolume.jl")
include("mesh.jl")
include("output.jl")
include("shockTube.jl")

################## Output ##################
println("Reading mesh")

### UnstructuredFVM from OpenFOAM Meshes ###

### Supersonic Wedge ###
# Freestream Conditions
P = 1
T = 1
U = [ 3, 0, 0 ]

# Choose a mesh
# meshPath = "OFunstructuredWedgeMesh"
# boundaryConditions = [ emptyBoundary, [], supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [], symmetryBoundary, [], zeroGradientBoundary, [], wallBoundary, [] ]
# meshPath = "OFwedgeMesh"
# boundaryConditions = [ supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [], symmetryBoundary, [], zeroGradientBoundary, [], wallBoundary, [], emptyBoundary, [] ]
# meshPath = "OFcoarseWedgeMesh"
# boundaryConditions = [ emptyBoundary, [], supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [], symmetryBoundary, [], zeroGradientBoundary, [], wallBoundary, [] ]
# meshPath = "OFmemesWedgeMesh"
# boundaryConditions = [ symmetryBoundary, [], emptyBoundary, [], supersonicInletBoundary, [P, T, U...], wallBoundary, [], zeroGradientBoundary, [], zeroGradientBoundary, [] ]

meshPath = "OFforwardStepMesh"
boundaryConditions = [ supersonicInletBoundary, [P, T, U..., 2.5], zeroGradientBoundary, [], symmetryBoundary, [], symmetryBoundary, [], wallBoundary, [], emptyBoundary, [] ]

# Solve
mesh = OpenFOAMMesh(meshPath)
cellPrimitives = initializeUniformSolution3D(mesh, P, T, U...)
unstructured3DFVM(mesh, meshPath, cellPrimitives, boundaryConditions, ShuOsher, initDt=0.0000001, endTime=4, outputInterval=0.5, targetCFL=0.95, silent=false, restart=false, gamma=1.4, R=0.714286, Cp=2.5)
# @time unstructured3DFVM(mesh, meshPath, cellPrimitives, boundaryConditions, initDt=0.0000001, endTime=0.0002, outputInterval=0.0002, targetCFL=0.5, silent=false, restart=false)


#### Simple convection ####
# mesh, cellPrimitives = initializeShockTube3DFVM(10)
# P = 100000
# T = 300
# U = [ 694.26, 0, 0 ]
# boundaryConditions = [ supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [], emptyBoundary, [] ]
# P, U, T, rho = unstructured3DFVM(mesh, initializeUniformSolution3D(mesh, P, T, U...), boundaryConditions, initDt=0.0000001, endTime=0.00002, targetCFL=0.1, silent=false)
# plotShockTubeResults_PyPlot(P, U, T, rho)
