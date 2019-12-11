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
# Freestream Conditions (Mach 2)
# P = 100000
# T = 300
# U = [ 694.44, 0, 0 ]
#
# meshPath = "OFunstructuredWedgeMesh"
# boundaryConditions = [ emptyBoundary, [], supersonicInletBoundary, [P, T, U..., 1005], zeroGradientBoundary, [], symmetryBoundary, [], zeroGradientBoundary, [], wallBoundary, [] ]
# meshPath = "OFwedgeMesh"
# boundaryConditions = [ supersonicInletBoundary, [P, T, U..., 1005], zeroGradientBoundary, [], symmetryBoundary, [], zeroGradientBoundary, [], wallBoundary, [], emptyBoundary, [] ]
# meshPath = "OFcoarseWedgeMesh"
# boundaryConditions = [ emptyBoundary, [], supersonicInletBoundary, [P, T, U..., 1005], zeroGradientBoundary, [], symmetryBoundary, [], zeroGradientBoundary, [], wallBoundary, [] ]
# meshPath = "OFmemesWedgeMesh"
# boundaryConditions = [ symmetryBoundary, [], emptyBoundary, [], supersonicInletBoundary, [P, T, U..., 1005], wallBoundary, [], zeroGradientBoundary, [], zeroGradientBoundary, [] ]

# mesh = OpenFOAMMesh(meshPath)
# cellPrimitives1 = initializeUniformSolution3D(mesh, P, T, U...)
# cellPrimitives2 = initializeUniformSolution3D(mesh, P, T, U...)
# unstructured3DFVM(mesh, meshPath, cellPrimitives1, boundaryConditions, initDt=0.00000001, endTime=0.0005, outputInterval=0.0005, targetCFL=0.5, silent=false)
# @time unstructured3DFVM(mesh, meshPath, cellPrimitives2, boundaryConditions, initDt=0.00000001, endTime=0.005, outputInterval=0.005, targetCFL=0.5, silent=true)

### Forward Step ###
# Freestream Conditions (Mach 3)
# P = 100000
# T = 300
# U = [ 1041.66, 0, 0 ]

# meshPath = "OFforwardStepMesh"
# boundaryConditions = [ supersonicInletBoundary, [P, T, U..., 1005], zeroGradientBoundary, [], symmetryBoundary, [], symmetryBoundary, [], wallBoundary, [], emptyBoundary, [] ]
# meshPath = "OFforwardStepFineMesh"
# boundaryConditions = [ supersonicInletBoundary, [P, T, U..., 1005], zeroGradientBoundary, [], symmetryBoundary, [], symmetryBoundary, [], wallBoundary, [], emptyBoundary, [] ]

# Solve
# mesh = OpenFOAMMesh(meshPath)
# cellPrimitives = initializeUniformSolution3D(mesh, P, T, U...)
# Endtime 0.01152
# unstructured3DFVM(mesh, meshPath, cellPrimitives, boundaryConditions, ShuOsher, initDt=0.00000001, endTime=0.03, outputInterval=0.000288, targetCFL=0.95, silent=false, restart=false)

#### Simple convection ####
# mesh, cellPrimitives = initializeShockTube3DFVM(10)
# P = 100000
# T = 300
# U = [ 694.26, 0, 0 ]
# boundaryConditions = [ supersonicInletBoundary, [P, T, U..., 1005], zeroGradientBoundary, [], emptyBoundary, [] ]
# P, U, T, rho = unstructured3DFVM(mesh, initializeUniformSolution3D(mesh, P, T, U...), boundaryConditions, initDt=0.0000001, endTime=0.00002, targetCFL=0.1, silent=false)
# plotShockTubeResults_PyPlot(P, U, T, rho)
