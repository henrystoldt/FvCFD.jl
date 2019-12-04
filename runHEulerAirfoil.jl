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

### Airfoil ###
# Freestream Conditions
P = 100000
T = 300
U = [ 277.7091, 6.059633, 0 ]

# Choose a mesh
# meshPath = "OFairfoilMesh"
# boundaryConditions = [ wallBoundary, [], emptyBoundary, [], supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [] ]
meshPath = "OFmemesAirfoilMesh"
boundaryConditions = [ wallBoundary, [], emptyBoundary, [], supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [] ]

# Solve
mesh = OpenFOAMMesh(meshPath)
cellPrimitives = initializeUniformSolution3D(mesh, P, T, U...)
unstructured3DFVM(mesh, meshPath, cellPrimitives, boundaryConditions, initDt=0.0000001, endTime=0.00002, outputInterval=0.00002, targetCFL=0.5, silent=false, restart=true)
# @time unstructured3DFVM(mesh, meshPath, cellPrimitives, boundaryConditions, initDt=0.0000001, endTime=0.0002, outputInterval=0.0002, targetCFL=0.5, silent=false, restart=false)


#### Simple convection ####
# mesh, cellPrimitives = initializeShockTube3DFVM(10)
# P = 100000
# T = 300
# U = [ 694.26, 0, 0 ]
# boundaryConditions = [ supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [], emptyBoundary, [] ]
# P, U, T, rho = unstructured3DFVM(mesh, initializeUniformSolution3D(mesh, P, T, U...), boundaryConditions, initDt=0.0000001, endTime=0.00002, targetCFL=0.1, silent=false)
# plotShockTubeResults_PyPlot(P, U, T, rho)
