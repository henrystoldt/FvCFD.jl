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
# P = 100000
# T = 300
# U = [ 277.7091, 6.059633, 0 ]

# Choose a mesh
# meshPath = "OFairfoilMesh"
# boundaryConditions = [ wallBoundary, [], emptyBoundary, [], supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [] ]
# meshPath = "OFmemesAirfoilMesh"
# boundaryConditions = [ wallBoundary, [], emptyBoundary, [], supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [] ]
#
# # Solve
# mesh = OpenFOAMMesh(meshPath)
# cellPrimitives = initializeUniformSolution3D(mesh, P, T, U...)
# unstructured3DFVM(mesh, meshPath, cellPrimitives, boundaryConditions, initDt=0.0000001, endTime=0.00002, outputInterval=0.00002, targetCFL=0.5, silent=false, restart=true)


### Delta Wing ###
# Freestream Conditions
P = 7669
T = 180
U = [ 539, 0, 0 ]

# Choose a mesh
meshPath = "OFmemesDeltaMesh"
boundaryConditions = [ supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [], symmetryBoundary, [], wallBoundary, [] ]

# Solve
mesh = OpenFOAMMesh(meshPath)
cellPrimitives = initializeUniformSolution3D(mesh, P, T, U...)
unstructured3DFVM(mesh, meshPath, cellPrimitives, boundaryConditions, initDt=0.00000001, endTime=0.001, outputInterval=0.00002, targetCFL=0.5, silent=false, restart=false)
