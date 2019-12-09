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
gamma = 1.4
R = 287.05
Cp = 1005

Umag = mag(U)
a = sqrt(gamma*R*T)
machNum = Umag/a
Pt = P*(1 + ((gamma-1)/2)*machNum^2)^(gamma/(gamma-1))
Tt = T*(1 + ((gamma-1)/2)*machNum^2)
UunitVec = normalize(U)

# Choose a mesh
meshPath = "OFairfoilMesh"
boundaryConditions = [ wallBoundary, [], emptyBoundary, [], subsonicInletBoundary, [Pt, Tt, UunitVec..., gamma, R, Cp], pressureOutletBoundary, P ]
# meshPath = "OFmemesAirfoilMesh"
# boundaryConditions = [ wallBoundary, [], emptyBoundary, [], subsonicInletBoundary, [Pt, Tt, UunitVec..., gamma, R, Cp], pressureOutletBoundary, P ]
#
# # Solve
# mesh = OpenFOAMMesh(meshPath)
# cellPrimitives = initializeUniformSolution3D(mesh, P, T, U...)
# unstructured3DFVM(mesh, meshPath, cellPrimitives, boundaryConditions, initDt=0.0000001, endTime=0.00002, outputInterval=0.00002, targetCFL=0.5, silent=false, restart=true)


# ### Delta Wing ###
# # Freestream Conditions
# P = 7669
# T = 180
# U = [ 539, 0, 0 ]
#
# # Choose a mesh
# meshPath = "OFmemesDeltaMesh"
# boundaryConditions = [ supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [], symmetryBoundary, [], wallBoundary, [] ]

### Cone Flare ###
# Freestream Conditions
# P = 100000
# T = 300
# U = [ 694.26, 0, 0 ]
#
# # Choose a mesh
# meshPath = "OFflaredConeMesh"
# boundaryConditions = [ supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [], wallBoundary, [] ]

# Solve
mesh = OpenFOAMMesh(meshPath)
cellPrimitives = initializeUniformSolution3D(mesh, P, T, U...)
unstructured3DFVM(mesh, meshPath, cellPrimitives, boundaryConditions, initDt=0.00000001, endTime=0.015, outputInterval=0.0002, targetCFL=0.5, silent=false, restart=false)
