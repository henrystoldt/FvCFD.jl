include("fvCFD.jl")
include("shockTube.jl")

###################################################### FORMERLY RUNAIRFOIL ###############################################################
# Freestream Conditions
# P = 100000
# T = 300
# U = [ 277.7091, 6.059633, 0 ]
# UunitVec = normalize(U)
# gamma = 1.4
# R = 287.05
# Cp = 1005

# a = sqrt(gamma*R*T)
# machNum = mag(U)/a
# Pt = P*(1 + ((gamma-1)/2)*machNum^2)^(gamma/(gamma-1))
# Tt = T*(1 + ((gamma-1)/2)*machNum^2)

# Choose a mesh
# meshPath = "../test/OFairfoilMesh"
# meshPath = "../test/OFmemesAirfoilMesh"
# boundaryConditions = [ wallBoundary, [], emptyBoundary, [], subsonicInletBoundary, [Pt, Tt, UunitVec..., gamma, R, Cp], pressureOutletBoundary, P ]
#
# Solve
# mesh = OpenFOAMMesh(meshPath)
# cellPrimitives = initializeUniformSolution3D(mesh, P, T, U...)
# solve(mesh, meshPath, cellPrimitives, boundaryConditions, LTSEuler, initDt=0.0000001, endTime=100, outputInterval=25, targetCFL=0.5, silent=false)
# solve(mesh, meshPath, cellPrimitives, boundaryConditions, initDt=0.0000001, endTime=0.00008, outputInterval=0.00002, targetCFL=0.5, silent=false)


# ### Delta Wing ###
# # Freestream Conditions
# P = 7669
# T = 180
# U = [ 539, 0, 0 ]
#
# # Choose a mesh
# meshPath = "../test/OFmemesDeltaMesh"
# boundaryConditions = [ supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [], symmetryBoundary, [], wallBoundary, [] ]

### Cone Flare ###
# Freestream Conditions
# P = 100000
# T = 300
# U = [ 694.26, 0, 0 ]
#
# # Choose a mesh
# meshPath = ../test/OFflaredConeMesh"
# boundaryConditions = [ supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [], wallBoundary, [] ]

# Solve
# mesh = OpenFOAMMesh(meshPath)
# cellPrimitives = initializeUniformSolution3D(mesh, P, T, U...)
# solve(mesh, meshPath, cellPrimitives, boundaryConditions, LTSEuler, initDt=0.00000001, endTime=18000, outputInterval=9.3, targetCFL=0.1, silent=false, restart=false)
