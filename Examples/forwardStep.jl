using FvCFD

# All file paths are relative to the repository's main directory, 'include' this script from there

### Freestream Conditions (Mach 3) ###
P = 1
T = 1
U = [ 3, 0, 0 ]

### Choose a mesh ###
# meshPath = "test/OFforwardStepMesh"
# boundaryConditions = [ supersonicInletBoundary, [P, T, U..., 1005], zeroGradientBoundary, [], symmetryBoundary, [], symmetryBoundary, [], wallBoundary, [], emptyBoundary, [] ]
meshPath = "test/OFforwardStepFineMesh"
boundaryConditions = [ supersonicInletBoundary, [P, T, U..., 2.5], zeroGradientBoundary, [], symmetryBoundary, [], symmetryBoundary, [], wallBoundary, [], emptyBoundary, [] ]

### Run ###
mesh = OpenFOAMMesh(meshPath)
cellPrimitives = initializeUniformSolution3D(mesh, P, T, U...)
solve(mesh, meshPath, cellPrimitives, boundaryConditions, ShuOsher, initDt=0.00000001, endTime=4, outputInterval=0.5, targetCFL=0.95, silent=false, restart=false, Cp=2.5, gamma=1.4, R=0.71428)