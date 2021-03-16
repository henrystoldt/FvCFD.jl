using FvCFD

# All file paths are relative to the repository's main directory, 'include' this script from there

### Freestream Conditions ###
P = 100000
T = 300
U = [ 694.26, 0, 0 ]
Cp = 1005

### Choose a mesh ###
meshPath = "test/OFflaredConeMesh"
boundaryConditions = [ supersonicInletBoundary, [P, T, U..., Cp], zeroGradientBoundary, [], wallBoundary, [] ]

### Solve ###
mesh = OpenFOAMMesh(meshPath)
cellPrimitives = initializeUniformSolution3D(mesh, P, T, U...)
solve(mesh, meshPath, cellPrimitives, boundaryConditions, LTSEuler, initDt=0.00000001, endTime=18000, outputInterval=9.3, targetCFL=0.1, silent=false, restart=false)