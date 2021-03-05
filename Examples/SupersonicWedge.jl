using fvCFD

### Freestream Conditions (Mach 2) ###
P = 100000
T = 300
U = [ 694.44, 0, 0 ]

### Select a mesh ###
meshPath = "../test/OFunstructuredWedgeMesh"
boundaryConditions = [ emptyBoundary, [], supersonicInletBoundary, [P, T, U..., 1005], zeroGradientBoundary, [], symmetryBoundary, [], zeroGradientBoundary, [], wallBoundary, [] ]
# meshPath = "../test/OFwedgeMesh"
# boundaryConditions = [ supersonicInletBoundary, [P, T, U..., 1005], zeroGradientBoundary, [], symmetryBoundary, [], zeroGradientBoundary, [], wallBoundary, [], emptyBoundary, [] ]
# meshPath = "../test/OFcoarseWedgeMesh"
# boundaryConditions = [ emptyBoundary, [], supersonicInletBoundary, [P, T, U..., 1005], zeroGradientBoundary, [], symmetryBoundary, [], zeroGradientBoundary, [], wallBoundary, [] ]
# meshPath = "../test/OFmemesWedgeMesh"
# boundaryConditions = [ symmetryBoundary, [], emptyBoundary, [], supersonicInletBoundary, [P, T, U..., 1005], wallBoundary, [], zeroGradientBoundary, [], zeroGradientBoundary, [] ]

### Run ###
mesh = OpenFOAMMesh(meshPath)
cellPrimitives = initializeUniformSolution3D(mesh, P, T, U...)

solve(mesh, meshPath, cellPrimitives, boundaryConditions, initDt=0.00000001, endTime=0.001, outputInterval=0.005, targetCFL=0.5, silent=false)