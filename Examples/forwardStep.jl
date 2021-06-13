using FvCFD

# All file paths are relative to the repository's main directory, 'include' this script from there

### Freestream Conditions (Mach 3) ###
P = 1
T = 1
U = [ 3, 0, 0 ]

### Define the fluid ###
Cp = 2.5
R = 0.71428
gamma = 1.4
fluid = Fluid(Cp, R, gamma)

### Choose a mesh ###
meshPath = "test/OFforwardStepMesh"
# meshPath = "test/OFforwardStepFineMesh"
boundaryConditions = [ supersonicInletBoundary, [P, T, U...], zeroGradientBoundary, [], symmetryBoundary, [], symmetryBoundary, [], wallBoundary, [], emptyBoundary, [] ]

### Run ###
mesh = OpenFOAMMesh(meshPath)
cellPrimitives = initializeUniformSolution3D(mesh, P, T, U...)
solve(mesh, meshPath, cellPrimitives, boundaryConditions, ShuOsher, initDt=0.00000001, endTime=4, outputInterval=0.5, targetCFL=0.95, silent=false, restart=false, fluid=fluid)
