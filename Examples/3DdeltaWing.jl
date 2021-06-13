using FvCFD

### Freestream Conditions ###
P = 7669
T = 180
U = [ 0, 0, 539 ]
Cp = 1005
R = 287.05
gamma = 1.4
fluid = Fluid(Cp, R, gamma)

### Choose a mesh ###
meshPath = "test/OFdeltaMesh"

boundaryConditions = [ 
    supersonicInletBoundary, [P, T, U...],
    zeroGradientBoundary, [], 
    zeroGradientBoundary, [], 
    symmetryBoundary, [], 
    zeroGradientBoundary, [], 
    wallBoundary, [] 
]

### Solve ###
mesh = OpenFOAMMesh(meshPath)
cellPrimitives = initializeUniformSolution3D(mesh, P, T, U...)
solve(mesh, meshPath, cellPrimitives, boundaryConditions, LTSEuler, initDt=0.00000001, endTime=3000, outputInterval=10, targetCFL=0.3, silent=false, restart=false, fluid=fluid)