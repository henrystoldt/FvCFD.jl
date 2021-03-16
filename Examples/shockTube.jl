using FvCFD
include("../test/shockTubeFunctions.jl")

### Select a mesh ###
# meshPath = "test/OFshockTube_100"
meshPath = "test/OFshockTube_400"

# Load it
OFmesh = OpenFOAMMesh(meshPath)
nCells = size(OFmesh.cells, 1)

### Initialize non-uniform initial conditions ###
_, cellPrimitives = initializeShockTube3DFVM(nCells...)

### Run ###
boundaryConditions = [ zeroGradientBoundary, [], emptyBoundary, [] ]
cellPrimitives = solve(OFmesh, meshPath, cellPrimitives, boundaryConditions, ShuOsher, initDt=0.0001, endTime=0.14267, outputInterval=0.14267, targetCFL=0.5, silent=false, createVTKOutput=true)