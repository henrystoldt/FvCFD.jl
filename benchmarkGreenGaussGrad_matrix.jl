using Profile
using ProfileView #Doesn't want to install on my work ubuntu desktop for some reason
using BenchmarkTools

include("shockTube.jl")
include("finiteDifference.jl")
include("finiteVolume.jl")
include("mesh.jl")

# Create mesh
meshPath = "OFShockTubeMesh"
OFmesh = OpenFOAMMesh(meshPath)
nCells = size(OFmesh.cells, 1)
mesh, cellPrimitives = initializeShockTube3DFVM(nCells...)
nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

# Create solution state
nVars = 5
nFluxes = 3*nVars
cellState = encodePrimitives3D(cellPrimitives)
cellFluxes = zeros(nCells, nFluxes)
fluxResiduals = zeros(nCells, nVars)
faceFluxes = zeros(nFaces, nFluxes)
solutionState = [ cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]
decodeSolution(solutionState)

# Create test data
fDeltas = faceDeltas(OFmesh, solutionState)

# Test
greenGaussGrad_matrix(OFmesh, fDeltas, false)
function run1000()
    for i in 1:1000
        greenGaussGrad_matrix(OFmesh, fDeltas, false)
    end
end
@benchmark greenGaussGrad_matrix(OFmesh, fDeltas, false)
