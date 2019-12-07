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
nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(OFmesh)

# Create solution state
nVars = 5
nFluxes = 3*nVars
cellState = encodePrimitives3D(cellPrimitives)
cellFluxes = zeros(nCells, nFluxes)
fluxResiduals = zeros(nCells, nVars)
faceFluxes = zeros(nFaces, nFluxes)
solutionState = [ cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]
decodeSolution(solutionState)
# Apply zero face fluxes to empty faces
zeroFlux = zeros(nFluxes)
for b in 1:nBoundaries
    if boundaryConditions[2*b-1] == emptyBoundary
        for f in OFmesh.boundaryFaces[b]
            faceFluxes[f] = zeroFlux
        end
    end
end


#### greenGaussGrad_matrix ####
# # Create test data
# fDeltas = faceDeltas(OFmesh, solutionState)
#
# # Test
# greenGaussGrad_matrix(OFmesh, fDeltas, false)
# function run1000()
#     for i in 1:1000
#         greenGaussGrad_matrix(OFmesh, fDeltas, false)
#     end
# end
# @profview run1000()
# @benchmark greenGaussGrad_matrix(OFmesh, fDeltas, false)

#### greenGaussGrad ####
P = cellPrimitives[:,1]
greenGaussGrad(OFmesh, false, P)[1]
function run1000()
    for i in 1:1000
        greenGaussGrad(OFmesh, false, P)
    end
end
# @profview run1000()
@benchmark greenGaussGrad(OFmesh, false, P)
