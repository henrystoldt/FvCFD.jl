# File used to benchmark the performance of individual functions and different versions of those functions to improve solver performance
# Useful for developers readOFBoundaryFile

# Install these packages first
using Profile
using ProfileView #Doesn't want to install on my work ubuntu desktop for some reason
using BenchmarkTools

include("shockTube.jl")
include("FvCFD.jl")

# Create mesh
meshPath = "test/OFshockTube_400"
OFmesh = OpenFOAMMesh(meshPath)
nCells = size(OFmesh.cells, 1)
mesh, cellPrimitives = initializeShockTube3DFVM(nCells...)
boundaryConditions = [ zeroGradientBoundary, [], emptyBoundary, [] ]
nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(OFmesh)

# Create solution state
nVars = 5
nFluxes = 3*nVars
cellState = encodePrimitives3D(cellPrimitives)
cellFluxes = zeros(nCells, nFluxes)
fluxResiduals = zeros(nCells, nVars)
faceFluxes = zeros(nFaces, nFluxes)
sln = SolutionState(cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes)
decodeSolution_3D(sln)

# Apply zero face fluxes to empty faces
zeroFlux = zeros(nFluxes)
for b in 1:nBoundaries
    if boundaryConditions[2*b-1] == emptyBoundary
        for f in OFmesh.boundaryFaces[b]
            sln.faceFluxes[f,:] = zeroFlux
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
P = sln.cellPrimitives[:,1]
P = reshape(P, nCells, :)
function run1000()
    for i in 1:1000
        greenGaussGrad(OFmesh, P)
    end
end
# @profview run1000()
output = zeros(400, 1, 3)
@benchmark greenGaussGrad(OFmesh, P)
# output = zeros(2001, 1)
# @benchmark linInterp_3D(OFmesh, P, output)
# @code_warntype greenGaussGrad(OFmesh, P, false)


#### LinInterp ####
# @code_warntype linInterp_3D(mesh, sln.cellFluxes, sln.faceFluxes)


#### Entire Code ####
# @profview solve(OFmesh, meshPath, cellPrimitives, boundaryConditions, initDt=0.00001, endTime=0.14267, outputInterval=0.14267, targetCFL=0.1, silent=true)
# @btime solve(OFmesh, meshPath, cellPrimitives, boundaryConditions, initDt=0.00001, endTime=0.005, outputInterval=0.14267, targetCFL=0.01, silent=true, createRestartFile=false, createVTKOutput=false)
# @code_warntype solve(OFmesh, meshPath, cellPrimitives, boundaryConditions, initDt=0.00001, endTime=0.00001, outputInterval=0.14267, targetCFL=0.01, silent=true, createRestartFile=false, createVTKOutput=false)