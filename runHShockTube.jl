using Plots
using Plots.PlotMeasures
using LaTeXStrings
using Profile
using ProfileView #Doesn't want to install on my work ubuntu desktop for some reason
using BenchmarkTools

pyplot()

include("shockTube.jl")
include("finiteDifference.jl")
include("finiteVolume.jl")
include("mesh.jl")

################## Output ##################

println("Meshing")
#### FDM or Structured FVM ###
# P, U, T, rho = macCormack1DFDM(initializeShockTubeFDM(nCells)..., initDt=0.000001, targetCFL=0.05, endTime=0.14267, Cx=0.5)
# P, T, U, rho = macCormack1DConservativeFDM(initializeShockTubeFDM(nCells)..., initDt=0.00001, targetCFL=0.1, endTime=0.14267, Cx=0.5)
# P, U, T, rho = upwind1DConservativeFDM(initializeShockTubeFDM(nCells)..., initDt=0.00001, endTime=0.14267, targetCFL=0.01, Cx=0.3)
# P, U, T, rho = structured1DFVM(initializeShockTube_StructuredFVM(nCells)..., forwardEuler, initDt=0.00001, endTime=0.14267, targetCFL=0.5, silent=true)
# P, U, T, rho = structured1DFVM(initializeShockTube_StructuredFVM(nCells)..., RK2_Mid, initDt=0.00001, endTime=0.14267, targetCFL=0.5, silent=false)
# P, U, T, rho = structured1DFVM(initializeShockTube_StructuredFVM(nCells)..., RK4, initDt=0.00001, endTime=0.14267, targetCFL=0.5, silent=true)
# P, U, T, rho = structured1DFVM(initializeShockTube_StructuredFVM(nCells)..., ShuOsher, initDt=0.00001, endTime=0.14267, targetCFL=0.5, silent=false)
# cellPrimitives = unstructured3DFVM(initializeShockTube3DFVM(nCells)..., ShuOsher, initDt=0.00001, endTime=0.14267, targetCFL=0.05, silent=false)
# xVel = U

### UnstructuredFVM from OpenFOAM Meshes ###
meshPath = "OFShockTubeMesh"
OFmesh = OpenFOAMMesh(meshPath)
nCells = size(OFmesh.cells, 1)
mesh, cellPrimitives = initializeShockTube3DFVM(nCells...)
# # Boundaries 1 and 2 are the ends, 3 is all the sides
boundaryConditions = [ zeroGradientBoundary, [], emptyBoundary, [] ]
cellPrimitives = unstructured3DFVM(OFmesh, meshPath, cellPrimitives, boundaryConditions, RK2_Mid, initDt=0.00001, endTime=0.14267, outputInterval=0.14267, targetCFL=0.3, silent=false)
# @profview unstructured3DFVM(OFmesh, meshPath, cellPrimitives, boundaryConditions, initDt=0.00001, endTime=0.14267, outputInterval=0.14267, targetCFL=0.1, silent=true)
# @btime unstructured3DFVM(OFmesh, meshPath, cellPrimitives, boundaryConditions, initDt=0.00001, endTime=0.05, outputInterval=0.14267, targetCFL=0.1, silent=true)

### Unstructured FVM ##
# @time cellPrimitives central_UnstructuredADFVM(initializeShockTubeFVM(nCells, silent=false)..., initDt=0.0000001, endTime=0.14267, targetCFL=0.1, Cx=0.5, silent=false)
# println("Formatting results")

# xVel = Array{Float64, 1}(undef, nCells)
# for i in 1:nCells
#     xVel[i] = U[i][1]
# end

P = cellPrimitives[:,1]
T = cellPrimitives[:,2]
xVel = cellPrimitives[:,3]
rho = zeros(nCells)
for i in 1:nCells
    rho[i] = idealGasRho(T[i], P[i])
end

println("Plotting results")
plotShockTubeResults_PyPlot(P, xVel, T, rho)
