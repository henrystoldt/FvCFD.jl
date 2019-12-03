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
nCells = 500

println("Meshing")
#### FDM or Structured FVM ###
# cellPrimitives macCormack1DFDM(initializeShockTubeFDM(nCells)..., initDt=0.000001, targetCFL=0.95, endTime=0.14267, Cx=0.1)
# cellPrimitives macCormack1DConservativeFDM(initializeShockTubeFDM(nCells)..., initDt=0.00001, targetCFL=0.95, endTime=0.14267, Cx=0.3)
# cellPrimitives upwind1DConservativeFDM(initializeShockTubeFDM(nCells)..., initDt=0.00001, endTime=0.14267, targetCFL=0.01, Cx=0.3)
# @time cellPrimitives structured1DFVM(initializeShockTube_StructuredFVM(nCells)..., forwardEuler, initDt=0.00001, endTime=0.14267, targetCFL=0.05, silent=true)
# @time cellPrimitives structured1DFVM(initializeShockTube_StructuredFVM(nCells)..., RK2_Mid, initDt=0.00001, endTime=0.14267, targetCFL=0.5, silent=true)
# @time cellPrimitives structured1DFVM(initializeShockTube_StructuredFVM(nCells)..., RK4, initDt=0.00001, endTime=0.14267, targetCFL=0.25, silent=true)
# cellPrimitives structured1DFVM(initializeShockTube_StructuredFVM(nCells)..., ShuOsher, initDt=0.00001, endTime=0.14267, targetCFL=0.99, silent=false)
# @time cellPrimitives unstructured3DFVM(initializeShockTube3DFVM(nCells)..., ShuOsher, initDt=0.00001, endTime=0.14267, targetCFL=0.95, silent=false)
# xVel = U

### UnstructuredFVM from OpenFOAM Meshes ###
OFmesh = OpenFOAMMesh("OFShockTubeMesh")
nCells = size(OFmesh[1], 1)
mesh, cellPrimitives = initializeShockTube3DFVM(nCells...)
# Boundaries 1 and 2 are the ends, 3 is all the sides
boundaryConditions = [ zeroGradientBoundary, [], zeroGradientBoundary, [], emptyBoundary, [] ]
@time cellPrimitives unstructured3DFVM(OFmesh, cellPrimitives, boundaryConditions, ShuOsher, initDt=0.00001, endTime=0.14267, targetCFL=0.95, silent=false)

### Unstructured FVM ###
# @time cellPrimitives central_UnstructuredADFVM(initializeShockTubeFVM(nCells, silent=false)..., initDt=0.0000001, endTime=0.14267, targetCFL=0.1, Cx=0.5, silent=false)
# println("Formatting results")

# xVel = Array{Float64, 1}(undef, nCells)
# for i in 1:nCells
#     xVel[i] = U[i][1]
# end

P, T, Ux, Uy, Uz = cellPrimitives
rho = zeros(nCells)
for i in 1:nCells
    rho[i] = idealGasRho(T[i], P[i])
end

println("Plotting results")
plotShockTubeResults_PyPlot(P, xVel, T, rho)
