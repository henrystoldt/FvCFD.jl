using Plots
using Plots.PlotMeasures
using LaTeXStrings
using Profile
# using ProfileView #Doesn't want to install on my work ubuntu desktop for some reason
using BenchmarkTools

pyplot()

include("shockTube.jl")
include("finiteVolume.jl")
include("mesh.jl")

################## Output ##################
println("Meshing")

### UnstructuredFVM from OpenFOAM Meshes ###
meshPath = "OFshockTube_100"
OFmesh = OpenFOAMMesh(meshPath)
nCells = size(OFmesh.cells, 1)

_, cellPrimitives = initializeShockTube3DFVM(nCells...)

boundaryConditions = [ zeroGradientBoundary, [], emptyBoundary, [] ]
cellPrimitives = unstructured3DFVM(OFmesh, meshPath, cellPrimitives, boundaryConditions, LTSEuler, initDt=0.0001, endTime=0.14267, outputInterval=0.14267, targetCFL=0.8, silent=false, createVTKOutput=false)
# @profview unstructured3DFVM(OFmesh, meshPath, cellPrimitives, boundaryConditions, initDt=0.00001, endTime=0.14267, outputInterval=0.14267, targetCFL=0.1, silent=true)
# @btime unstructured3DFVM(OFmesh, meshPath, cellPrimitives, boundaryConditions, initDt=0.00001, endTime=0.005, outputInterval=0.14267, targetCFL=0.01, silent=true, createRestartFile=false, createVTKOutput=false)

P = cellPrimitives[:,1]
T = cellPrimitives[:,2]
xVel = cellPrimitives[:,3]
rho = zeros(nCells)
for i in 1:nCells
    rho[i] = idealGasRho(T[i], P[i])
end

println("Plotting results")
plotShockTubeResults_PyPlot(P, xVel, T, rho)
