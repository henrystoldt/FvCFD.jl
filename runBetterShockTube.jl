using Plots
using Plots.PlotMeasures
#using PyPlot
using LaTeXStrings
using Profile
#using ProfileView

pyplot()

#plotly()

include("shockTube.jl")
#include("finiteDifference.jl")
include("finiteVolumeRoe.jl")


################## Output ##################
nCells = 500
nCells = 15

#### FDM or Structured FVM ###
# P, U, T, rho = macCormack1DFDM(initializeShockTubeFDM(nCells)..., initDt=0.000001, targetCFL=0.95, endTime=0.14267, Cx=0.1)
# P, U, T, rho = macCormack1DConservativeFDM(initializeShockTubeFDM(nCells)..., initDt=0.00001, targetCFL=0.95, endTime=0.14267, Cx=0.3)
# P, U, T, rho = upwind1DConservativeFDM(initializeShockTubeFDM(nCells)..., initDt=0.00001, endTime=0.14267, targetCFL=0.01, Cx=0.3)
#@time P, U, T, rho = upwindFVMRoe1D(initializeShockTubeFVM(nCells)..., initDt=0.00001, endTime=0.14267, targetCFL=0.2, verbose=true)
@time P, U, T, rho = upwindFVMRoe1D(initializeShockTubeFVM(nCells)..., initDt=0.00001, endTime=0.01603, targetCFL=0.2, verbose=true)

#@time P, U, T, rho = JST_Structured1DFVM(initializeShockTube_StructuredFVM(nCells)..., forwardEuler, initDt=0.00001, endTime=0.14267, targetCFL=0.5, silent=true)
#=
xVel = U[:][4]
println(xVel)
println(typeof(U))

for i in 1:nCells
    xVel[i] = U[i][1]
    #println(xVel[i])
end
=#

#print(size(U))

xVel = U
#xVel2 = U2
#xVel3 = U3
#xVel4 = U4

### Unstructured FVM ###
# @time P, U, T, rho = central_UnstructuredADFVM(initializeShockTubeFVM(nCells, silent=false)..., initDt=0.0000001, endTime=0.14267, targetCFL=0.1, Cx=0.5, silent=false)
# println("Formatting results")

# xVel = Array{Float64, 1}(undef, nCells)
# for i in 1:nCells
#     xVel[i] = U[i][1]
# end

println("Vel = ", xVel)

println("Plotting results")
plotShockTubeResults_PyPlot(P, xVel, T, rho)
#plotShockTubeResults_PyPlot(P2, xVel2, T2, rho2)
#plotShockTubeResults_PyPlot(P3, xVel3, T3, rho3)
#plotShockTubeResults_PyPlot(P4, xVel4, T4, rho4)
