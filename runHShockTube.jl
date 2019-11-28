using Plots
using Plots.PlotMeasures
using LaTeXStrings

pyplot()

include("shockTube.jl")
include("finiteDifference.jl")
include("finiteVolume.jl")


################## Output ##################
nCells = 500

#### FDM ###
# P, U, T, rho = macCormack1DFDM(initializeShockTubeFDM(nCells)..., initDt=0.00000001, endTime=0.14267)
# P, U, T, rho = macCormack1DConservativeFDM(initializeShockTubeFDM(nCells)..., initDt=0.00001, endTime=0.1, Cx=0.3)
# P, U, T, rho = upwind1DConservativeFDM(initializeShockTubeFDM(nCells)..., initDt=0.00001, endTime=0.1, targetCFL=0.01, Cx=0.3)

# xVel = U

### FVM ###
P, U, T, rho = upwindFVM(initializeShockTubeFVM(nCells, silent=false)..., initDt=0.0000001, endTime=0.14267, targetCFL=0.1, Cx=0.5, silent=false)
println("Formatting results")

xVel = Array{Float64, 1}(undef, nCells)
for i in 1:nCells
    xVel[i] = U[i][1]
end

println("Plotting results")
plotShockTubeResults_PyPlot(P, xVel, T, rho)