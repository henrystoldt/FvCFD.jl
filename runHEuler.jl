include("HEuler.jl")


################## Output ##################
nCells = 500
# P, U, T, rho = macCormack1DFDM(initializeShockTubeFDM(nCells)..., initDt=0.00000001, endTime=0.14267)
# P, U, T, rho = macCormack1DConservativeFDM(initializeShockTubeFDM(nCells)..., initDt=0.00001, endTime=0.1, Cx=0.3)
# P, U, T, rho = upwind1DConservativeFDM(initializeShockTubeFDM(nCells)..., initDt=0.00001, endTime=0.1, targetCFL=0.01, Cx=0.3)
# xVel = U

P, U, T, rho = upwindFVM(initializeShockTubeFVM(nCells)..., initDt=0.0000001, endTime=0.14267, targetCFL=0.1, Cx=0.3)
xVel = Array{Float64, 1}(undef, nCells)
for i in 1:nCells
    xVel[i] = U[i][1]
end

plotShockTubeResults(P, xVel, T, rho)