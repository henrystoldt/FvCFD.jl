include("Assignment3Q1.jl")
using Plots
using Plots.PlotMeasures
plotly()

function dXdt(X, Y)
    return 1 - (2 + 1)*X + X*X*Y
end

function dYdt(X, Y)
    return 2*X - X*X*Y
end

function solve_ModifiedEuler(vars, ddts, dt, nTimeSteps)
    #Setup
    nVars = length(vars)
    timeStepCounter = 2
    results = Array{Float64, 2}(undef, nVars+1, nTimeSteps+1)

    for i in 1:nVars
        results[i,1] = vars[i]
    end
    results[nVars+1,1] = 0

    predictor1 = Array{Float64, 1}(undef, nVars)
    predictor2 = Array{Float64, 1}(undef, nVars)
    while timeStepCounter <= nTimeSteps+1
        for i in 1:nVars
            predictor1[i] = results[i, timeStepCounter-1] + dt*ddts[i](results[1:nVars, timeStepCounter-1]...)
        end
        
        for i in 1:nVars
            predictor2[i] = results[i, timeStepCounter-1] + dt*ddts[i](predictor1...)
            results[i, timeStepCounter] = (predictor1[i] + predictor2[i]) / 2
        end
        
        results[nVars+1, timeStepCounter] = results[nVars+1, timeStepCounter-1] + dt
        timeStepCounter += 1
    end

    return results
end

ddts = [ dXdt, dYdt ]
initConditions = [ 1.001, 2 ]
plots = []
timeStepSizes = [ 0.01, 0.03, 0.09, 0.27 ]

for i in 1:4
    dt = timeStepSizes[i]

    eulerResults = solve_ExplicitEuler(initConditions, ddts, dt, 10000)
    times = eulerResults[3,:]
    explicitEulerPlot = plot(times, eulerResults[1,:], label="X", title="Explicit Euler, dt=$dt", xlabel="t", bottom_margin=15mm, left_margin=10mm)
    explicitEulerPlot = plot!(times, eulerResults[2,:], label="Y")
    push!(plots, explicitEulerPlot)

    modifiedEulerResults = solve_ModifiedEuler(initConditions, ddts, dt, 10000)
    times = modifiedEulerResults[3,:]
    modifiedEulerPlot = plot(times, modifiedEulerResults[1,:], label="X", title="Modified Euler, dt=$dt", xlabel="t", bottom_margin=15mm, left_margin=10mm)
    modifiedEulerPlot = plot!(times, modifiedEulerResults[2,:], label="Y")
    push!(plots, modifiedEulerPlot)
end


plot(plots..., layout=(4, 2), size=(1720, 880), window_title="Assignment 3 Q2", legend=false)
gui()