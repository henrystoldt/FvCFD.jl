include("ODESolvers.jl")
using Plots
using Plots.PlotMeasures
plotly()

########## Functions specific to assignment Question 2 ###########
function dXdt(X, Y)
    return 1 - (2 + 1)*X + X*X*Y
end

function dYdt(X, Y)
    return 2*X - X*X*Y
end

############### Script to drive above code for Q2, output results ###############
function showQ2Plots()
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
end