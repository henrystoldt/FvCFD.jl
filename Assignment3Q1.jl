include("LinearSolvers.jl")
using Plots
plotly()

function dxdt(x, y, z)
    return -10*x + 10*y
end

function dydt(x, y, z)
    return -x*z + 28*x - y
end

function dzdt(x, y, z)
    return x*y - (8/3)*z
end

function dzdt2(x, y, z)
    return x*y - (8/3 + 0.001)*z
end

function solve_ExplicitEuler(vars, ddts, dt, nTimeSteps)
    timeStepCounter = 2
    results = Array{Float64, 2}(undef, length(vars)+1, nTimeSteps+1)

    for i in 1:length(vars)
        results[i,1] = vars[i]
    end
    results[length(vars)+1,1] = 0

    while timeStepCounter <= nTimeSteps+1
        for i in 1:length(vars)
            results[i, timeStepCounter] = results[i, timeStepCounter-1] + dt*ddts[i](results[1:3, timeStepCounter-1]...)
        end
        results[length(vars)+1, timeStepCounter] = results[length(vars)+1, timeStepCounter-1] + dt

        timeStepCounter += 1
    end

    return results
end

function solve_ImplicitEuler(vars, ddts, dt, nTimeSteps)

end

function solve_RK4(vars, ddts, dt, nTimeSteps)
    timeStepCounter = 2
    results = Array{Float64, 2}(undef, length(vars)+1, nTimeSteps+1)

    for i in 1:length(vars)
        results[i, 1] = vars[i]
    end
    results[length(vars)+1, 1] = 0

    while timeStepCounter <= nTimeSteps+1
        dy1 = Array{Float64, 1}(undef, length(vars))
        for i in 1:length(vars)
            dy1[i] = ddts[i](results[1:3, timeStepCounter-1]...) * dt
        end
        
        y2 = results[1:3, timeStepCounter-1] + dy1/2
        dy2 = Array{Float64, 1}(undef, length(vars))
        for i in 1:length(vars)
            dy2[i] = ddts[i](y2...) * dt
        end

        y3 = results[1:3, timeStepCounter-1] + dy2/2
        dy3 = Array{Float64, 1}(undef, length(vars))
        for i in 1:length(vars)
            dy3[i] = ddts[i](y3...) * dt
        end

        y4 = results[1:3, timeStepCounter-1] + dy3
        dy4 = Array{Float64, 1}(undef, length(vars))
        for i in 1:length(vars)
            dy4[i] = ddts[i](y4...) * dt
        end

        for i in 1:length(vars)
            results[i, timeStepCounter] = results[i, timeStepCounter-1] + (1/6)*(dy1[i] + 2*dy2[i] + 2*dy3[i] + dy4[i])
        end

        results[length(vars)+1, timeStepCounter] = results[length(vars)+1, timeStepCounter-1] + dt
        timeStepCounter += 1
    end

    return results
end

function plotResults(results, plotTitle)
    times = results[4, :]
    plot(times, results[1,:], label="x", title=plotTitle, size=(1720, 880))
    plot!(times, results[2,:], label="y")
    plot!(times, results[3,:], label="z")
    gui()
end

initConditions = [ 0.001, 0.001, 0.002 ]
ddts1 = [ dxdt, dydt, dzdt ]
ddts2 = [ dxdt, dydt, dzdt2 ]

eulerResults = solve_ExplicitEuler(initConditions, ddts1, 0.0025, 10000)
plotResults(eulerResults, "Explicit Euler 1")

eulerResults2 = solve_ExplicitEuler(initConditions, ddts2, 0.0025, 10000)
plotResults(eulerResults2, "Explicit Euler 2")

RK4Results = solve_RK4(initConditions, ddts1, 0.0025, 10000)
plotResults(RK4Results, "RK4 1")

RK4Results2 = solve_RK4(initConditions, ddts2, 0.0025, 10000)
plotResults(RK4Results2, "RK4 2")

