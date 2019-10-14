include("LinearSolvers.jl")
include("NonLinearSolvers.jl")
using Plots
using Plots.PlotMeasures
plotly()

############## ddt functions defined in Assignment 3, Q1 ###############
function dxdt(x, y, z)
    return -10*x + 10*y
end

function dydt(x, y, z)
    return -x*z + 28*x - y
end

function dzdt(x, y, z)
    return x*y - (8/3)*z
end

# Modified function with b = 8/3 + 0.001
function dzdt2(x, y, z)
    return x*y - (8/3 + 0.001)*z
end

############ Solvers ##############
function solve_ExplicitEuler(vars, ddts, dt, nTimeSteps)
    #Setup
    nVars = length(vars)
    timeStepCounter = 2
    results = Array{Float64, 2}(undef, nVars+1, nTimeSteps+1)

    for i in 1:nVars
        results[i,1] = vars[i]
    end
    results[nVars+1,1] = 0

    while timeStepCounter <= nTimeSteps+1
        for i in 1:nVars
            results[i, timeStepCounter] = results[i, timeStepCounter-1] + dt*ddts[i](results[1:nVars, timeStepCounter-1]...)
        end
        results[nVars+1, timeStepCounter] = results[nVars+1, timeStepCounter-1] + dt

        timeStepCounter += 1
    end

    return results
end

function solve_ImplicitEuler(vars, ddts, dt, nTimeSteps)
    #Setup
    nVars = length(vars)
    timeStepCounter = 2
    results = Array{Float64, 2}(undef, nVars+1, nTimeSteps+1)

    for i in 1:nVars
        results[i,1] = vars[i]
    end
    results[nVars+1,1] = 0

    augmentedMatrix = Array{Float64, 2}(undef, length(ddts), nVars+1)

    while timeStepCounter <= nTimeSteps+1
        #Predict solution with explicit Euler
        initGuesses = Array{Float64, 1}(undef, nVars)
        for i in 1:nVars
            initGuesses[i] = results[i, timeStepCounter-1] + dt*ddts[i](results[1:nVars, timeStepCounter-1]...)
        end
    
        #Refine solution using Newton-Raphson iterations to get backward Euler solution
        residualTolerance = 0.00001
        maxRes = 1
        maxIter = 100
        iterCount = 0
        while maxRes > residualTolerance && iterCount < maxIter
            augmentedMatrix[:,1:nVars] = ddx(ddts, initGuesses)*dt

            for i in 1:nVars
                augmentedMatrix[i, i] -= 1
            end

            maxRes = 0
            for i in 1:nVars
                residual = -(results[i, timeStepCounter-1] + ddts[i](initGuesses...)*dt - initGuesses[i])
                if abs(residual) > maxRes
                    maxRes = abs(residual)
                end
                augmentedMatrix[i, nVars+1] = residual
            end
            initGuesses += Solve_GaussElim!(augmentedMatrix)
            iterCount += 1
        end

        if iterCount == maxIter
            println("ERROR: timestep $timeStepCounter not converged, maxRes = $maxRes")
        end

        results[1:nVars, timeStepCounter] = initGuesses
        results[nVars+1, timeStepCounter] = results[nVars+1, timeStepCounter-1] + dt

        timeStepCounter += 1
    end

    return results
end

function solve_RK4(vars, ddts, dt, nTimeSteps)
    #Setup
    nVars = length(vars)
    timeStepCounter = 2
    results = Array{Float64, 2}(undef, nVars+1, nTimeSteps+1)

    for i in 1:nVars
        results[i, 1] = vars[i]
    end
    results[nVars+1, 1] = 0

    while timeStepCounter <= nTimeSteps+1
        dy1 = Array{Float64, 1}(undef, nVars)
        for i in 1:nVars
            dy1[i] = ddts[i](results[1:nVars, timeStepCounter-1]...) * dt
        end
        
        y2 = results[1:3, timeStepCounter-1] + dy1/2
        dy2 = Array{Float64, 1}(undef, nVars)
        for i in 1:nVars
            dy2[i] = ddts[i](y2...) * dt
        end

        y3 = results[1:3, timeStepCounter-1] + dy2/2
        dy3 = Array{Float64, 1}(undef, nVars)
        for i in 1:nVars
            dy3[i] = ddts[i](y3...) * dt
        end

        y4 = results[1:3, timeStepCounter-1] + dy3
        dy4 = Array{Float64, 1}(undef, nVars)
        for i in 1:nVars
            dy4[i] = ddts[i](y4...) * dt
        end

        for i in 1:nVars
            results[i, timeStepCounter] = results[i, timeStepCounter-1] + (1/6)*(dy1[i] + 2*dy2[i] + 2*dy3[i] + dy4[i])
        end

        results[nVars+1, timeStepCounter] = results[nVars+1, timeStepCounter-1] + dt
        timeStepCounter += 1
    end

    return results
end

function showQ1Plots()
    function plotResultsVsT(results, plotTitle)
        times = results[4, :]
        a = plot(times, results[1,:], label="x", title=plotTitle, size=(1720, 880), legend=true, xlabel="t", bottom_margin=15mm, left_margin=10mm)
        a = plot!(times, results[2,:], label="y")
        a = plot!(times, results[3,:], label="z")
        return a
    end

    function plotResults(res, res2, title)
        p1 = plotResultsVsT(res, "b=8/3")
        p2 = plotResultsVsT(res2, "b=8/3+0.001")

        pxy = plot(res[1,:], res[2,:], title="X vs Y, b=8/3", xlabel="x", ylabel="y", bottom_margin=15mm, left_margin=10mm)
        pyz = plot(res[2,:], res[3,:], title="Y vs Z, b=8/3", xlabel="y", ylabel="z", bottom_margin=15mm, left_margin=10mm)
        pxz = plot(res[1,:], res[3,:], title="X vs Z, b=8/3", xlabel="x", ylabel="z", bottom_margin=15mm, left_margin=10mm)

        pxy2 = plot(res2[1,:], res2[2,:], title="X vs Y, b=8/3+0.001", xlabel="x", ylabel="y", bottom_margin=15mm, left_margin=10mm)
        pyz2 = plot(res2[2,:], res2[3,:], title="Y vs Z, b=8/3+0.001", xlabel="y", ylabel="z", bottom_margin=15mm, left_margin=10mm)
        pxz2 = plot(res2[1,:], res2[3,:], title="X vs Z, b=8/3+0.001", xlabel="x", ylabel="z", bottom_margin=15mm, left_margin=10mm)

        plot(p1, p2, pxy, pxy2, pyz, pyz2, pxz, pxz2, layout=(2, 4), size=(1720, 880), window_title=title, legend=false, plot_title="Explicit Euler")
        gui()
    end

    initConditions = [ 0.001, 0.001, 0.002 ]
    # Original time derivative functions
    ddts1 = [ dxdt, dydt, dzdt ]
    # Modified set, with b = 8/3 + 0.001
    ddts2 = [ dxdt, dydt, dzdt2 ]

    eulerResults = solve_ExplicitEuler(initConditions, ddts1, 0.0025, 10000)
    eulerResults2 = solve_ExplicitEuler(initConditions, ddts2, 0.0025, 10000)
    plotResults(eulerResults, eulerResults2, "Explicit Euler")

    backwardEulerResults = solve_ImplicitEuler(initConditions, ddts1, 0.0025, 10000)
    backwardEulerResults2 = solve_ImplicitEuler(initConditions, ddts2, 0.0025, 10000)
    plotResults(backwardEulerResults, backwardEulerResults2, "Implicit Euler")

    RK4Results = solve_RK4(initConditions, ddts1, 0.0025, 10000)
    RK4Results2 = solve_RK4(initConditions, ddts2, 0.0025, 10000)
    plotResults(RK4Results, RK4Results2, "RK4")
end

showQ1Plots()