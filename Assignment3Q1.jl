include("ODESolvers.jl")
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

########### Script to drive above solvers and display output ##############
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