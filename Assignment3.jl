include("Assignment3Q1.jl")
include("Assignment3Q2.jl")
include("Assignment3Q3.jl")


############################ Q1 #################################
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
plotResults(eulerResults, eulerResults2, "Q1: Explicit Euler")

backwardEulerResults = solve_ImplicitEuler(initConditions, ddts1, 0.0025, 10000)
backwardEulerResults2 = solve_ImplicitEuler(initConditions, ddts2, 0.0025, 10000)
plotResults(backwardEulerResults, backwardEulerResults2, "Q1: Implicit Euler")

RK4Results = solve_RK4(initConditions, ddts1, 0.0025, 10000)
RK4Results2 = solve_RK4(initConditions, ddts2, 0.0025, 10000)
plotResults(RK4Results, RK4Results2, "Q1: RK4")



############################ Q2 #################################
ddts = [ dXdt, dYdt ]
initConditions = [ 1.001, 2 ]
plots = []
timeStepSizes = [ 0.01, 0.04, 0.16, 0.64 ]

for i in 1:4
    dt = timeStepSizes[i]

    eulerResults = solve_ExplicitEuler(initConditions, ddts, dt, floor(Int32, 250/dt))
    times = eulerResults[3,:]
    explicitEulerPlot = plot(times, eulerResults[1,:], label="X", title="Explicit Euler, dt=$dt", xlabel="t", bottom_margin=15mm, left_margin=10mm)
    explicitEulerPlot = plot!(times, eulerResults[2,:], label="Y")
    push!(plots, explicitEulerPlot)

    modifiedEulerResults = solve_ModifiedEuler(initConditions, ddts, dt, floor(Int32, 250/dt))
    times = modifiedEulerResults[3,:]
    modifiedEulerPlot = plot(times, modifiedEulerResults[1,:], label="X", title="Modified Euler, dt=$dt", xlabel="t", bottom_margin=15mm, left_margin=10mm)
    modifiedEulerPlot = plot!(times, modifiedEulerResults[2,:], label="Y")
    push!(plots, modifiedEulerPlot)
end

plot(plots..., layout=(4, 2), size=(1720, 880), window_title="Assignment 3 Q2", legend=false)
gui()



############################ Q3 #################################
nNodes = 800
z, T = solve_Q3(nNodes, 0, 1, 0, 0.00001, false)
eta = Array{Float64, 1}(undef, nNodes+2)
for i in 1:nNodes+2
    eta[i] = ((i-1)/(nNodes+1)) * 20
end

plot(eta, z, label="zeta", size=(860, 440), window_title="Assignment 3 Q3")
plot!(eta, T, label="T")
gui()