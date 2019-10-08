include("NonLinearSolvers.jl")

functions = getFunctions()
xInit = getInitEstimates()


results = solve_NonLinear!(functions, xInit, printResiduals=true)
println("")
println("Results: $results")