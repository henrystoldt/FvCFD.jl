using Printf

include("NonLinearSolvers.jl")

functions = getFunctions()
println("")
xInit = getInitEstimates()

println("")
results = solve_NonLinear!(functions, xInit, printResiduals=true)
println("")
println("Results:")
for i in 1:length(results)
    xi = results[i]
    @printf("x%.0f = %.4f\n", i, xi)
end