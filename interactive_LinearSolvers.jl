include("LinearSolvers.jl")

matrix = getMatrix()
# getMatrix2 = getMatrix()
matrix2 = copy(matrix)
matrix3 = copy(matrix)

println("")
println("")

# println("TriDiag: ")
# @time x = Solve_Diag!(matrix)
# println(x)

println("Gauss Elimination: ")
@time x = Solve_GaussElim!(matrix)
println(x)
println("")

# println("Gauss Elimination: ")
# @time x = Solve_GaussElim!(matrix2)
# println(x)
# println("")

print("Gauss-Seidel: ")
@time x = Solve_GaussSeidel!(matrix2)
println("Result: $x")
println("")

print("SOR: ")
@time x = Solve_SOR!(matrix3, 1.1)
println("Result: $x")