include("LinearSolvers.jl")
include("NonLinearSolvers.jl")
using Printf

################ Q1 #################
println("Question 1: ")
matrix1 = Array{Float64, 2}([
1 -1 4 0 2 9 19;
0 5 -2 7 8 4 2;
1 0 5 7 3 -2 13;
6 -1 2 3 0 8 -7;
-4 2 0 5 -5 3 -9;
0 7 -1 5 4 -2 2; ])

matrix1a = copy(matrix1)
matrix1b = copy(matrix1)

println("Gauss Elimination: ")
x = Solve_GaussElim!(matrix1)
println("x = $x")

println("Gauss-Seidel: ")
x = Solve_GaussSeidel!(matrix1a)
println("x = $x")

println("SOR: ")
x = Solve_SOR!(matrix1b, 1.5)
println("x = $x")
println("")

################ Q2 #################
println("Question 2: ")
matrix2Diag = Array{Float64, 2}([
0 0 0 0 3 1 0 0 2  1;
0 0 0 1 3 1 0 0 2 -1;
0 0 0 1 3 1 0 0 0 2;
0 0 0 1 3 1 0 0 0 0;
2 0 0 1 3 1 0 0 0 3;
2 0 0 1 3 0 0 0 0 -4;
])

matrix2 = Array{Float64, 2}([
3 1 0 0 2 0 1;
1 3 1 0 0 2 -1;
0 1 3 1 0 0 2;
0 0 1 3 1 0 0;
2 0 0 1 3 1 3;
0 2 0 0 1 3 -4;
])

println("N-Diagonal Solver:")
x = Solve_Diag!(matrix2Diag)
println("x = $x")

println("Gauss Elimination:")
x = Solve_GaussElim!(matrix2)
println("x = $x")
println("")

################ Q3 #################
println("Question 3: ")
# Functions f1, f2, f3 from the assignment are defined in NonLinearSolver.jl

x = solve_NonLinear!([f1, f2, f3], [0.0, 0.0, 0.0])
for i in 1:length(x)
    xi = x[i]
    @printf("x%.0f = %.4f\n", i, xi)
end
