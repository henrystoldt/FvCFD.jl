include("LinearSolvers.jl")

matrix = getMatrix()

println("")
println("")

println("Enter number to choose a solver:
1. Gauss Elimination
2. Gauss-Seidel
3. SOR (omega=1.5)
4. Generalized Thomas Algorithm")

choice = parse(Int32, chomp(readline()))

if choice == 1
    println("Gauss Elimination: ")
    @time x = Solve_GaussElim!(matrix)
    println("x = $x")
elseif choice == 2
    print("Gauss-Seidel: ")
    @time x = Solve_GaussSeidel!(matrix)
    println("x = $x")
elseif choice == 3
    print("SOR: ")
    @time x = Solve_SOR!(matrix, 1.5)
    println("x = $x")
else
    println("TriDiag: ")
    @time x = Solve_Diag!(matrix)
    println("x = $x")
end