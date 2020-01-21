using Test
include("LinearSolvers.jl")

# For numbers of a magnitude less than ~1000
function almostEqual(iterable1, iterable2, allowedError=0.0000000000001)
    if size(iterable1, 1) != size(iterable2, 1)
        return false
    end

    for i in 1:size(iterable1, 1)
        if abs(iterable1[i] - iterable2[i]) > allowedError
            return false
        end
    end

    return true
end

@testset "Gauss Elimination" begin
    matrix1 = Array{Float64, 2}([ 
        3 1 0 0 2 0 1;
        1 3 1 0 0 2 -1;
        0 1 3 1 0 0 2;
        0 0 1 3 1 0 0;
        2 0 0 1 3 1 3;
        0 2 0 0 1 3 -4 
    ])
    solution1 = [ 118, -127, 44, 21, -107, 103 ] / 13.0
    @test almostEqual(Solve_GaussElim!(matrix1), solution1)

    matrix2 = Array{Float64, 2}([
        -2.25 1 0 0 0 0 0 0
        1 -2.25 1 0 0 0 0 0
        0 1 -2.25 1 0 0 0 0
        0 0 1 -2.25 1 0 0 0
        0 0 0 1 -2.25 1 0 0
        0 0 0 0 1 -2.25 1 0
        0 0 0 0 0 1 -2.25 -100
    ])
    solution2 = [ 1638400/833049.0, 409600/92561.0, 6656000/833049.0, 25600/1889.0, 18745600/833049.0, 3432000/92561.0, 50752400/833049.0 ]
    @test almostEqual(Solve_GaussElim!(matrix2), solution2)

    matrix3 = Array{Float64, 2}([
        1 -1 4 0 2 9 19;
        0 5 -2 7 8 4 2;
        1 0 5 7 3 -2 13;
        6 -1 2 3 0 8 -7;
        -4 2 0 5 -5 3 -9;
        0 7 -1 5 4 -2 2;
    ])
    solution3 = [ -349097/198146.0, 88792/99073.0, 401437/99073.0, -160214/99073.0, 404597/198146.0, 30085/198146.0 ]
    @test almostEqual(Solve_GaussElim!(matrix3), solution3)

end;

@testset "Iterative Solvers" begin
    matrix2 = Array{Float64, 2}([ 
        4 -1 0 1 0 100;
        -1 4 -1 0 1 100;
        0 -1 4 -1 0 100;
        1 0 -1 4 -1 100;
        0 1 0 -1 4 100; 
    ])
    solution2 = [ 25.0, 250/7.0, 300/7.0, 250/7.0, 25 ]
    @test almostEqual(Solve_GaussSeidel!(matrix2), solution2, 0.00001)
    @test almostEqual(Solve_SOR!(matrix2, 1.1), solution2, 0.00001)

end;

@testset "Diagonal Solver" begin
    matrix3 = Array{Float64, 2}([
        0 -2.25 1 0;
        1 -2.25 1 0;
        1 -2.25 1 0;
        1 -2.25 1 0;
        1 -2.25 1 0;
        1 -2.25 1 0;
        1 -2.25 0 -100;
    ])
    solution3 = [ 1638400/833049.0, 409600/92561.0, 6656000/833049.0, 25600/1889.0, 18745600/833049.0, 3432000/92561.0, 50752400/833049.0 ]
    @test almostEqual(Solve_Diag!(matrix3), solution3)

    matrix4 = Array{Float64, 2}([
        0 0 0 0 3 1 0 0 2  1;
        0 0 0 1 3 1 0 0 2 -1;
        0 0 0 1 3 1 0 0 0 2;
        0 0 0 1 3 1 0 0 0 0;
        2 0 0 1 3 1 0 0 0 3;
        2 0 0 1 3 0 0 0 0 -4;
    ])
    solution4 = [ 118, -127, 44, 21, -107, 103 ] / 13.0
    @test almostEqual(Solve_Diag!(matrix4), solution4)
end;