using Test
include("NonLinearSolvers.jl")
include("test.jl")

@testset "Functions 1, 2, 3" begin
    x = [1,1,1]
    res1 = exp(2) - 1 + 4
    @test res1 == f1(x...)

    res2 = 1 - 1 - 1
    @test res2 == f2(x...)

    res3 = 1 - sin(1)
    @test res3 == f3(x...)
end

@testset "Partial Derivatives" begin
    x = [1.0, 1.0, 1.0]
    
    res1 = Array{Float64, 2}(undef, 1, 3)
    res1[1,1] = 2*exp(2)
    res1[1,2] = -1
    res1[1,3] = 0
    ddxs = ddx([f1], x)
    @test almostEqual(ddxs, res1, 0.0000001)

    res2 = Array{Float64, 2}(undef, 1, 3)
    res2[1,1] = 0
    res2[1,2] = 1
    res2[1,3] = -2
    ddxs = ddx([f2], x)
    @test almostEqual(ddxs, res2, 0.0000001)

    res3 = Array{Float64, 2}(undef, 1, 3)
    res3[1,1] = -0.54030230586
    res3[1,2] = 0
    res3[1,3] = 1
    ddxs = ddx([f3], x)
    @test almostEqual(ddxs, res3, 0.0000001)

    allResults = Array{Float64, 2}(undef, 3, 3)
    res = [ res1, res2, res3 ]
    for i in 1:3
        for a in 1:3
            allResults[i,a] = res[i][a]
        end
    end 
    ddxs = ddx([f1, f2, f3], x)
    @test almostEqual(ddxs, allResults, 0.0000001)
end

@testset "Calc residuals" begin
    result = calcResiduals([f1, f2, f3], [1.0, 1.0, 1.0])
    res1 = exp(2) - 1 + 4
    res2 = 1 - 1 - 1
    res3 = 1 - sin(1)
    solution = [ res1, res2, res3 ]
    @test almostEqual(result, solution)
end 

@testset "Nonlinear Solvers" begin
    function fTest1(x2, x3)
        return 6*cos(x2) + 8*cos(x3) + 4*cos(3.8397) - 10
    end
    function gTest1(x2, x3)
        return 6*sin(x2) + 8*sin(x3) + 4*sin(3.8397)
    end    
    result = solve_NonLinear!([fTest1, gTest1], [0.523598776, 0.0 ])
    solution = [ 0.558770302, -0.076288115 ]
    @test almostEqual(result, solution, 0.0001)
end