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