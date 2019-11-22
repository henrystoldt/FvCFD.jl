using Test
include("HEuler.jl")
include("test.jl")

@testset "Vector Utilities" begin
    v1 = (1, 2, 3)
    @test mag(v1) == 3.7416573867739413
    v2 = (2, 3, 4)
    @test dot(v1, v2) == 20.0
end;

@testset "Shock Tube initialization" begin
    nCells = 4
    
    # FDM
    dx = [ 0.25, 0.25, 0.25, 0.25 ]
    P = [ 1, 1, 0.1, 0.1 ]
    T = [ 0.00348432, 0.00348432, 0.00278746, 0.00278746 ]
    U = [ 0, 0, 0, 0 ]
    dx2, P2, T2, U2 = initializeShockTubeFDM(4)
    @test almostEqual(dx2, dx)
    @test almostEqual(P2, P)
    @test almostEqual(T2, T)
    @test almostEqual(U2, U)

    # FVM
    U = [ [0.0,0.0,0.0],  [0.0,0.0,0.0],  [0.0,0.0,0.0],  [0.0,0.0,0.0] ]
    cells = [ [4,1], [1,2], [2,3], [3,5] ]
    faces = [ [1,2], [2,3], [3,4], [1,], [4,] ]
    fAVecs = [ [0.01,0,0], [0.01,0,0], [0.01,0,0], [0.01,0,0], [0.01,0,0] ]
    cVols = [ 0.0025, 0.0025, 0.0025, 0.0025 ]
    bdryFaces = [ [4,], [5,] ]
    mesh1 = [ cells, faces, fAVecs, bdryFaces, cVols ]
    mesh2, P2, T2, U2 = initializeShockTubeFVM(4)
    @test almostEqual(P2, P)
    @test almostEqual(T2, T)
    for i in 1:4
        @test size(mesh1[i], 1) == size(mesh2[i], 1)
        for a in 1:size(mesh1[i], 1)
            @test almostEqual(mesh1[i][a], mesh2[i][a])
        end
    end
    @test almostEqual(mesh1[5], mesh2[5])
end;

@testset "Linear interpolation" begin
    @test 1 == 1
end;