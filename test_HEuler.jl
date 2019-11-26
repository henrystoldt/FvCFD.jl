using Test
include("finiteVolume.jl")
include("finiteDifference.jl")
include("shockTube.jl")
include("vectorFunctions.jl")
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
    cVols = [ 0.0025, 0.0025, 0.0025, 0.0025 ]
    cCenters = [ [0.125,0.0,0.0], [0.375,0.0,0.0], [0.625,0.0,0.0], [0.875,0.0,0.0] ]
    faces = [ [1,2], [2,3], [3,4], [-1,1], [4,-1] ]
    fAVecs = [ [0.01,0,0], [0.01,0,0], [0.01,0,0], [0.01,0,0], [0.01,0,0] ]
    fCenters = [ [0.25,0.0,0.0], [0.5,0.0,0.0], [0.75,0.0,0.0], [0.0,0.0,0.0], [1.0, 0.0, 0.0] ]
    boundaryFaces = [ [4,], [5,] ]
    mesh1 = [ cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces ]
    
    mesh2, P2, T2, U2 = initializeShockTubeFVM(4)
    @test almostEqual(P2, P)
    @test almostEqual(T2, T)
    for i in [ 1,3,4,5,6,7 ]
        @test size(mesh1[i], 1) == size(mesh2[i], 1)
        for a in 1:size(mesh1[i], 1)
            @test almostEqual(mesh1[i][a], mesh2[i][a])
        end
    end
    @test almostEqual(mesh1[2], mesh2[2])
end;

@testset "Decode Primitives" begin
    rho = 2
    xMom = 1
    eV2 = 3
    U = 0.5
    e = 1.5 - 0.125
    T = e/(1005 - 287.05)
    P = rho*287.05*T
    res1 = [ P, T, U ]
    res2 = decodePrimitives(rho, xMom, eV2, 287.05, 1005)
    for i in 1:3
        @test almostEqual(res1[i], res2[i])
    end
end;

@testset "Encode Primitives" begin
    rho = 2
    xMom = 1
    eV2 = 3
    U = 0.5
    e = 1.5 - 0.125
    T = e/(1005 - 287.05)
    P = rho*287.05*T
    res1 = [ rho, xMom, eV2 ]
    res2 = encodePrimitives(P, T, U, 287.05, 1005)
    for i in 1:3
        @test almostEqual(res1[i], res2[i])
    end
end;

@testset "Green-Gauss Gradient" begin
    mesh, P, T, U = initializeShockTubeFVM(4)
    faceVals = [ 1, 2, 3, 4, 5 ]
    gradient = [ [0.0,0.0,0.0], [4.0,0.0,0.0], [4.0,0.0,0.0], [0.0,0.0,0.0] ]
    grad2 = greenGaussGrad(mesh, faceVals)[1]
    for i in 1:4
        @test almostEqual(gradient[i], grad2[i])
    end
end;

@testset "Linear Interpolation" begin
    mesh, P, T, U = initializeShockTubeFVM(4)
    cellVals = [ 1, 2, 3, 4 ]
    faceVals = [ 1.5, 2.5, 3.5, 0, 0 ]
    faceVals2 = linInterp(mesh, cellVals)[1]
    @test almostEqual(faceVals, faceVals2)
end;

@testset "Upwind Interpolation" begin
    mesh, P, T, U = initializeShockTubeFVM(4)
    cellVel = [ [-1,0,0], [-1,0,0], [-1,0,0], [-1,0,0] ]
    cellVals = [ 1, 2, 3, 4 ]
    faceVals = [ 2, 3, 4, 0, 0 ]
    faceVals2 = upwindInterp(mesh, cellVel, cellVals)[1]
    @test almostEqual(faceVals, faceVals2)

    cellVel = [ [1,0,0], [1,0,0], [1,0,0], [1,0,0] ]
    faceVals = [ 1, 2, 3, 0, 0 ]
    faceVals2 = upwindInterp(mesh, cellVel, cellVals)[1]
    @test almostEqual(faceVals, faceVals2)

    cellVel = [ [1,0,0], [1,0,0], [-1,0,0], [-1,0,0] ]
    faceVals = [ 1, 2.5, 4, 0, 0 ]
    faceVals2 = upwindInterp(mesh, cellVel, cellVals)[1]
    @test almostEqual(faceVals, faceVals2)
end;

@testset "Upwind FVM" begin
    nCells = 4
    P = [ 1, 0.99855553, 0.087046989, 0.1 ]
    U = [ 0, 0.090015665, 0.720126353, 0 ]
    T = [ 0.00348432, 0.00347929, 0.0024263, 0.00278746 ]
    
    P2, U2, T2, rho2 = upwindFVM(initializeShockTubeFVM(nCells, Pratio=0.1)..., initDt=0.051, endTime=0.05, Cx=0)
    xVel = []
    for i in 1:4
        push!(xVel, U2[i][1])
    end

    @test almostEqual(P, P2, 6)
    @test almostEqual(U, xVel, 6)
    @test almostEqual(T, T2, 6)
end;