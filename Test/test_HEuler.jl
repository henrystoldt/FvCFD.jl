using Test
include("../mesh.jl")
include("../finiteVolume.jl")
include("../1D/finiteDifference.jl")
include("../shockTube.jl")
include("../vectorFunctions.jl")
include("test.jl")
include("testMeshes.jl")

@testset "Vector Utilities" begin
    v1 = [1, 2, 3]
    @test mag(v1) == 3.7416573867739413

    v2 = [2, 3, 4]
    @test dot(v1, v2) == 20.0

    @test almostEqual(normalize(v1), [0.267261242, 0.534522484, 0.801783726], 1e-9)

    v1 = [ 1.0, 2.0, 3.0 ]
    v2 = [ 2.0, 2.0, 2.0 ]
    crossProduct = [ -2.0, 4.0, -2.0 ]
    crossProduct2 = crossProd(v1, v2)
    @test almostEqual(crossProduct, crossProduct2)
end;

@testset "Shock Tube Initialization" begin
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
    #### Simple uniform gradient test case ####
    mesh, cellPrimitives = initializeShockTube3DFVM(4)
    faceVals = zeros(5,1)
    faceVals[:,1] = [ 2, 3, 4, 1, 5 ] # Faces are numbered such that boundary faces come last, spatially values increase from one to five from left to right
    expectedGradient = [ [4.0,0.0,0.0], [4.0,0.0,0.0], [4.0,0.0,0.0], [4.0,0.0,0.0] ] # Expect uniform gradient
    gradient = greenGaussGrad(mesh, faceVals, true)
   
    for i in 1:4
        @test almostEqual(expectedGradient[i], gradient[i,1,:])
    end

    ########## Case from Moukalled pg. 249 ##########
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = meshMoukalled813()
    # Adjust properties to match fudged values in Moukalled
    cCenters = [ [1.75, 2.0, 0.5], [4.25, 3.5, 0.5] ]
    cells = [ [1, 2, 3, 4, 5 ], [1], [2], [3], [4], [5] ]
    cVols = [ 8.625, 1, 1, 1, 1, 1 ]
    faces = [ [1,2], [1,3], [1,4], [1,5], [1,6] ]
    fCenters = fCenters[1:5]
    mesh = [ cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces ]

    function samplePhiDistribution(x, y, z)
        return x^2 + y^2 + x^2*y^2
    end
    phiFaceVals = []
    for i in 1:5
        push!(phiFaceVals, samplePhiDistribution(fCenters[i]...))
    end

    grad = greenGaussGrad(mesh, true, phiFaceVals)[1]
    @test almostEqual(grad[1], [ 20.63111, 17.31454, 0.0 ], 1e-5)
end;

@testset "Laplacian" begin
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = meshMoukalled813()
    # Adjust properties to match fudged values in Moukalled
    cCenters = [ [1.75, 2.0, 0.5], [4.25, 3.5, 0.5], [0,0,0], [0,0,0], [0,0,0], [0,0,0] ]
    cells = [ [1, 2, 3, 4, 5 ], [1], [2], [3], [4], [5] ]
    cVols = [ 8.625, 1, 1, 1, 1, 1 ]
    faces = [ [1,2] ]
    fCenters = fCenters[1:5]
    mesh = [ cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces ]

    gradVals = [ [20.63111, 17.31454, 0.0], [112.625, 133.4375, 0.0], [0,0,0], [0,0,0], [0,0,0] ]
    faceGrad = linInterp(mesh, gradVals)[1]

    Sf = [ 2.5, 1.0, 0.0 ]
    n = normalize(Sf)
    CF = (cCenters[2] - cCenters[1])
    e = normalize(CF)

    #### Test Calculation of Ef ####
    Ef2 = laplacian_Ef("None", Sf, n, CF, e)
    Ef = Sf
    @test Ef == Ef2
    Ef2 = laplacian_Ef("MinCorr", Sf, n, CF, e)
    Ef = [ 2.279, 1.368, 0.0 ]
    @test almostEqual(Ef2, Ef, 3)
    Ef2 = laplacian_Ef("OrthogCorr", Sf, n, CF, e)
    Ef = [ 2.309, 1.385, 0.0 ]
    @test almostEqual(Ef2, Ef, 3)
    Ef2 = laplacian_Ef("OverRelax", Sf, n, CF, e)
    Ef = [ 2.339, 1.403, 0.0 ]
    @test almostEqual(Ef2, Ef, 3)


    #### Test Calculation of Laplacian ####
    # Note that values here are the negative of those presented in Moukalled, since we are interested in only the Laplacian,
    # not the diffusive flux, which is related to the negative Laplacian
    dPhi = 251.578125 - 19.3125
    orthoGrad = dPhi / mag(CF)

    faceIntegral = 0.9122*dPhi - 13.014
    faceIntegral2 = laplacian_FaceFlux(mesh, "MinCorr", 1, faceGrad[1], orthoGrad)
    @test almostEqual(faceIntegral, faceIntegral2, 1)
    faceIntegral = 0.924*dPhi - 16.294
    faceIntegral2 = laplacian_FaceFlux(mesh, "OrthogCorr", 1, faceGrad[1], orthoGrad)
    @test almostEqual(faceIntegral, faceIntegral2, 1)
    faceIntegral = 0.936*dPhi - 19.649
    faceIntegral2 = laplacian_FaceFlux(mesh, "OverRelax", 1, faceGrad[1], orthoGrad)
    @test almostEqual(faceIntegral, faceIntegral2, 1)

    mesh, P, T, U = initializeShockTubeFVM(4)
    Vals = [ 1, 2, 4, 8 ]
    faceGrads = [ 4, 8, 16, 0, 0 ]
    lapl = [ 0.0, 16.0, 32.0, 0.0 ]
    lapl2 = laplacian(mesh, "None", Vals)[1]
    for i in 1:4
        @test almostEqual(lapl[i], lapl2[i])
    end
end;

@testset "Linear Interpolation" begin
    mesh, P, T, U = initializeShockTubeFVM(4)
    cellVals = [ 1, 2, 3, 4 ]
    faceVals = [ 1.5, 2.5, 3.5, 0, 0 ]
    faceVals2 = linInterp(mesh, cellVals)[1]
    @test almostEqual(faceVals, faceVals2)

    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = meshMoukalled813()
    # Adjust properties to match fudged values in Moukalled
    cCenters = [ [1.75, 2.0, 0.5], [4.25, 3.5, 0.5], [0,0,0], [0,0,0], [0,0,0], [0,0,0] ]
    cells = [ [1, 2, 3, 4, 5 ], [1], [2], [3], [4], [5] ]
    cVols = [ 8.625, 1, 1, 1, 1, 1 ]
    faces = [ [1,2] ]
    fCenters = fCenters[1:5]
    mesh = [ cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces ]

    gradVals = [ [20.63111, 17.31454, 0.0], [112.625, 133.4375, 0.0], [0,0,0], [0,0,0], [0,0,0] ]
    faceGrad = linInterp(mesh, gradVals)[1]
    @test almostEqual(faceGrad[1], [ 66.628055, 75.37602, 0.0 ])
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

@testset "Central + AD FVM" begin
    nCells = 4
    P = [ 1, 0.99855553, 0.087046989, 0.1 ]
    U = [ 0, 0.090015665, 0.720126353, 0 ]
    T = [ 0.00348432, 0.00347929, 0.0024263, 0.00278746 ]

    P2, U2, T2, rho2 = central_UnstructuredADFVM(initializeShockTubeFVM(nCells, Pratio=0.1)..., initDt=0.051, endTime=0.05, Cx=0)
    xVel = []
    for i in 1:4
        push!(xVel, U2[i][1])
    end

    @test almostEqual(P, P2, 6)
    @test almostEqual(U, xVel, 6)
    @test almostEqual(T, T2, 6)
end;
