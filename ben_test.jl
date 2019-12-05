using Test
include("finiteVolumeRoe.jl")
#include("finiteDifference.jl")
include("shockTube.jl")
include("vectorFunctions.jl")
include("test.jl")
include("testMeshes.jl")

@testset "Vector Utilities" begin
    v1 = [1, 2, 3]
    @test mag(v1) == 3.7416573867739413

    v2 = [2, 3, 4]
    @test dot(v1, v2) == 20.0

    @test almostEqual(normalize(v1), [0.267261242, 0.534522484, 0.801783726], 9)

    v1 = [ 1.0, 2.0, 3.0 ]
    v2 = [ 2.0, 2.0, 2.0 ]
    crossProduct = [ -2.0, 4.0, -2.0 ]
    crossProduct2 = crossProd(v1, v2)
    @test almostEqual(crossProduct, crossProduct2)
end;

@testset "Shock Tube Initialization" begin
    nCells = 4

    # FVM
    U = [ [0.0,0.0,0.0],  [0.0,0.0,0.0],  [0.0,0.0,0.0],  [0.0,0.0,0.0] ]
    P = [ 1, 1, 0.1, 0.1 ]
    T = [ 0.00348432, 0.00348432, 0.00278746, 0.00278746 ]
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

    rho = 0.999825974
    xMom = 0
    eV2 = 2.501132206
    U = 0.0
    T = 0.00348432
    P = 1

    res1 = [rho, xMom, eV2]
    res2 = encodePrimitives(P, T, U, 287.05, 1005)
    for i in 1:3
        @test almostEqual(res1[i], res2[i], 1e-8)
    end
end;

#=

@testset "Green-Gauss Gradient" begin
    mesh, P, T, U = initializeShockTubeFVM(4)
    faceVals = [ 1, 2, 3, 4, 5 ]
    gradient = [ [0.0,0.0,0.0], [4.0,0.0,0.0], [4.0,0.0,0.0], [0.0,0.0,0.0] ]
    grad2 = greenGaussGrad(mesh, true, faceVals)[1]
    for i in 1:4
        @test almostEqual(gradient[i], grad2[i])
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
    @test almostEqual(grad[1], [ 20.63111, 17.31454, 0.0 ], 5)
end;
=# #Green Gauss


#TODO: Add tests with velo not = 0
@testset "Limiters" begin
    nCells = 4
    Q1 = [1.0, 0.8, 0.6] #Linear
    Q2 = [1.0, 0.8, 0.8] #No grad


    S1_vL = 0.8888888888888888
    S2_L_vL = 0.00
    S2_R_vL = 1.6

    #Test van Leer Limiter
    s1 = vanLeer(Q1[1], Q1[2], Q1[3], dx=0.05)
    s2_l = vanLeer(Q2[1], Q2[2], Q2[3], dx=0.05)
    s2_r = vanLeer(Q2[1], Q2[2], Q2[3], dx=0.05, direc=1)

    @test almostEqual(S1_vL, s1, 1e-6)
    @test almostEqual(S2_L_vL, s2_l, 1e-6)
    @test almostEqual(S2_R_vL, s2_r, 1e-6)

    #Test van Albeda when I figure out what delta to use


end

@testset "MUSCL Interp" begin

    #mesh, P, T, U = initializeShockTubeFVM(4)
    #rho, xMom, eV2 = encodePrimitives(P, T, U, 287.05, 1005)

    rho1 = [1.0, 1.0, 0.5, 0.5]
    rho2 = [1.0, 0.9, 0.6, 0.6]

    QL_1 = 1.0
    QR_1 = 0.5

    QL_2 = 0.863265306
    QR_2 = 0.6

    #Interp to left side
    @test almostEqual(leftFaceVals(rho1[3],rho1[2],rho1[1], 0.05), QL_1)
    @test almostEqual(leftFaceVals(rho2[3], rho2[2], rho2[1], 0.05), QL_2, 1e-6)
    #Interp to right side
    @test almostEqual(rightFaceVals(rho1[4], rho1[3], rho1[2], 0.05), QR_1)
    @test almostEqual(rightFaceVals(rho2[4], rho2[3], rho2[2], 0.05), QR_2)

end

@testset "MUSCL Difference" begin

    mesh, P, T, U = initializeShockTubeFVM(4)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells = size(cells, 1)

    rho = Array{Float64, 1}(undef, nCells)
    xMom = Array{Float64, 1}(undef, nCells)
    eV2 = Array{Float64, 1}(undef, nCells)

    #for i in 1:nCells
    #    rho[i], xMom[i], eV2[i] = encodePrimitives(P[i], T[i], U[i], 287.05, 1005)
    #end

    rho = [0.999825974, 0.999825974, 0.124978067, 0.124978067]
    xMom = [0.0, 0.0, 0.0, 0.0]
    eV2  = [2.501132206, 2.501132206, 0.25011322, 0.25011322]


    #println(rho, xMom, eV2)

    CONS_L3 = [0.999825974, 0.0, 2.501132206]
    CONS_R3 = [0.124978067, 0.0, 0.25011322]

    PRIMS_L3 = [0.0, 1.000452882, 3.502629974]
    PRIMS_R3 = [0.0, 0.100045288, 16.81336771]


    cl, pl, cr, pr = musclDifference(mesh, rho, xMom, eV2, 0.05)

    #print(cl[:,3], typeof(cl))

    @test almostEqual(cl[2,:], CONS_L3, 1e-6)
    @test almostEqual(cr[2,:], CONS_R3)
    @test almostEqual(pl[2,:], PRIMS_L3, 1e-6)
    @test almostEqual(pr[2,:], PRIMS_R3, 1e-6)

end

@testset "Roe Averaging" begin

    res = [0.353491609, 0.0, 6.979444402, 1.670861382]

    rho_l = 0.999825974
    rho_r = 0.124978067
    u_l = 0.0
    u_r = 0.0
    H_l = 3.502629974
    H_r = 16.81336771

    n = 1


    res2 = roeAveraged(n, rho_l, rho_r, u_l, u_r, H_l, H_r)

    @test almostEqual(res, res2, 1e-6)

    #Add a test with velo??
end

@testset "Flux Vectors" begin
    #Hand calcs results
    FL = [0.0 1.000452882 0.0]
    FR = [0.0 0.100045288 0.0]
    F1 = [0.0 0.0 0.0]
    F2 = [-0.269427976 -0.450203797 -1.880683125]
    F3 = [-0.269427976 0.450203797 -1.880683125]

    F_tot = [0.269427976, 0.550249085, 1.880683125]

    #Declare variables for inputs
    nFaces = Array{Float64, 1}(undef, 1)

    CONS_L3 = Array{Float64,2}(undef, 1, 3)
    CONS_R3 = Array{Float64,2}(undef, 1, 3)

    PRIMS_L3 = Array{Float64,2}(undef, 1, 3)
    PRIMS_R3 = Array{Float64,2}(undef, 1, 3)

    CONS_L3 = [[0.999825974 0.0 2.501132206]; [0.0 0.0 0.0]]
    CONS_R3 = [[0.124978067 0.0 0.25011322]; [0.0 0.0 0.0]]

    #println(CONS_L3[1,3])
    #println(typeof(CONS_R3))

    PRIMS_L3 = [[0.0 1.000452882 3.502629974]; [0.0 0.0 0.0]]
    PRIMS_R3 = [[0.0 0.100045288 16.81336771]; [0.0 0.0 0.0]]

    f_l = fluxVector(CONS_L3, PRIMS_L3, nFaces)
    f_r = fluxVector(CONS_R3, PRIMS_R3, nFaces)

    @test almostEqual(f_l, FL, 1e-6)
    @test almostEqual(f_r, FR, 1e-6)

    rhoR = 0.353491609
    uR = 0.0
    HR = 6.980281542
    aR = 1.670961585
    deltaRho = -0.874847907
    deltaU = 0.0
    deltaP = -0.900407594
    deltaA = 1.409475552

    f_1, f_2, f_3 = eigenFluxVectors(rhoR, uR, HR, aR, deltaRho, deltaU, deltaP, deltaA, nFaces)


    res = findFluxes(f_l, f_r, f_1, f_2, f_3, nFaces)


    @test almostEqual(f_1, F1, 1e-6)
    @test almostEqual(f_2, F2, 1e-6)
    @test almostEqual(f_3, F3, 1e-6)

    for i in 1:3
        @test almostEqual(res[i], F_tot[i], 1e-6)
    end
    #println("Res: ", res, "\n")
    #println("F_tot: ", F_tot)

end



@testset "Roe FVM" begin
    nCells = 4
    P = [ 1, 0.848275045, 0.241434252, 0.1 ]
    U = [ 0, 0.095689296, 0.503421285, 0 ]
    T = [ 0.00348432, 0.003140533, 0.004702554, 0.00278746 ]

    #P2, U2, T2, rho2 = central_UnstructuredADFVM(initializeShockTubeFVM(nCells, Pratio=0.1)..., initDt=0.051, endTime=0.05, Cx=0)
    P2, U2, T2, rho2 = upwindFVMRoe1D(initializeShockTubeFVM(nCells)..., initDt=0.051, endTime=0.05)

    xVel = []
    for i in 1:4
        push!(xVel, U2[i][1])
    end

    @test almostEqual(P, P2, 1e-6)
    @test almostEqual(U, xVel, 1e-6)
    @test almostEqual(T, T2, 1e-6)
end;
