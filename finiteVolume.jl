using Printf
include("constitutiveRelations.jl")
include("vectorFunctions.jl")
include("timeDiscretizations.jl")
include("mesh.jl")
include("output.jl")
include("dataStructures.jl")

__precompile__()

######################### Initialization ###########################
# Returns cellPrimitives matrix for uniform solution
function initializeUniformSolution3D(mesh, P, T, Ux, Uy, Uz)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

    cellPrimitives = zeros(nCells, 5)
    for c in 1:nCells
        cellPrimitives[c, :] = [ P, T, Ux, Uy, Uz ]
    end

    return cellPrimitives
end

######################### CFL #######################################
# function maxCFL3D(mesh::Mesh, sln::SolutionState, dt, gamma=1.4, R=287.05)
#     nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
#
#     maxCFL = 0.0
#     for c in 1:nCells
#         a = sqrt(gamma * R * sln.cellPrimitives[c,2])
#
#         CFL = 0.0
#         for d in 1:3
#             cDeltaT = (abs(sln.cellPrimitives[c,d+2]) + a)*dt
#             CFL += cDeltaT / mesh.cellSizes[c, d]
#         end
#         maxCFL = max(maxCFL, CFL)
#     end
#
#     return maxCFL
# end

# Calculates CFL at each cell. Expects sln.cellState, sln.cellPrimitives and sln.faceFluxes to be up to date
function CFL!(CFL, mesh::Mesh, sln::SolutionState, dt=1, gamma=1.4, R=287.05)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

    fill!(CFL, 0.0)
    #faceRhoT = linInterp_3D(mesh, hcat(sln.cellState[:,1], sln.cellPrimitives[:,2]))
    sln.faceState[1:nFaces-nBdryFaces,1] = linInterpToFace_3D(mesh, sln, 1, 1)
    sln.facePrimitives[1:nFaces-nBdryFaces,2] = linInterpToFace_3D(mesh, sln, 3, 2)

    for f in 1:nFaces
        ownerCell = mesh.faces[f][1]
        neighbourCell = mesh.faces[f][2]

        #faceRho = faceRhoT[f, 1]
        #faceT = faceRhoT[f, 2]
        faceRho = sln.faceState[f,1]
        faceT = sln.facePrimitives[f,2]


        # Use cell center values on the boundary
        #if neighbourCell == -1
        #    faceRho = sln.cellState[ownerCell, 1]
        #    faceT = sln.cellPrimitives[ownerCell, 2]
        #end

        @views faceVel = sln.faceFluxes[f, 1:3] ./ faceRho
        @views flux = abs(dot(faceVel, mesh.fAVecs[f]))*dt

        if faceT < 0.0 #Empty boundary will always have a T value of 0
            println("Found it!")
        end

        a = sqrt(gamma * R * faceT)
        flux += mag(mesh.fAVecs[f])*a*dt

        if faceT != 0
            CFL[ownerCell] += flux
        end


        if neighbourCell > -1
            CFL[neighbourCell] += flux
        end
    end

    CFL ./= (2 .* mesh.cVols)
end

######################### Gradient Computation #######################
#TODO: make all these function return a single array if you pass in a single value
#TODO: Add variable stencils for WLS, instead of just face neighbours
#TODO: Document
function leastSqGrad(mesh::Mesh, sln::SolutionState, slnPointer, slnColPointer=-1)
    # slnPointer MUST be 1:3 for boundaries to work
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    slnArray = [sln.cellState, sln.cellFluxes, sln.cellPrimitives, sln.fluxResiduals, sln.faceState, sln.faceFluxes, sln.facePrimitives]

    matrix = slnArray[slnPointer]
    faceMatrix = slnArray[slnPointer+4]

    if slnColPointer != -1
        matrix = matrix[:,slnColPointer]
        faceMatrix = faceMatrix[:,slnColPointer]
    end

    nVars = size(matrix, 2)

    grad = zeros(nCells, nVars, 3)

    L11 = zeros(nCells,nVars)
    L12 = zeros(nCells,nVars)
    L13 = zeros(nCells,nVars)
    L22 = zeros(nCells,nVars)
    L23 = zeros(nCells,nVars)
    L33 = zeros(nCells,nVars)
    L1f = zeros(nCells,nVars)
    L2f = zeros(nCells,nVars)
    L3f = zeros(nCells,nVars)

    @fastmath for f in 1:nFaces-nBdryFaces

        c1 = mesh.faces[f][1]
        c2 = mesh.faces[f][2]

        dx = mesh.cCenters[c2][1] - mesh.cCenters[c1][1]
        dy = mesh.cCenters[c2][2] - mesh.cCenters[c1][2]
        dz = mesh.cCenters[c2][3] - mesh.cCenters[c1][3]

        weight = 1 / sqrt(dx*dx + dy*dy + dz*dz)

        wdx = weight * dx
        wdy = weight * dy
        wdz = weight * dz

        for v in 1:nVars
            dv = matrix[c2,v] - matrix[c1,v]
            wdv = weight * dv

            computeLeastSqMat!(L11,L12,L13,L22,L23,L33,L1f,L2f,L3f, wdx, wdy, wdz, wdv, c1, c2, v)
        end
    end


    # Deal with boundary faces, and add them to the matrix vectors
    @fastmath for f in nFaces-nBdryFaces+1:nFaces
        c1 = mesh.faces[f][1]
        c2 = -1

        # Use a pseudo "mirror" cell for each boundary cell
        dx = (mesh.fCenters[f][1] - mesh.cCenters[c1][1]) * 2
        dy = (mesh.fCenters[f][2] - mesh.cCenters[c1][2]) * 2
        dz = (mesh.fCenters[f][3] - mesh.cCenters[c1][3]) * 2

        weight = 1 / sqrt(dx^2 + dy^2 + dz^2)

        wdx = weight * dx
        wdy = weight * dy
        wdz = weight * dz


        for v in 1:nVars
            #For all variables, just linearly interpolate what the value  would be at the pseudocell
            #I THINK that should work for everything
            slope = faceMatrix[f,v] - matrix[c1,v]
            pseudoVal = matrix[c1,v] + 2*slope # No distances b/c pseudocell is mirrored, so face is always exactly halfway

            dv = pseudoVal - matrix[c1,v]
            wdv = weight * dv

            computeLeastSqMat!(L11,L12,L13,L22,L23,L33,L1f,L2f,L3f, wdx, wdy, wdz, wdv, c1, c2, v)
        end
    end

    L = zeros(3,3)
    R = zeros(3,1)

    for c in 1:nCells
        for v in 1:nVars
            fillLeastSqMat!(L, R, L11[c,v], L12[c,v], L13[c,v], L22[c,v], L23[c,v], L33[c,v], L1f[c,v], L2f[c,v], L3f[c,v])
            x = inv(L) * R

            grad[c,v, :] = x
        end

    end

    return grad
end

# Updates the elements in the least squares matrix
function computeLeastSqMat!(L11, L12, L13, L22, L23, L33, L1f, L2f, L3f, wdx, wdy, wdz, wdv, c1, c2, v)
    if c2 != -1
        L11[c1,v] += wdx^2
        L12[c1,v] += wdx * wdy
        L13[c1,v] += wdx * wdz
        L22[c1,v] += wdy^2
        L23[c1,v] += wdy * wdz
        L33[c1,v] += wdz^2

        L1f[c1,v] += wdx * wdv
        L2f[c1,v] += wdy * wdv
        L3f[c1,v] += wdz * wdv

        L11[c2,v] += wdx^2
        L12[c2,v] += wdx * wdy
        L13[c2,v] += wdx * wdz
        L22[c2,v] += wdy^2
        L23[c2,v] += wdy * wdz
        L33[c2,v] += wdz^2

        L1f[c2,v] += wdx * wdv
        L2f[c2,v] += wdy * wdv
        L3f[c2,v] += wdz * wdv
    else
        L11[c1,v] += wdx^2
        L12[c1,v] += wdx * wdy
        L13[c1,v] += wdx * wdz
        L22[c1,v] += wdy^2
        L23[c1,v] += wdy * wdz
        L33[c1,v] += wdz^2

        L1f[c1,v] += wdx * wdv
        L2f[c1,v] += wdy * wdv
        L3f[c1,v] += wdz * wdv
    end
end

function fillLeastSqMat!(L, R, L11, L12, L13, L22, L23, L33, R1, R2, R3)
    L[1,1] = L11
    L[1,2] = L12
    L[1,3] = L13
    L[2,1] = L12
    L[2,2] = L22
    L[2,3] = L23
    L[3,1] = L13
    L[3,2] = L23
    L[3,3] = L33

    R[1,1] = R1
    R[2,1] = R2
    R[3,1] = R3
end



#=
    Takes the gradient of (scalar) data provided in matrix form (passed into arugment 'matrix'):
    Cell      x1      x2      x3
    Cell 1    x1_1    x2_1    x3_1
    Cell 2    x1_2    x2_2    x3_2
    ...
    (where x1, x2, and x3 are arbitrary scalars)
    and output a THREE-DIMENSIONAL gradient arrays of the following form
    Cell      x1          x2          x3
    Cell 1    grad(x1)_1  grad(x2)_1  grad(x3)_1
    Cell 2    grad(x1)_2  grad(x2)_2  grad(x3)_2
    ...
    Where each grad(xa)_b is made up of THREE elements for the (x,y,z) directions

    Ex. Gradient @ cell 1 of P would be: greenGaussGrad(mesh, P)[1, 1, :]
=#
function greenGaussGrad(mesh::Mesh, matrix::AbstractArray{Float64, 2}, valuesAtFaces=false)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    nVars = size(matrix, 2)

    #Interpolate values to faces if necessary
    if valuesAtFaces != true
        #faceVals = linInterp_3D(mesh, matrix)
        println("Interpolating inside GG grad is depracated! Please only input matrix of face values!")
        #faceVals = linInterpToFace_3D(mesh, sln, 3, 1)
    else
        faceVals = matrix
    end

    # Create matrix to hold gradients
    grad = zeros(nCells, nVars, 3)

    # Integrate fluxes from each face
    @inbounds @fastmath for f in eachindex(mesh.faces)
        ownerCell = mesh.faces[f][1]
        neighbourCell = mesh.faces[f][2]

        for v in 1:nVars
            for d in 1:3
                # Every face has an owner
                grad[ownerCell, v, d] += mesh.fAVecs[f][d] * faceVals[f, v]

                # Boundary faces don't - could split into two loops
                if neighbourCell > -1
                    grad[neighbourCell, v, d] -= mesh.fAVecs[f][d] * faceVals[f, v]
                end
            end
        end
    end

    # Divide integral by cell volume to obtain gradients
    @inbounds @fastmath for c in 1:nCells
        for v in 1:nVars
            for d in 1:3
                grad[c,v,d] /= mesh.cVols[c]
            end
        end
    end

    return grad
end


############ Limiters ###################

function vanAlbeda(q_ip, q_i, q_im, dx; direc=0)
    #=Impplements the vanAlbeda flux limiter, with a delta value of dx^3
    =#
    delta = dx.*dx.*dx
    #delta = dx
    #delta = 1.0
    if direc ==0
        r = (q_i - q_im)/(q_ip - q_i + delta)
    else
        r = (q_ip - q_i)/(q_i - q_im + delta)
        #r = (q_i - q_im)/(q_ip - q_i + delta)
    end
    if r > 0
        #s = 2*r / (r^2 + 1)
        s = (r^2 + r)/(r^2 + 1)
    else
        s = 0
    end
    #s = (2 * (q_ip-q_i)*(q_i-q_im) + delta )/( (q_ip-q_i)^2 * (q_i-q_im)^2 + delta )
    if s > 0.05
        s = 0.05
    end
    return s
end



####################### Face value interpolation ######################

function MUSCL_difference(mesh::Mesh, sln::SolutionState; sigma=-1)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    println("THeres $nCells cells, $nFaces faces, $nBdryFaces bdryfaces")

    nVars = size(sln.cellState, 2)

    interpOwner = zeros(nFaces-nBdryFaces, nVars)
    interpNeighbour = zeros(nFaces-nBdryFaces, nVars)

    sln.faceState[1:nFaces-nBdryFaces,:] = linInterpToFace_3D(mesh, sln, 1)

    for v in 1:nVars
        #@views P = sln.cellPrimitives[:,1]
        #P = reshape(P, nCells, :)


        facesCons = zeros(nFaces,1)
        facesCons[:,:] = sln.faceState[:,v]

        gradCons = greenGaussGrad(mesh, facesCons, true)

        for f in 1:nFaces-nBdryFaces
            ownerCell = mesh.faces[f][1]
            neighbourCell = mesh.faces[f][2]
            d = mesh.cCenters[neighbourCell] .- mesh.cCenters[ownerCell]

            # 1. Find variable at owner/neighbour cells
            qOwner = sln.cellState[ownerCell, v]
            qNeighbour = sln.cellState[neighbourCell, v]

            # 2. Calculate conserved quants at 'virtual' far-owner and far-neighbour cells using the pressure gradient (2nd-order)
            @views farOwnerQ = qOwner - dot(d, gradCons[ownerCell, 1, :])
            @views farNeighbourQ = qNeighbour + dot(d, gradCons[neighbourCell, 1, :])

            dx = 0.1
            sLeft = vanAlbeda(qNeighbour, qOwner, farOwnerQ, dx)
            sRight = vanAlbeda(qOwner, qNeighbour, farNeighbourQ, dx)

            # if f > 153 && f < 155 && v == 2
            #     println("This is face $f, var $v, and limiters are: $sLeft, $sRight")
            #     println("Owners is: $qOwner, $farOwnerQ")
            #     println("Neighbours is: $qNeighbour, $farNeighbourQ")
            # end

            #TODO: Fix MUSCL differencing!
            interpOwner[f,v] = qOwner + 0.5*dot(d,gradCons[ownerCell,1,:])
            interpNeighbour[f,v] = qNeighbour - 0.5*dot(d, gradCons[neighbourCell,1,:])
            # interpOwner[f,v] = qOwner + sLeft*0.25*( (1-sigma*sLeft)*(qOwner-farOwnerQ) + (1+sigma*sLeft)*(qNeighbour-qOwner) ) * qOwner
            # interpNeighbour[f,v] = qNeighbour - sRight*0.25*( (1-sigma*sRight)*(farNeighbourQ-qNeighbour) + (1+sigma*sRight)*(qNeighbour-qOwner) )*qNeighbour
        end
    end
    return interpOwner, interpNeighbour
end

#=
    Interpolates to all INTERIOR faces
    Arbitrary value matrix interpolation
    Example Input Matrix:
    Cell      x1       x2       x3
    Cell 1    x1_c1    x2_c1    x3_c1
    Cell 2    x1_c2    x2_c2    x3_c2
    ...
    Outputs a matrix of the following form
    Cell      x1       x2       x3
    Face 1    x1_f1    x2_f1    x3_f1
    Face 2    x1_f2    x2_f2    x3_f2
    ...
=#
# Rewrite so that it accepts the sln structure, and then writes results into sln.faceState or sln.faceFluxes (or whatever required??)
#Always writes to sln face. No longer mutates inputs.
function linInterpToFace_3D(mesh::Mesh, sln::SolutionState, slnPointer, slnPointerCols=-1)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    slnArray = [sln.cellState, sln.cellFluxes, sln.cellPrimitives, sln.fluxResiduals, sln.faceState, sln.faceFluxes, sln.facePrimitives]

    matrix = slnArray[slnPointer]

    if slnPointerCols != -1
        matrix = matrix[:,slnPointerCols]
    end

    nVars = size(matrix, 2)

    # Make a new matrix if one is not passed in
    #targetSize = size(sln_array[slnPointer],1)

    faceVals = zeros(nFaces-nBdryFaces, nVars)

    # Boundary face fluxes must be set separately
    @inbounds @fastmath for f in 1:nFaces-nBdryFaces
        # Find value at face using linear interpolation
        c1 = mesh.faces[f][1]
        c2 = mesh.faces[f][2]

        #TODO: Precompute these distances?
        c1Dist = mag(mesh.cCenters[c1] .- mesh.fCenters[f])
        c2Dist = mag(mesh.cCenters[c2] .- mesh.fCenters[f])
        totalDist = c1Dist + c2Dist

        for v in 1:nVars
            faceVals[f, v] = matrix[c1, v]*(c2Dist/totalDist) + matrix[c2, v]*(c1Dist/totalDist)
        end
    end

    return faceVals
end
#=
function linInterp_3D(mesh::Mesh, matrix::AbstractArray{Float64, 2}, faceVals=zeros(2,2))
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    nVars = size(matrix, 2)

    # Make a new matrix if one is not passed in
    if faceVals == zeros(2,2)
        faceVals = zeros(nFaces, nVars)
    end

    # Boundary face fluxes must be set separately
    @inbounds @fastmath for f in 1:nFaces-nBdryFaces
        # Find value at face using linear interpolation
        c1 = mesh.faces[f][1]
        c2 = mesh.faces[f][2]

        #TODO: Precompute these distances?
        c1Dist = mag(mesh.cCenters[c1] .- mesh.fCenters[f])
        c2Dist = mag(mesh.cCenters[c2] .- mesh.fCenters[f])
        totalDist = c1Dist + c2Dist

        for v in 1:nVars
            faceVals[f, v] = matrix[c1, v]*(c2Dist/totalDist) + matrix[c2, v]*(c1Dist/totalDist)
        end
    end

    return faceVals
end
=#

# Similar to linInterp_3D. Instead of linearly interpolating, selects the maximum value of the two adjacent cell centers as the face value
function maxInterp(mesh::Mesh, vars...)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

    nVars = size(vars, 1)
    faceVals = zeros(nFaces, nVars)

    # Boundary face fluxes must be set separately
    @inbounds @fastmath for f in 1:nFaces-nBdryFaces
        c1 = mesh.faces[f][1]
        c2 = mesh.faces[f][2]

        for v in 1:nVars
            faceVals[f, v] = max(vars[v][c1], vars[v][c2])
        end
    end

    return faceVals
end

# Similar to linInterp_3D. Instead of linearly interpolating, calculates the change (delta) of each variable across each face (required for JST method)
function faceDeltas(mesh::Mesh, sln::SolutionState)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    nVars = size(sln.cellState, 2)

    faceDeltas = zeros(nFaces, nVars)

    # Boundary face fluxes must be set separately (faceDelta is zero at all possible boundary conditions right now)
    @inbounds @fastmath for f in 1:nFaces-nBdryFaces
        ownerCell = mesh.faces[f][1]
        neighbourCell = mesh.faces[f][2]

        for v in 1:nVars
            faceDeltas[f, v] = sln.cellState[neighbourCell, v] - sln.cellState[ownerCell, v]
        end
    end

    return faceDeltas
end

# TODO: MUSCL Interp

######### Flux vectors

function findFluxVectors(mesh::Mesh, sln::SolutionState, roeAveraged, eig1, eig2, eig3, ownerPrims, neighbourPrims, ownerCons, neighbourCons; gamma=1.4)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    #nVars = size(sln.cellFluxes, 1)
    nVars = 5

    nIntFaces = nFaces-nBdryFaces

    FL = zeros(nIntFaces, nVars)
    FR = zeros(nIntFaces, nVars)

    F1 = zeros(nIntFaces, nVars)
    F2 = zeros(nIntFaces, nVars)
    F3 = zeros(nIntFaces, nVars)

    # cons: rho, rhoU, rhoV, rhoW, eV2
    # prims: P, T, U, V, W
    # roe: rho, U, V, W, H

    for f in 1:nFaces-nBdryFaces

        fL_1 = [ownerCons[f,2], ownerCons[f,3], ownerCons[f,4]]
        fL_2 = [ownerCons[f,2]*ownerPrims[f,3]+ownerPrims[f,1], ownerCons[f,2]*ownerPrims[f,4], ownerCons[f,2]*ownerPrims[f,5]]
        fL_3 = [ownerCons[f,3]*ownerPrims[f,3], ownerCons[f,3]*ownerPrims[f,4]+ownerPrims[f,1], ownerCons[f,3]*ownerPrims[f,5]]
        fL_4 = [ownerCons[f,4]*ownerPrims[f,3], ownerCons[f,4]*ownerPrims[f,4], ownerCons[f,4]*ownerPrims[f,5] + ownerPrims[f,1]]
        fL_5 = [ownerPrims[f,3]*( ownerCons[f,5]+ownerPrims[f,1]  ), ownerPrims[f,4]*( ownerCons[f,5]+ownerPrims[f,1]  ), ownerPrims[f,5]*( ownerCons[f,5]+ownerPrims[f,1]  )  ]
        fL = [fL_1, fL_2, fL_3, fL_4, fL_5]
        normFVec = mesh.fAVecs[f] / mag(mesh.fAVecs[f])


        # FL[f,1] = dot(normFVec, fL_1)
        # FL[f,2] = dot(normFVec, fL_2)
        # FL[f,3] = dot(normFVec, fL_3)
        # FL[f,4] = dot(normFVec, fL_4)
        # FL[f,5] = dot(normFVec, fL_5)

        fR_1 = [neighbourCons[f,2], neighbourCons[f,3], neighbourCons[f,4]]
        fR_2 = [neighbourCons[f,2]*neighbourPrims[f,3]+neighbourPrims[f,1], neighbourCons[f,2]*neighbourPrims[f,4], neighbourCons[f,2]*neighbourPrims[f,5]]
        fR_3 = [neighbourCons[f,3]*neighbourPrims[f,3], neighbourCons[f,3]*neighbourPrims[f,4]+neighbourPrims[f,1], neighbourCons[f,3]*neighbourPrims[f,5]]
        fR_4 = [neighbourCons[f,4]*neighbourPrims[f,3], neighbourCons[f,4]*neighbourPrims[f,4], neighbourCons[f,4]*neighbourPrims[f,5] + neighbourPrims[f,1]]
        fR_5 = [neighbourPrims[f,3]*( neighbourCons[f,5]+neighbourPrims[f,1]  ), neighbourPrims[f,4]*( neighbourCons[f,5]+neighbourPrims[f,1]  ), neighbourPrims[f,5]*( neighbourCons[f,5]+neighbourPrims[f,1]  )  ]
        fR = [fR_1, fR_2, fR_3, fR_4, fR_5]
        # FR[f,1] = dot(normFVec, fR_1)
        # FR[f,2] = dot(normFVec, fR_2)
        # FR[f,3] = dot(normFVec, fR_3)
        # FR[f,4] = dot(normFVec, fR_4)
        # FR[f,5] = dot(normFVec, fR_5)

        # println("Owner vecs are:")
        # display(ownerCons[f,:])
        # display(ownerPrims[f,:])
        # println("And neighbour vecs are:")
        # display(neighbourCons[f,:])
        # display(neighbourPrims[f,:])
        # println("And flux vecs")
        #
        # display(normFVec)
        # display(normFVec[2])


        sound = sqrt( (gamma-1)*(roeAveraged[f,5] - 0.5*( mag(roeAveraged[f,2:4])^2 ) ) )
        UBar = (neighbourPrims[f,3]-ownerPrims[f,3])*normFVec[1] + (neighbourPrims[f,4]-ownerPrims[f,4])*normFVec[2] + (neighbourPrims[f,5]-ownerPrims[f,5])*normFVec[3]
        URoe = dot(normFVec, roeAveraged[f,2:4])

        f_11 = (neighbourCons[f,1]-ownerCons[f,1])-(neighbourPrims[f,1]-ownerPrims[f,1])/(sound^2) #TODO: Some uncertainty about the sign +/- here
        f_12 = roeAveraged[f,2]*( f_11 ) + roeAveraged[f,1]*( (neighbourPrims[f,3]-ownerPrims[f,3])-normFVec[1]*UBar )
        f_13 = roeAveraged[f,3]*( f_11 ) + roeAveraged[f,1]*( (neighbourPrims[f,4]-ownerPrims[f,4])-normFVec[2]*UBar )
        f_14 = roeAveraged[f,4]*( f_11 ) + roeAveraged[f,1]*( (neighbourPrims[f,5]-ownerPrims[f,5])-normFVec[3]*UBar )
        f_15 = 0.5 * (mag(roeAveraged[f,2:4])^2) * (f_11) + roeAveraged[f,1]*( roeAveraged[f,2]*(neighbourPrims[f,3]-ownerPrims[f,3]) + roeAveraged[f,3]*(neighbourPrims[f,4]-ownerPrims[f,4]) + roeAveraged[f,4]*(neighbourPrims[f,5]-ownerPrims[f,5]) - URoe * UBar )
        F_1 = [f_11, f_12, f_13, f_14, f_15]

        # F1[f,1] = abs(eig1[f]) * f_11
        # F1[f,2] = abs(eig1[f]) * f_12
        # F1[f,3] = abs(eig1[f]) * f_13
        # F1[f,4] = abs(eig1[f]) * f_14
        # F1[f,5] = abs(eig1[f]) * f_15


        f2Pre = (neighbourPrims[f,1]-ownerPrims[f,1])/(2*sound^2) + (roeAveraged[f,1]*UBar)/(2*sound) #TODO: Another place with +/- uncertainty

        f_21 = 1
        f_22 = roeAveraged[f,2] + normFVec[1] * sound
        f_23 = roeAveraged[f,3] + normFVec[2] * sound
        f_24 = roeAveraged[f,4] + normFVec[3] * sound
        f_25 = roeAveraged[f,5] + URoe * sound
        F_2 = [f_21, f_22, f_23, f_24, f_25]

        f3Pre = (neighbourPrims[f,1]-ownerPrims[f,1])/(2*sound^2) - (roeAveraged[f,1]*UBar)/(2*sound)

        f_31 = 1
        f_32 = roeAveraged[f,2] - normFVec[1]*sound
        f_33 = roeAveraged[f,3] - normFVec[2]*sound
        f_34 = roeAveraged[f,4] - normFVec[3]*sound
        f_35 = roeAveraged[f,5] - URoe*sound
        F_3 = [f_31, f_32, f_33, f_34, f_35]


        for i in 1:5
            FL[f,i] = dot(normFVec, fL[i])
            FR[f,i] = dot(normFVec, fR[i])
            F1[f,i] = abs(eig1[f]) * F_1[i]
            F2[f,i] = abs(eig2[f]) * f2Pre * F_2[i]
            F3[f,i] = abs(eig3[f]) * f3Pre * F_3[i]
        end

        # display(fR_1)
        # display(fR_2)
        # display(fR_3)
        # display(fR_4)
        # display(fR_5)
        # display(F3[f,:])

        # println("$breakdown")
    end


    return FL, FR, F1, F2, F3

end


function findEigenvalues(roe, left, right; gamma=1.4, K=1)
    nFaces = size(roe,1)

    sound = zeros(nFaces)

    eig1 = zeros(nFaces)
    eig2 = zeros(nFaces)
    eig3 = zeros(nFaces)

    for f in 1:nFaces
        #display(roe[f,5])
        #println("Face is $f")
        sound = sqrt( (gamma-1)*(roe[f,5] - 0.5*( mag(roe[f,2:4])^2 ) ) )

        eig1[f] = mag( roe[f,2:4] )
        eig2[f] = mag( roe[f,2:4] ) + sound
        eig3[f] = mag( roe[f,2:4] ) - sound

        epsilon = K * max( ( mag(right[f,3:5])-mag(left[f,3:5]) ), 0)

        if epsilon > abs(eig1[f])
            eig1[f] = ( eig1[f]^2 + epsilon^2 )/(2*epsilon)
            println("Replaced eig1")
        end

        if epsilon > abs(eig2[f])
            eig2[f] = ( eig2[f]^2 + epsilon^2 )/(2*epsilon)
            println("Replaced eig2")
        end

        if epsilon > abs(eig3[f])
            println("Replaced eig3 for face $f, sound is $sound !")
            display(eig3[f])
            eig3[f] = ( eig3[f]^2 + epsilon^2 )/(2*epsilon)
            println("Is now")
            display(eig3[f])
        end
    end

    return eig1, eig2, eig3
end




######################### Convective Term Things #######################
# Calculates eps2 and eps4, the second and fourth-order artificial diffusion coefficients used in the JST method
function unstructured_JSTEps(mesh::Mesh, sln::SolutionState, k2=0.5, k4=(1/32), c4=1, gamma=1.4, R=287.05)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

    # Calc Pressure Gradient
    @views P = sln.cellPrimitives[:,1]
    P = reshape(P, nCells, :)
    sln.facePrimitives[1:nFaces-nBdryFaces,1] = linInterpToFace_3D(mesh, sln, 3, 1)

    facesP = zeros(nFaces,1)
    facesP[:,:] = sln.facePrimitives[:,1]

    gradP = greenGaussGrad(mesh, facesP, true)
    gradP = reshape(gradP, nCells, 3)

    sj = zeros(nCells) # 'Sensor' used to detect shock waves and apply second-order artificial diffusion to stabilize solution in their vicinity
    sjCount = zeros(nCells) # Store the number of sj's calculated for each cell, cell-center value will be the average of all of them
    @inbounds @fastmath for f in 1:nFaces-nBdryFaces
        # At each internal face, calculate sj, rj, eps2, eps4
        ownerCell = mesh.faces[f][1]
        neighbourCell = mesh.faces[f][2]
        d = mesh.cCenters[neighbourCell] .- mesh.cCenters[ownerCell]

        # 1. Find pressure at owner/neighbour cells
        oP = P[ownerCell]
        nP = P[neighbourCell]

        # 2. Calculate pressures at 'virtual' far-owner and far-neighbour cells using the pressure gradient (2nd-order)
        @views farOwnerP = nP - 2*dot(d, gradP[ownerCell, :])
        @views farNeighbourP = oP + 2*dot(d, gradP[neighbourCell, :])

        # 3. With the known and virtual values, can calculate sj at each cell center.
        sj[ownerCell] += (abs( nP - 2*oP + farOwnerP )/ max( abs(nP - oP) + abs(oP - farOwnerP), 0.0000000001))^2
        sjCount[ownerCell] += 1
        sj[neighbourCell] += (abs( oP - 2*nP + farNeighbourP )/ max( abs(farNeighbourP - nP) + abs(nP - oP), 0.0000000001))^2
        sjCount[neighbourCell] += 1
    end

    rj = zeros(nCells) # 'Spectral radius' -> maximum possible speed of wave propagation relative to mesh
    @inbounds @fastmath for c in 1:nCells
        @views rj[c] = mag(sln.cellPrimitives[c,3:5]) +  sqrt(gamma * R * sln.cellPrimitives[c,2]) # Velocity magnitude + speed of sound
        sj[c] /= sjCount[c] # Average the sj's computed at each face for each cell
    end

    # Values of rj, sj at faces is the maximum of their values at the two adjacent cell centers
    rjsjF = maxInterp(mesh, rj, sj) # column one is rj, column two is sj, both at face centers

    # Calculate eps2 and eps4
    eps2 = zeros(nFaces)
    eps4 = zeros(nFaces)
    for f in 1:nFaces-nBdryFaces
        eps2[f] = k2 * rjsjF[f,2] * rjsjF[f,1]
        eps4[f] = max(0, k4*rjsjF[f,1] - c4*eps2[f])
    end

    return eps2, eps4
end

#=
    Inputs: Expects that sln.cellState, sln.cellPrimitives and sln.cellFluxes are up-to-date
    Outputs: Updates sln.faceFluxes and sln.cellResiduals
    Returns: Updated sln.cellResiduals

    Applies classical JST method: central differencing + JST artificial diffusion. Each face treated as a 1D problem
    http://aero-comlab.stanford.edu/Papers/jst_2015_updated_07_03_2015.pdf -  see especially pg.5-6
=#
function unstructured_JSTFlux(mesh::Mesh, sln::SolutionState, boundaryConditions, gamma, R)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    nVars = size(sln.cellState, 2)

    #### 1. Apply boundary conditions ####
    for b in 1:nBoundaries
        bFunctionIndex = 2*b-1
        boundaryConditions[bFunctionIndex](mesh, sln, b, boundaryConditions[bFunctionIndex+1])
    end

    #### 2. Centrally differenced fluxes ####
    #linInterp_3D(mesh, sln.cellFluxes, sln.faceFluxes)
    #sln.faceFluxes = linInterp_3D(mesh, sln.cellFluxes, sln.faceFluxes)
    sln.faceFluxes[1:nFaces-nBdryFaces,:] = linInterpToFace_3D(mesh, sln, 2, -1)

    #### 3. Add JST artificial Diffusion ####
    fDeltas = faceDeltas(mesh, sln) # TODO: Figure out how to deal with boundaries, as currently they are not included here
    fDGrads = greenGaussGrad(mesh, fDeltas, true)
    eps2, eps4 = unstructured_JSTEps(mesh, sln, 0.5, (1.2/32), 1, gamma, R)

    diffusionFlux = zeros(nVars)
    @inbounds @fastmath for f in 1:nFaces-nBdryFaces
        ownerCell = mesh.faces[f][1]
        neighbourCell = mesh.faces[f][2]
        d = mesh.cCenters[neighbourCell] .- mesh.cCenters[ownerCell]

        @views fD = fDeltas[f,:]
        @views farOwnerfD = fD .- dot(d, fDGrads[ownerCell,:,:])
        @views farNeighbourfD = fD .+ dot(d, fDGrads[ownerCell,:,:])

        diffusionFlux = eps2[f]*fD - eps4[f]*(farNeighbourfD - 2*fD + farOwnerfD)

        # Add diffusion flux in component form
        unitFA = normalize(mesh.fAVecs[f])
        for v in 1:nVars
            i1 = (v-1)*3+1
            i2 = i1+2
            sln.faceFluxes[f,i1:i2] .-= (diffusionFlux[v] .* unitFA)
        end
    end


    #### 4. Integrate fluxes at in/out of each cell (sln.faceFluxes) to get change in cell center values (sln.fluxResiduals) ####
    return integrateFluxes_unstructured3D(mesh, sln, boundaryConditions)
end


######## ROE Flux ########################################
#=
    Inputs: Expects that sln.cellState, sln.cellPrimitives and sln.cellFluxes are up-to-date
    Outputs: Updates sln.faceFluxes and sln.cellResiduals
    Returns: Updated sln.cellResiduals

    Applies classical JST method: central differencing + JST artificial diffusion. Each face treated as a 1D problem
    http://aero-comlab.stanford.edu/Papers/jst_2015_updated_07_03_2015.pdf -  see especially pg.5-6
=#
function unstructured_ROEFlux(mesh::Mesh, sln::SolutionState, boundaryConditions, gamma, R)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    nVars = size(sln.cellState, 2)

    # println("Cells for face 154")
    # display(mesh.faces[154])
    # #display(mesh.cell)
    # println("Cell prims on 154 sides")
    # display(sln.cellPrimitives[77,:])
    # display(sln.cellPrimitives[90,:])
    # println("Cell state on 154 sides")
    # display(sln.cellState[77,:])
    # display(sln.cellState[90,:])
    # println("Face 2 owner is:")
    # display(mesh.cells[1])

    # println("Face 1 and 2 area vecs:")
    # display(mesh.fAVecs[1,:])
    # display(mesh.fAVecs[2,:])

    #### 0. Reset all face fluxes
    nLoop = size(sln.faceFluxes, 2)
    for f in 1:nFaces-nBdryFaces
        for flux in 1:nLoop
            sln.faceFluxes[f,flux] = 0
        end
    end

    #### 1. Apply boundary conditions ####
    for b in 1:nBoundaries
        bFunctionIndex = 2*b-1
        boundaryConditions[bFunctionIndex](mesh, sln, b, boundaryConditions[bFunctionIndex+1])
    end

    #### 2. MUSCL difference fluxes ####
    consOwner, consNeighbour = MUSCL_difference(mesh, sln)

    # println("Cons at face 154:")
    # display(consOwner[154,:])
    # display(consNeighbour[154,:])

    #ownerPrims = zeros(nFaces-nBdryFaces, 5)
    #neighbourPrims = zeros(nFaces-nBdryFaces, 5)

    ownerPrims = decodePrimitives3D(consOwner)
    neighbourPrims = decodePrimitives3D(consNeighbour)

    # println("Prims at face 154:")
    # display(ownerPrims[154,:])
    # display(neighbourPrims[154,:])

    #### 3. Find ROE averaged quantites

    roeAveraged = calculateROEAveraged(ownerPrims, neighbourPrims, consOwner, consNeighbour)

    # println("Roe averaged props at face 154")
    # display(roeAveraged[154,:])

    #### 4. Find (and correct) eigenvalues

    eig1, eig2, eig3 = findEigenvalues(roeAveraged, ownerPrims, neighbourPrims)


    #### 5. Find flux vectors (at interior faces???)

    FL, FR, F1, F2, F3 = findFluxVectors(mesh, sln, roeAveraged, eig1, eig2, eig3, ownerPrims, neighbourPrims, consOwner, consNeighbour) #TODO: Pretty sure this fcn isn't treating boundaries correctly

    # println("Theres are $nVars vars")
    # println("$breakdown")

    # println("These are face fluxes at 3492:")
    # display(sln.faceFluxes[3492,:])
    # display(mesh.fAVecs[3492,:])



    @inbounds @fastmath for f in 1:nFaces-nBdryFaces
        ownerCell = mesh.faces[f][1]
        neighbourCell = mesh.faces[f][2]
        d = mesh.cCenters[neighbourCell] .- mesh.cCenters[ownerCell]

        # @views fD = fDeltas[f,:]
        # @views farOwnerfD = fD .- dot(d, fDGrads[ownerCell,:,:])
        # @views farNeighbourfD = fD .+ dot(d, fDGrads[ownerCell,:,:])

        faceFluxes = 0.5 * (FL[f,:] + FR[f,:] - (F1[f,:]+F2[f,:]+F3[f,:]))

        #diffusionFlux = eps2[f]*fD - eps4[f]*(farNeighbourfD - 2*fD + farOwnerfD)
        if f == 2
            println("Face fluxes for face 2 are:")
            display(faceFluxes)
            # println("Individ vectors are:")
            # display(FL[f,:])
            # display(FR[f,:])
            # display(F1[f,:])
            # display(F2[f,:])
            # display(F3[f,:])
            # holder = dot(sln.faceFluxes[f, 1:3], mesh.fAVecs[2])
            # println("Mass flux flow: $holder")
        end

        # Add diffusion flux in component form
        unitFA = normalize(mesh.fAVecs[f])
        for v in 1:nVars
            i1 = (v-1)*3+1
            i2 = i1+2
            sln.faceFluxes[f,i1:i2] .+= (faceFluxes[v] .* unitFA)
        end
    end

    holder = dot(sln.faceFluxes[2, 1:3], mesh.fAVecs[2])
    println("Mass flux flow: $holder")
    display(sln.faceFluxes[2,1:3])



    #### 4. Integrate fluxes at in/out of each cell (sln.faceFluxes) to get change in cell center values (sln.fluxResiduals) ####
    return integrateFluxes_unstructured3D(mesh, sln, boundaryConditions)
end


######################### TimeStepping #################################
# Calculate sln.cellPrimitives and sln.cellFluxes from sln.cellState
function decodeSolution_3D(sln::SolutionState, R=287.05, Cp=1005)
    nCells = size(sln.cellState, 1)
    for c in 1:nCells
        # Updates cell primitives
        @views decodePrimitives3D!(sln.cellPrimitives[c,:], sln.cellState[c,:], R, Cp)
        # Updates mass, xMom, eV2 x,y,z-direction fluxes
        @views calculateFluxes3D!(sln.cellFluxes[c, :], sln.cellPrimitives[c,:], sln.cellState[c,:])
    end
end

#=
    Calculates sln.fluxResiduals from sln.faceFluxes
=#
function integrateFluxes_unstructured3D(mesh::Mesh, sln::SolutionState, boundaryConditions)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

    # Recomputing flux balances, so wipe existing values
    fill!(sln.fluxResiduals, 0)
    nVars = size(sln.fluxResiduals, 2)

    #### Flux Integration ####
    @inbounds @fastmath for f in eachindex(mesh.faces)
        ownerCell = mesh.faces[f][1]
        neighbourCell = mesh.faces[f][2]

        for v in 1:nVars
            i1 = (v-1)*3 + 1
            i2 = i1+2
            @views flow = dot(sln.faceFluxes[f, i1:i2], mesh.fAVecs[f])

            if f == 2
                println("Owner: Flow going from face $f and flux vars $v is: (neg is in)")
                display(flow)
                #display(ownerCell)
            end
            #
            # if neighbourCell == 1
            #     println("Neighbour: Flow is leaving cell 1 through face $f and var $v at a rate of (pos is in)")
            #     display(flow)
            # end

            # Subtract from owner cell
            sln.fluxResiduals[ownerCell, v] -= flow
            # Add to neighbour cell
            if neighbourCell > -1
                sln.fluxResiduals[neighbourCell, v] += flow
            end
        end
    end

    println("Cell flux resids:")
    display(sln.fluxResiduals[1,:])

    # Divide by cell volume
    for c in 1:nCells
        for v in 1:nVars
            sln.fluxResiduals[c,v] /= mesh.cVols[c]
        end
    end

    return sln.fluxResiduals
end

######################### Boundary Conditions #########################
# The job of the BC functions is to calculate/enforce face fluxes at boundary faces.
# The JST method only calculates fluxes at internal faces, then these conditions are applied to calculate them at the boundaries

# InletConditions: [ Static Pressure, Static Temperture, Ux, Uy, Uz, Cp ]
function supersonicInletBoundary!(mesh::Mesh, sln::SolutionState, boundaryNumber, inletConditions)
    P, T, Ux, Uy, Uz, Cp = inletConditions
    rho = idealGasRho(T, P)
    xMom, yMom, zMom = [Ux, Uy, Uz] .* rho
    e = calPerfectEnergy(T, Cp)
    eV2 = rho*(e + (mag([Ux, Uy, Uz])^2)/2)

    boundaryFluxes = Vector{Float64}(undef, 15)
    calculateFluxes3D!(boundaryFluxes, [P, T, Ux, Uy, Uz],  [rho, xMom, yMom, zMom, eV2])
    currentBoundary = mesh.boundaryFaces[boundaryNumber]
    @inbounds for face in currentBoundary
        sln.faceFluxes[face,:] = boundaryFluxes
        sln.faceState[face,:] = [rho, xMom, yMom, zMom, eV2]
        sln.facePrimitives[face,:] = [P, T, Ux, Uy, Uz]
    end
end

# InletConditions: [ totalPressure, totalTemp, nx, ny, nz, gamma, R, Cp ]
# Where n is the unit vector representing the direction of inlet velocity
# Using method from FUN3D solver
# TODO: Allow for variation of R and gamma
function subsonicInletBoundary!(mesh, sln::SolutionState, boundaryNumber, inletConditions)
    Pt, Tt, nx, ny, nz, gamma, R, Cp = inletConditions

    boundaryFluxes = Vector{Float64}(undef, 15)
    primitives = Vector{Float64}(undef, 5)
    state = Vector{Float64}(undef, 5)
    velUnitVector = [ nx, ny, nz ]
    currentBoundary = mesh.boundaryFaces[boundaryNumber]

    @inbounds @fastmath for face in currentBoundary
        ownerCell = mesh.faces[face][1]
        @views adjustedVelocity = dot(velUnitVector, sln.cellPrimitives[ownerCell, 3:5])
        primitives[2] = Tt - (gamma-1)/2 * (adjustedVelocity^2)/(gamma*R)
        machNum = abs(adjustedVelocity) / sqrt(gamma * R * primitives[2])
        primitives[1] = Pt*(1 + (gamma-1)/2 * machNum^2)^(-gamma/(gamma-1))
        primitives[3:5] = adjustedVelocity .* velUnitVector

        # Calculate state variables
        state[1] = idealGasRho(primitives[2], primitives[1], R)
        state[2:4] .= primitives[3:5] .* state[1]
        e = calPerfectEnergy(primitives[2], Cp, R)
        state[5] = state[1]*(e + (mag(primitives[3:5])^2)/2 )

        calculateFluxes3D!(boundaryFluxes, primitives, state)
        sln.faceFluxes[face, :] = boundaryFluxes
        sln.faceState[face, :] = state
        sln.facePrimitives[face, :] = primitives
    end
end

function pressureOutletBoundary!(mesh, sln::SolutionState, boundaryNumber, outletPressure)
    nFluxes = size(sln.cellFluxes, 2)

    # Directly extrapolate cell center flux to boundary (zero gradient between the cell center and the boundary)
    currentBoundary = mesh.boundaryFaces[boundaryNumber]
    @inbounds @fastmath for face in currentBoundary
        ownerCell = mesh.faces[face][1]
        # Mass fluxes are unaffected by pressure boundary
        for flux in 1:nFluxes
            sln.faceFluxes[face, flux] = sln.cellFluxes[ownerCell, flux]
        end

        origP = sln.cellPrimitives[ownerCell, 1]
        # Adjust the momentum fluxes containing pressure
        for flux in [4, 8, 12]
            sln.faceFluxes[face, flux] += outletPressure - origP
        end
        #Adjust the energy fluxes
        for d in 1:3
            sln.faceFluxes[face, 12+d] += sln.cellPrimitives[ownerCell, 2+d]*(outletPressure - origP)
        end


        #Use cell center values for non-pressure. Alternative is to interpolate and use fluxes?
        prims = zeros(1,5)
        prims[1,:] = sln.cellPrimitives[ownerCell, :]
        prims[1,1] = outletPressure
        sln.faceState[face,:] = encodePrimitives3D(prims)
        sln.facePrimitives[face,:] = prims
    end
end

function zeroGradientBoundary!(mesh::Mesh, sln::SolutionState, boundaryNumber, _)
    nFluxes = size(sln.cellFluxes, 2)

    # Directly extrapolate cell center flux to boundary (zero gradient between the cell center and the boundary)
    currentBoundary = mesh.boundaryFaces[boundaryNumber]
    @inbounds for face in currentBoundary
        ownerCell = mesh.faces[face][1] #One of these will be -1 (no cell), the other is the boundary cell we want
        for flux in 1:nFluxes
            sln.faceFluxes[face, flux] = sln.cellFluxes[ownerCell, flux]
        end

        #Zero grad, so all cell center values should be the same at the face
        prims = zeros(1,5)
        prims[1,:] = sln.cellPrimitives[ownerCell,:]
        sln.faceState[face, :] = encodePrimitives3D(prims)
        sln.facePrimitives[face, :] = prims
    end
end

function wallBoundary!(mesh::Mesh, sln::SolutionState, boundaryNumber, _)
    currentBoundary = mesh.boundaryFaces[boundaryNumber]
    @inbounds for f in currentBoundary
        ownerCell = max(mesh.faces[f][1], mesh.faces[f][2]) #One of these will be -1 (no cell), the other is the boundary cell we want

        unitFA = normalize(mesh.fAVecs[f])

        faceP = sln.cellPrimitives[ownerCell, 1]
        faceT = sln.cellPrimitives[ownerCell, 2]
        # Momentum flux is Pressure in each of the normal directions (dot product)
        sln.faceFluxes[f, 4] = faceP #* -unitFA[1]
        sln.faceFluxes[f, 8] = faceP #* -unitFA[2]
        sln.faceFluxes[f, 12] = faceP #* -unitFA[3]

        # Mass Flux is zero
        sln.faceFluxes[f, 1:3] .= 0.0
        # Energy Flux is zero
        sln.faceFluxes[f, 13:15] .= 0.0

        # For states, want to find conserved quantities with normal velocity = 0
        #unitFA = normalize(mesh.fAVecs[f])
        cellVelo = sln.cellPrimitives[ownerCell, 3:5]

        faceNormalVelo = zeros(3)
        faceTangVelo = zeros(3)
        for d in 1:3
            faceNormalVelo[d] = dot(cellVelo, unitFA) * mesh.fAVecs[f][d]
            faceTangVelo[d] = cellVelo[d] - faceNormalVelo[d]
        end

        prims = zeros(1,5)
        prims[1,1] = faceP
        prims[1,2] = faceT
        prims[1,3:5] = faceTangVelo

        sln.faceState[f, :] = encodePrimitives3D(prims)
        sln.facePrimitives[f, :] = prims

    end
end
symmetryBoundary! = wallBoundary!

function emptyBoundary!(mesh::Mesh, sln::SolutionState, boundaryNumber, _)
    return
end

######################### Numerical Integration #########

# xFcn should be the gradients at every cell
function gaussQuadrature(mesh::AdaptMesh, sln::SolutionState, xFcn, points=15, )
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    if points == 15
        eta = [0.0000,-0.2011940939974345, 0.2011940939974345,-0.3941513470775634, 0.3941513470775634, -0.5709721726085388, 0.5709721726085388, -0.7244177313601701, 0.7244177313601701, -0.8482065834104272, 0.8482065834104272, -0.9372733924007060, 0.9372733924007060, -0.9879925180204854, 0.9879925180204854]
        w = [0.2025782419255613, 0.1984314853271116, 0.1984314853271116, 0.1861610000155622, 0.1861610000155622, 0.1662692058169939, 0.1662692058169939, 0.1395706779261543, 0.1395706779261543, 0.1071592204671719, 0.1071592204671719, 0.0703660474881081, 0.0703660474881081, 0.0307532419961173, 0.0307532419961173]
    else
        println("Need to add abscissa and weights for other points!")
    end

    sum = zeros(nCells)

    #Need to set it to a mesh value at initiliaze so that it doesn't mistakenly just stick with a non mesh value
    xMin = mesh.fPoints[1,1]
    xMax = mesh.fPoints[1,1]
    yMin = mesh.fPoints[1,2]
    yMax = mesh.fPoints[1,2]
    zMin = mesh.fPoints[1,3]
    zMax = mesh.fPoints[1,3]

    # For each cell, go through each
    for c in 1:nCells
        nFacesCel = size(mesh.cells[c], 1)

        xMin, xMax, yMin, yMax, zMin, zMax =setCellScalingLims!(xMin, xMax, yMin, yMax, zMin, zMax, mesh, c)

        limits = [xMin, yMin, zMin, xMax, yMax, zMax]

        for d in 1:3
            for i in 1:points
                dim = convertEtaToDim(eta[i], limits[d], limits[d+3])

                # ?????????????????????????
                sum[c] += xFcn[c,1,d]*(d_d)
            end
        end




    end

end

function setCellScalingLims!(xMin, xMax, yMin, yMax, zMin, zMax, mesh::AdaptMesh, c)
    for f in mesh.cells[c]
        for p in mesh.fPoints[f]
            if mesh.fPointLocs[p,1] < xMin
                xMin = mesh.fPointLocs[p,1]
            elseif mesh.fPointLocs[p,1] > xMax
                xMax = mesh.fPointLocs[p,1]
            end
            if mesh.fPointLocs[p,2] < yMin
                yMin = mesh.fPointLocs[p,2]
            elseif mesh.fPointLocs[p,2] > yMax
                yMax = mesh.fPointLocs[p,2]
            end
            if mesh.fPointLocs[p,3] < zMin
                zMin = mesh.fPointLocs[p,3]
            elseif mesh.fPointLocs[p,3] > zMax
                zMax = mesh.fPointLocs[p,3]
            end
        end
    end

    return xMin, xMax, yMin, yMax, zMin, zMax
end

function convertEtaToDim(eta, l, u)

end

function convertDimToEta()
end

function findCellsToAdapt(error; adaptPercent=0.05)
    nCells = size(error, 1)
    nCellsAdapt = floor(Int, nCells * adaptPercent)

    nAdaptList = zeros(Int, nCellsAdapt)

    for i in 1:nCellsAdapt
        (x, index) = findmax(error)
        nAdaptList[i] = index
        error[index] = 0
    end
    return nAdaptList
end

function interpSlnToNewMesh(sln::SolutionState, newSln::SolutionState, oldNCells, nAdaptList, newCellsList)
    #oldNCells = size(mesh.cells,1)

    # Then loop through all the old cells (full loop, but exclude adapted cells) and map prims and conserved to new sln file
    #   need a way to map the old cell numbers to the new ones? I think it's newCellIndex = oldIndex- cellsDeletedInFrontOfIt
    for i in 1:oldNCells
        if !any(j->(j==i), nAdaptList) #Skip the cells we deleted for now, they will have a different treatment

            iNew = i
            counter = 0
            for j in 1:size(nAdaptList,1)
                if i >= nAdaptList[j]
                    iNew -= 1
                    counter += 1
                end
            end


            newSln.cellState[iNew,:] = sln.cellState[i,:]
            newSln.cellFluxes[iNew,:] = sln.cellFluxes[i,:]
            newSln.cellPrimitives[iNew,:] = sln.cellPrimitives[i,:]
            newSln.fluxResiduals[iNew,:] = sln.fluxResiduals[i,:]
        end
    end

    for i in 1:size(nAdaptList,1)
        oldCell = nAdaptList[i]
        for j in 1:size(newCellsList[i],1)
            newCell = newCellsList[i][j]
            newSln.cellState[newCell,:] = sln.cellState[oldCell,:]
            newSln.cellFluxes[newCell,:] = sln.cellFluxes[oldCell,:]
            newSln.cellPrimitives[newCell,:] = sln.cellPrimitives[oldCell,:]
            newSln.fluxResiduals[newCell,:] = sln.fluxResiduals[oldCell,:]
        end
    end


    # Map to faces done from BCs outside this, and linInterp_3D (after BCs)

    return newSln
end



######################### Solvers #######################
#=
    This is where the magic happens!
    Do CFD!

    Arguments:
        mesh: see dataStructuresDefinitions.md / dataStructures.jl
        meshPath:           (string) path to folder where mesh is stored
        cellPrimitives:     initial cell-center conditions, see dataStructuresDefinitions.md / dataStructures.jl
        boundaryConditions: Array of alternating boundary condition function references and associated boundary condition parameters:
            Ex: [ emptyBoundary, [], supersonicInletBoundary, [P, T, Ux, Uy, Uz, Cp], zeroGradientBoundary, [], symmetryBoundary, [], zeroGradientBoundary, [], wallBoundary, [] ]
            Order of BC's must match the order of boundaries defined in the mesh's 'boundaries' file
        timeIntegrationFn:  The desired time integration function from timeDiscretizations.jl (should update sln.cellState, sln.cellPrimitives, and sln.cellFluxes based on their old values and sln.cellResiduals (obtained from fluxFunction))
        fluxFunction:       The desired function to calculate sln.cellResiduals from sln.cellState, sln.cellPrimitives, and sln.cellFluxes
        initDt:             Initial time step (s)
        endTime:            Simulation End Time (s). Start Time = 0
        outputInterval:     Writes solution/restart files every outputInterval seconds of simulation time
        targetCFL:          Target maximum CFL in the computation domain (time step will be adjusted based on this value)
        gamma/R/Cp:         Fluid properties
        silent:             Controls whether progress is written to console (can slow down simulation slightly for very simple cases)
        restart:            Controls whether restarting from a restart file
        createRestartFile:  Controls whether to write restart files. These are overwritten every time they are outputted.
        createVTKOutput:    Controls whether to write .vtk output. These are not overwritten every time they are outputted.
        restartFiles:       (string) path to restart file to read/write from

    Returns:
        sln.cellPrimitives at end of simulation

    Outputs:
        restart file (if specified)
        .vtk files (is specified)
=#
function unstructured3DFVM(mesh::Mesh, meshPath, cellPrimitives::Matrix{Float64}, boundaryConditions; timeIntegrationFn=forwardEuler,
        fluxFunction=unstructured_JSTFlux, initDt=0.001, endTime=0.14267, outputInterval=0.01, adaptInterval=0.005,targetCFL=0.2, gamma=1.4, R=287.05, Cp=1005,
        silent=true, restart=false, createRestartFile=true, createVTKOutput=true, restartFile="JuliaCFDRestart.txt")

    if !silent
        println("Initializing Simulation")
    end
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    nDims = 3

    # Each dimension adds one momentum equation
    nVars = 2+nDims
    # Each dimension adds a flux for each conserved quantity
    nFluxes = nVars*nDims

    initP = cellPrimitives[1,1]
    initT = cellPrimitives[1,2]
    initU = cellPrimitives[1,3:5]

    if restart
        if !silent
            println("Reading restart file: $restartFile")
        end

        cellPrimitives = readRestartFile(restartFile)

        # Check that restart file has the same nnumber of cells as the mesh
        nCellsRestart = size(cellPrimitives, 1)
        if nCellsRestart != nCells
            throw(ArgumentError("Number of cells in restart file ($nCellsRestart) does not match the present mesh ($nCells). Please provide a matching mesh and restart file."))
        end
    end

    # rho, xMom, total energy from P, T, Ux, Uy, Uz
    cellState = encodePrimitives3D(cellPrimitives, R, Cp)
    cellFluxes = zeros(nCells, nFluxes)
    fluxResiduals = zeros(nCells, nVars)
    faceState = zeros(nFaces, nVars)
    faceFluxes = zeros(nFaces, nFluxes)
    facePrimitives = zeros(nFaces, nVars)
    # Initialize solution state
    sln = SolutionState(cellState, cellFluxes, cellPrimitives, fluxResiduals, faceState, faceFluxes, facePrimitives)

    # Calculates cell fluxes, primitives from cell state
    decodeSolution_3D(sln, R, Cp)

    if !silent
        println("Starting iterations")
    end

    dt = initDt
    if timeIntegrationFn==LTSEuler
        dt = zeros(nCells)
    end

    ## Inital solve of boundary conditions, to prevent CFL from failing on them
    for b in 1:nBoundaries
        bFunctionIndex = 2*b-1
        boundaryConditions[bFunctionIndex](mesh, sln, b, boundaryConditions[bFunctionIndex+1])
    end


    currTime = 0
    timeStepCounter = 0
    nextOutputTime = outputInterval
    nextAdaptTime = adaptInterval
    writeOutputThisIteration = false
    adaptMeshThisIteration = false
    CFLvec = zeros(nCells)
    while currTime < endTime
        ############## Timestep adjustment #############
        if timeIntegrationFn == LTSEuler
            CFL = 1.0
            if timeStepCounter < 10
                # Ramp up CFL linearly in the first ten time steps to reach the target CFL
                CFL = (timeStepCounter+1) * targetCFL / 10
            else
                CFL = targetCFL
            end

            # Store the target CFL for the present time step in the first element of the dt vector.
            # The time discretization function will determine the actual local time step based on this target CFL
            dt[1] = CFL
        else
            CFL!(CFLvec, mesh, sln, dt, gamma, R)

            CFL = maximum(CFLvec)

            # If CFL too high, attempt to preserve stability by cutting timestep size in half
            if CFL > targetCFL*1.01
                dt *= targetCFL/(2*CFL)
            # Otherwise slowly approach target CFL
            else
                dt *= ((targetCFL/CFL - 1)/10+1)
            end

            # Cut timestep to hit endtime/writeTimes as necessary
            if (endTime - currTime) < dt
                dt = endTime - currTime
            elseif (nextOutputTime - currTime) < dt
                dt = nextOutputTime - currTime
                writeOutputThisIteration = true
            end

            if (nextAdaptTime - currTime) < dt
                adaptMeshThisIteration = true
            end
        end

        ############## Take a timestep #############
        sln = timeIntegrationFn(mesh, fluxFunction, sln, boundaryConditions, gamma, R, Cp, dt)


        if timeIntegrationFn == LTSEuler
            currTime += CFL
            if (nextOutputTime - currTime) < CFL
                writeOutputThisIteration = true
            end
        else
            currTime += dt
        end
        timeStepCounter += 1

        if !silent
            @printf("Timestep: %5.0f, simTime: %9.4g, Max CFL: %9.4g \n", timeStepCounter, currTime, CFL)
        end

        if writeOutputThisIteration
            updateSolutionOutput(sln.cellPrimitives, restartFile, meshPath, createRestartFile, createVTKOutput)
            writeOutputThisIteration = false
            nextOutputTime += outputInterval
        end

        if adaptMeshThisIteration
            #updateSolutionOutput(sln.cellPrimitives, restartFile, meshPath, createRestartFile, createVTKOutput)
            oldNCells = nCells
            println("Adapting mesh...")

            gradP_LS = leastSqGrad(mesh, sln, 3, 1)

            gradP_mag = zeros(nCells)
            for i in 1:nCells
                gradP_mag[i] = mag(gradP_LS[i, 1, :])
            end

            nAdaptList = findCellsToAdapt(gradP_mag, adaptPercent=0.15)

            origAdaptList = copy(nAdaptList)

            adaptMesh = createMeshAdaptStruct(mesh, meshPath)

            println("Writing new mesh...")

            adaptedMesh, newCellsList, facesData = adaptNewMesh(adaptMesh, nAdaptList, boundaryConditions) # NOTE: adaptedMesh needs to be written to OF file and then read back in to be used
            if size(adaptedMesh.cells,1) > 1095
                println("1095 debug")
                println("Faces bdry indices")
                display(facesData.bdryIndices)
                println("Faces in cell 1095")
                display(adaptedMesh.cells[1095])
                for f in adaptedMesh.cells[1095]
                    display(adaptMesh.fPoints[f])
                end
            end

            # println("Cell 1020")
            # display(adaptedMesh.cells[1010:end])
            # #println("$breakdown")
            # nCells = size(adaptedMesh.cells, 1)
            # counter = 0
            # maxCell = 0
            # faces = adaptedMesh.faces
            # for face in faces
            #     for i in 1:2
            #         if face[i] > nCells
            #             counter += 1
            #         end
            #         if face[i] > maxCell
            #             maxCell = face[i]
            #         end
            #     end
            # end
            # println("Extent of problem: $counter count and $maxCell max")
            # println("FacesData: $facesData")

            newPath = writeNewOpenFOAMMesh(adaptedMesh, facesData)

            newMesh = OpenFOAMMesh(newPath)

            #Interp to new cells

            println("Interpolating solution to new mesh...")

            newNCells = size(newMesh.cells,1)
            newNFaces = size(newMesh.faces,1)

            # rho, xMom, total energy from P, T, Ux, Uy, Uz
            newCellPrims = initializeUniformSolution3D(newMesh, initP, initT, initU...)



            cellState = encodePrimitives3D(newCellPrims)
            cellFluxes = zeros(newNCells, nFluxes)
            fluxResiduals = zeros(newNCells, nVars)
            faceState = zeros(newNFaces, nVars)
            faceFluxes = zeros(newNFaces, nFluxes)
            facePrimitives = zeros(newNFaces, nVars)
            # Initialize solution state
            newSln = SolutionState(cellState, cellFluxes, newCellPrims, fluxResiduals, faceState, faceFluxes, facePrimitives)

            #Replace all the old variables
            meshPath = newPath
            sln = interpSlnToNewMesh(sln, newSln, oldNCells, origAdaptList, newCellsList)
            mesh = newMesh
            CFLvec = zeros(size(sln.cellState,1))

            nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

            println("After adaptation, we have $nCells cells, and $nFaces faces")

            # holder1 = size(mesh.cells,1)
            # holder2 = size(mesh.faces,1)
            # #holder3 = size(mesh.fPoints,1)
            # #holder4 = size(mesh.fPLocs,1)
            # #holder5 = size(nList,1)
            # #holder6 = cells[631:635]
            #
            # display(mesh.faces[3931])
            # println("Debug statement: cells $holder1, faces $holder2,")# fPoints $holder3, fPLocs $holder4, nList $holder5")
            #
            # display(mesh.cells[1000:end])
            #println("$breakdown")
            # nCells = size(sln.cellPrimitives,1)
            # println("New solution size: $nCells")


            #TODO: I don't think this interpolation to the faces is required - check if it's redundant
            # Apply BCs to sln
            #### 1. Apply boundary conditions ####
            for b in 1:nBoundaries
                bFunctionIndex = 2*b-1
                boundaryConditions[bFunctionIndex](mesh, sln, b, boundaryConditions[bFunctionIndex+1])
            end

            #### 2. Interp all values to faces ####
            sln.faceState[1:nFaces-nBdryFaces,:] = linInterpToFace_3D(mesh,sln, 1, -1)
            sln.faceFluxes[1:nFaces-nBdryFaces,:] = linInterpToFace_3D(mesh, sln, 2, -1)
            sln.facePrimitives[1:nFaces-nBdryFaces,:] = linInterpToFace_3D(mesh, sln, 3, -1)


            #Interp from cells to faces

            #println("Wrote a new mesh!")
            #println("$breakdown")

            #println("Returned from adapting!")
            #print("$breakdown")

            #sln.cellPrimitives[:,1] = gradP_LS[:,:,1]
            #sln.cellPrimitives[:,2] = gradP_GG[:,:,1]
            #sln.cellPrimitives[:,3] = 100 * (gradP_LS[:,:,1]-gradP_GG[:,:,1]) ./ gradP_LS[:,:,1]

            createRestartFile = false


            adaptMesh = createMeshAdaptStruct(mesh, meshPath)

            if nCells > 1095
                println("1095 debug")
                println("Faces bdry indices")
                display(facesData.bdryIndices)
                println("Faces in cell 1095")
                display(mesh.cells[1095])
                for f in mesh.cells[1095]
                    display(adaptMesh.fPoints[f])
                end
            end




            updateSolutionOutput(sln.cellPrimitives, restartFile, meshPath, createRestartFile, createVTKOutput)

            println("Wrote new sln file!")

            #println("This is how you $break")

            firstPassFlag = false
            adaptMeshThisIteration = false
            nextAdaptTime += adaptInterval


        end
    end

    updateSolutionOutput(sln.cellPrimitives, restartFile, meshPath, createRestartFile, createVTKOutput)
    return sln.cellPrimitives
end
