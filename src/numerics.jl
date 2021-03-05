######################### Gradient Computation #######################
#TODO: LSQ Gradient
# Non-functional
function leastSqGrad(mesh::Mesh, matrix::AbstractArray{Float64, 2}, stencil=zeros(2,2))
    # Stencil should be a list of lists, with each sublist containing the cells contained in the stencil of the main cell


    #cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    #bdryFaceIndices = Array(nFaces-nBdryFaces:nFaces)

    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    nVars = size(matrix, 2) # Returns the length of the second dimension of "matrix"

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

        dx = mesh.cCenter[c2][1] - mesh.cCenter[c1][1]
        dy = mesh.cCenter[c2][2] - mesh.cCenter[c1][2]
        dz = mesh.cCenter[c2][3] - mesh.cCenter[c1][3]

        weight = 1 / sqrt(dx*dx + dy*dy + dz*dz)

        wdx = weight * dx
        wdy = weight * dy
        wdz = weight * dz

        for v in 1:nVars
            dv = matrix[c2,v] - matrix[c1,v]
            wdv = weight * dv

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

        end

    end

    # Deal with boundary faces, and add them to the matrix vectors



    return grad
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
    if valuesAtFaces
        faceVals = matrix
    else
        faceVals = linInterp_3D(mesh, matrix)
    end

    # Create matrix to hold gradients
    grad = zeros(nCells, nVars, 3)

    # Integrate fluxes from each face
    @fastmath for f in eachindex(mesh.faces)
        ownerCell = mesh.faces[f][1]
        neighbourCell = mesh.faces[f][2]

        @inbounds for v in 1:nVars
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
    @simd for c in 1:nCells
        for v in 1:nVars
            for d in 1:3
                grad[c,v,d] /= mesh.cVols[c]
            end
        end
    end

    return grad
end

####################### Face value interpolation ######################
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
function linInterp_3D(mesh::Mesh, matrix::AbstractArray{Float64, 2}, faceVals=nothing)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    nVars = size(matrix, 2)

    # Make a new matrix if one is not passed in
    if faceVals == nothing
        faceVals = zeros(nFaces, nVars)
    end

    # Boundary face fluxes must be set separately
    for f in 1:nFaces-nBdryFaces
        # Find value at face using linear interpolation
        c1 = mesh.faces[f][1]
        c2 = mesh.faces[f][2]

        #TODO: Precompute these distances?
        c1Dist = 0
        c2Dist = 0
        for i in 1:3
            c1Dist += (mesh.cCenters[c1][i] - mesh.fCenters[f][i])^2
            c2Dist += (mesh.cCenters[c2][i] - mesh.fCenters[f][i])^2
        end
        totalDist = c1Dist + c2Dist

        for v in 1:nVars
            faceVals[f, v] = matrix[c1, v]*(c2Dist/totalDist) + matrix[c2, v]*(c1Dist/totalDist)
        end
    end

    return faceVals
end

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

            # Subtract from owner cell
            sln.fluxResiduals[ownerCell, v] -= flow
            # Add to neighbour cell
            if neighbourCell > -1
                sln.fluxResiduals[neighbourCell, v] += flow
            end
        end
    end

    # Divide by cell volume
    for c in 1:nCells
        for v in 1:nVars
            sln.fluxResiduals[c,v] /= mesh.cVols[c]
        end
    end

    return sln.fluxResiduals
end
