using Printf
include("constitutiveRelations.jl")
include("vectorFunctions.jl")
include("timeDiscretizations.jl")
include("JST_structured_1D.jl")
include("mesh.jl")
include("output.jl")

__precompile__()

######################### Initialization ########################
# Returns cellPrimitives matrix
function initializeUniformSolution3D(mesh, P, T, Ux, Uy, Uz)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

    cellPrimitives = zeros(nCells, 5)
    for c in 1:nCells
        cellPrimitives[c, :] = [ P, T, Ux, Uy, Uz ]
    end

    return cellPrimitives
end

######################### CFL ########################
# TODO: Generalize cell size
# TODO: 3D definition of CFL: Sum up in all directions
function CFL(U, T, dt, dx, gamma=1.4, R=287.05)
    return (abs(U[1]) + sqrt(gamma * R * T)) * dt / dx
end

function maxCFL3D(mesh, solutionState, dt, gamma=1.4, R=287.05)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState

    maxCFL = 0
    for c in 1:nCells
        cDeltaT = ([ abs(cellPrimitives[c,i]) for i in 3:5 ] .+ sqrt(gamma * R * cellPrimitives[c,2])) .* dt
        #TODO: Precompute
        maxCoords = [ -1000000.0, -1000000.0, -1000000.0 ]
        minCoords = [ 1000000.0, 1000000.0, 1000000.0 ]
        for f in cells[c]
            for d in 1:3
                maxCoords[d] = max(maxCoords[d], fCenters[f][d])
                minCoords[d] = min(minCoords[d], fCenters[f][d])
            end
        end
        CFL = 0.0
        for d in 1:3
            dx = (maxCoords[d] - minCoords[d])
            if dx > 0
                CFL += cDeltaT[d] / dx
            end
        end
        maxCFL = max(maxCFL, CFL)
    end

    return maxCFL
end

######################### Gradient Computation #######################
#TODO: make all these function return a single array if you pass in a single value
function leastSqGrad(mesh, values...)
    #TODO
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    bdryFaceIndices = Array(nFaces-nBdryFaces:nFaces)

    result = []
    for vals in values
    end
    return result
end

# Pass in values that have already been interpolated to faces to avoid re-interpolating
# Returns array of 3-D vectors
function greenGaussGrad(mesh, valuesAtFaces=false, values...)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    bdryFaceIndices = Array(nFaces-nBdryFaces:nFaces)

    result = []
    for vals in values
        #Interpolate values to faces
        if valuesAtFaces != true
            faceVals = linInterp(mesh, vals)[1]
        else
            faceVals = vals
        end

        # Initialize gradients to zero
        grad = Array{Array{Float64, 1}, 1}(undef, nCells)
        for c in 1:nCells
            grad[c] = [0, 0, 0]
        end

        # Integrate fluxes from each face
        for f in 1:nFaces-nBdryFaces
            faceIntegral = fAVecs[f] .* faceVals[f]

            ownerCell = faces[f][1]
            neighbourCell = faces[f][2]

            grad[ownerCell] += faceIntegral
            grad[neighbourCell] -= faceIntegral
        end

        # Divide integral by cell volume to obtain gradients
        for c in 1:nCells
            grad[c] = grad[c] / cVols[c]
        end

        # Set boundary gradients to zero
        for f in nFaces-nBdryFaces+1:nFaces
            for cell in faces[f]
                if cell != -1
                    grad[cell] = [0, 0, 0]
                end
            end
        end

        push!(result, grad)
    end
    return result
end

# Will take the gradient of (scalar) data provided in matrix form
# Cell      x1      x2      x3
# Cell 1    x1_1    x2_1    x3_1
# Cell 2    x1_2    x2_2    x3_2
# ...
# and output a three-dimensional gradient matrix of the following form
# Cell      x1          x2          x3
# Cell 1    grad(x1)_1  grad(x2)_1  grad(x3)_1
# Cell 2    grad(x1)_2  grad(x2)_2  grad(x3)_2
# ...
# Where each grad(xa)_b is made up of three elements for the (x,y,z) directions
function greenGaussGrad_matrix(mesh, matrix, valuesAtFaces=false)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    bdryFaceIndices = Array(nFaces-nBdryFaces:nFaces)
    nVars = size(matrix, 2)

    # Create matrix to hold gradients
    grad = zeros(nCells, nVars, 3)

    #Interpolate values to faces
    if valuesAtFaces != true
        faceVals = linInterp_3D(mesh, matrix)
    else
        faceVals = matrix
    end

    # Integrate fluxes from each face
    faceIntegral = zeros(nVars, 3)
    for f in 1:nFaces
        for v in 1:nVars
            faceIntegral[v,:] = fAVecs[f] .* faceVals[f, v]
        end

        ownerCell = faces[f][1]
        neighbourCell = faces[f][2]

        if ownerCell > -1
            grad[ownerCell, :, :] += faceIntegral
        end
        if neighbourCell > -1
            grad[neighbourCell, :, :] -= faceIntegral
        end
    end

    # Divide integral by cell volume to obtain gradients
    for c in 1:nCells
        grad[c,:,:] ./= cVols[c]
    end

    return grad
end

######################### Laplacian Computation #######################
# Calculates the vector Ef for use in the laplacian calculation, as described in Moukalled pg. 241-244
function laplacian_Ef(nonOrthoCorrection, Sf, n, CF, e)
    if nonOrthoCorrection == "None"
        return Sf
    elseif nonOrthoCorrection == "MinCorr"
        return dot(e, Sf) .* e
    elseif nonOrthoCorrection == "OrthogCorr"
        return mag(Sf) .* e
    elseif nonOrthoCorrection == "OverRelax"
        return (dot(Sf, Sf) / dot(e, Sf)) .* e
    end
end

function laplacian_FaceFlux(mesh, nonOrthoCorrection, face, interpolatedFaceGrad, orthoFaceGrad)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh

    ownerCell = faces[face][1]
    neighbourCell = faces[face][2]

    Sf = fAVecs[face]
    n = normalize(Sf)
    CF = cCenters[neighbourCell] .- cCenters[ownerCell]
    e = normalize(CF)

    Ef = laplacian_Ef(nonOrthoCorrection, Sf, n, CF, e)
    Tf = Sf .- Ef

    # Orthogonal contribution + non-ortho correction term
    return mag(Ef)*orthoFaceGrad + dot(interpolatedFaceGrad, Tf)
end

# Used as orthogonal contribution in calculation of laplacian
function ortho_FaceGradient(mesh, values...)
    result = []

    # Interpret mesh
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    bdryFaceIndices = Array(nFaces-nBdryFaces:nFaces)

    for vals in values
        orthoGrad = Array{Float64, 1}(undef, nFaces)

        # Integrate fluxes from each face
        for f in 1:nFaces-nBdryFaces
            ownerCell = faces[f][1]
            neighbourCell = faces[f][2]

            #TODO: Precompute these distances
            distance = mag(cCenters[neighbourCell] .- cCenters[ownerCell])
            orthoGrad[f] = (vals[neighbourCell] .- vals[ownerCell]) ./ distance
        end

        # Set boundary gradients to zero
        for f in nFaces-nBdryFaces+1:nFaces
            orthoGrad[f] = 0.0
        end

        push!(result, orthoGrad)
    end
    return result
end

# Non-ortho correction options are: (all from Moukalled)
#   None
#   MinCorr (Section 8.6.2)
#   OrthogCorr (Section 8.6.3)
#   OverRelax (Seciton 8.6.4)
function laplacian(mesh, nonOrthoCorrection="None", values...)
    result = []

    # Interpret mesh
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    bdryFaceIndices = Array(nFaces-nBdryFaces:nFaces)

    for vals in values
        ###### Precomputations ######
        # Compute gradients
        grads = greenGaussGrad(mesh, false, vals)[1]
        orthoFaceGrads = ortho_FaceGradient(mesh, vals)[1]

        #Interpolate gradients to faces
        faceGrads = linInterp(mesh, grads)[1]

        ###### Laplacian ######
        # Initialize to zero
        lapl = Array{Float64, 1}(undef, nCells)
        for c in 1:nCells
            lapl[c] = 0.0
        end

        # Integrate gradient fluxes from each face
        for f in 1:nFaces-nBdryFaces
            faceIntegral = laplacian_FaceFlux(mesh, nonOrthoCorrection, f, faceGrads[f], orthoFaceGrads[f])

            ownerCell = faces[f][1]
            neighbourCell = faces[f][2]
            lapl[ownerCell] += faceIntegral
            lapl[neighbourCell] -= faceIntegral
        end

        # Divide integral by cell volume to obtain gradients
        for c in 1:nCells
            lapl[c] = lapl[c] / cVols[c]
        end

        # Set boundary laplacians to zero
        for f in nFaces-nBdryFaces+1:nFaces
            for cell in faces[f]
                if cell != -1
                    lapl[cell] = 0.0
                end
            end
        end

        push!(result, lapl)
    end
    return result
end

####################### Face value interpolation ####################
# Interpolates to all INTERIOR faces
function upwindInterp(mesh, U, values...)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nFaces = size(faces, 1)
    nBdryFaces = size(boundaryFaces, 1)

    # Compute face fluxes to see if they're positive or not
    faceVels = linInterp(mesh, U)[1]

    fFluxes = Array{Float64, 1}(undef, nFaces)
    for i in 1:nFaces-nBdryFaces
        fFluxes[i] = dot(fAVecs[i], faceVels[i])
    end

    result = []
    for vals in values
        fVals = []

        # Use sign of flux to choose between owner or neighbour node values
        for i in 1:nFaces-nBdryFaces
            if fFluxes[i] > 0
                push!(fVals, vals[faces[i][1]])
            elseif fFluxes[i] < 0
                push!(fVals, vals[faces[i][2]])
            else
                push!(fVals, (vals[faces[i][1]] + vals[faces[i][2]])/2)
            end
        end

        for i in 1:nBdryFaces
            push!(fVals, 0)
        end
        push!(result, fVals)
    end
    return result
end

# Array versions
# Handles scalar or vector-valued variables
function linInterp(mesh, values...)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

    result = []
    for vals in values
        fVals = []

        ############ Default Value ############
        # Construct a default value of the appropriate dimensionality
        # defaultVal = 0.0 for 1D, [ 0.0, 0.0 ] for 2D etc...
        defaultVal = []
        if size(vals[1], 1) ==1
            defaultVal = 0
        else
            for a in vals[1]
                push!(defaultVal, 0.0)
            end
        end

        ########### Do Interpolations ###########
        for i in 1:nFaces
            if i > nFaces-nBdryFaces
                # Faces on boundaries not treated, must be set externally
                push!(fVals, defaultVal)

            else
                # Find value at face using linear interpolation
                c1 = faces[i][1]
                c2 = faces[i][2]

                #TODO: Precompute these distances
                c1Dist = mag(cCenters[c1] .- fCenters[i])
                c2Dist = mag(cCenters[c2] .- fCenters[i])
                totalDist = c1Dist + c2Dist
                push!(fVals, vals[c1].*(c2Dist/totalDist) .+ vals[c2].*(c1Dist/totalDist))
            end
        end

        push!(result, fVals)
    end
    return result
end

# Flux interpolation
function linInterp_3D(mesh, solutionState::Array{Array{Float64, 2}, 1})
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    nFluxes = size(cellFluxes, 2)

    # Boundary face fluxes must be set separately
    for f in 1:nFaces-nBdryFaces
        # Find value at face using linear interpolation
        c1 = faces[f][1]
        c2 = faces[f][2]

        #TODO: Precompute these distances
        c1Dist = mag(cCenters[c1] .- fCenters[f])
        c2Dist = mag(cCenters[c2] .- fCenters[f])
        totalDist = c1Dist + c2Dist

        faceFluxes[f, :] .= cellFluxes[c1, :].*(c2Dist/totalDist) .+ cellFluxes[c2, :].*(c1Dist/totalDist)
    end
end

# Arbitrary value matrix interpolation
function linInterp_3D(mesh, matrix::Array{Float64, 2})
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    nVars = size(matrix, 2)

    faceVals = zeros(nFaces, nVars)

    # Boundary face fluxes must be set separately
    for f in 1:nFaces-nBdryFaces
        # Find value at face using linear interpolation
        c1 = faces[f][1]
        c2 = faces[f][2]

        #TODO: Precompute these distances
        c1Dist = mag(cCenters[c1] .- fCenters[f])
        c2Dist = mag(cCenters[c2] .- fCenters[f])
        totalDist = c1Dist + c2Dist

        faceVals[f, :] = matrix[c1, :].*(c2Dist/totalDist) .+ matrix[c2, :].*(c1Dist/totalDist)
    end

    return faceVals
end

function avgInterp(mesh, solutionState)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState
    nFluxes = size(cellFluxes, 2)

    # Boundary face fluxes must be set separately
    for f in 1:nFaces-nBdryFaces
        c1 = faces[f][1]
        c2 = faces[f][2]

        for v in 1:nFluxes
            faceFluxes[f, v] = (cellFluxes[c1, v] + cellFluxes[c2, v])/2
        end
    end
end

# Returns matrix of values
function maxInterp(mesh, vars...)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

    nVars = size(vars, 1)
    faceVals = zeros(nFaces, nVars)

    # Boundary face fluxes must be set separately
    for f in 1:nFaces-nBdryFaces
        c1 = faces[f][1]
        c2 = faces[f][2]

        for v in 1:nVars
            faceVals[f, v] = max(vars[v][c1], vars[v][c2])
        end
    end

    return faceVals
end

function faceDeltas(mesh, solutionState)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState
    nVars = size(cellState, 2)

    faceDeltas = zeros(nFaces, nVars)

    # Boundary face fluxes must be set separately (faceDelta is zero at all possible boundary conditions right now)
    for f in 1:nFaces-nBdryFaces
        ownerCell = faces[f][1]
        neighbourCell = faces[f][2]

        faceDeltas[f, :] = cellState[neighbourCell, :] .- cellState[ownerCell, :]
    end

    return faceDeltas
end

# TODO: TVD Interp

######################### Convective Term Things #######################
# Returns the fractional portion of the maccormack aritificial diffusivity term (Eq. 6.58 in Anderson).
# Result must still be multiplied by (nextU - 2U + U) for each flux variable.
function macCormackAD_S(mesh, C, P, lapl_P)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells = size(cells, 1)

    S = Array{Float64, 1}(undef, nCells)
    for i in 1:nCells
        avgDimension = 1/nCells
        S[i] = C * abs(lapl_P[i]*avgDimension^2) / (lapl_P[i]*avgDimension^2 + 4*P[i])
    end

    return S
end

function unstructured_JSTEps(mesh, solutionState, k2=0.5, k4=(1/32), c4=1, gamma=1.4, R=287.05)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState

    P = cellPrimitives[:,1]
    gradP = greenGaussGrad(mesh, false, P)[1]

    sj = zeros(nCells)
    rj = zeros(nCells)
    sjCount = zeros(nCells)
    for f in 1:nFaces-nBdryFaces
        # Calculate sj, rj, eps2, eps4
        ownerCell = faces[f][1]
        neighbourCell = faces[f][2]
        d = cCenters[neighbourCell] .- cCenters[ownerCell]

        oP = P[ownerCell]
        nP = P[neighbourCell]
        farOwnerP = nP - 2*dot(d, gradP[ownerCell])
        farNeighbourP = oP + 2*dot(d, gradP[neighbourCell])
        sj[ownerCell] += (abs( nP - 2*oP + farOwnerP )/ max( abs(nP - oP) + abs(oP - farOwnerP), 0.0000000001))^2
        sjCount[ownerCell] += 1
        sj[neighbourCell] += (abs( oP - 2*nP + farNeighbourP )/ max( abs(farNeighbourP - nP) + abs(nP - oP), 0.0000000001))^2
        sjCount[neighbourCell] += 1
    end

    for c in 1:nCells
        rj[c] = mag(cellPrimitives[c,3:5]) +  sqrt(gamma * R * cellPrimitives[c,2]) # Velocity magnitude + speed of sound
        sj[c] /= sjCount[c] # Average the sj's computed by each face for each cell
    end

    rjsjF = maxInterp(mesh, rj, sj) # column one is rj, column two is sj, both at face centers

    eps2 = zeros(nFaces)
    eps4 = zeros(nFaces)
    for f in 1:nFaces-nBdryFaces
        eps2[f] = k2 * rjsjF[f,2] * rjsjF[f,1]
        eps4[f] = max(0, k4*rjsjF[f,1] - c4*eps2[f])
    end

    return eps2, eps4
end

# Requires correct cellState and cellPrimitives as input
# Classical JST, central differencing + artificial diffusion. Each face treated as 1D
function unstructured_JSTFlux(mesh, solutionState, boundaryConditions)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState
    nVars = size(cellState, 2)

    #### Boundaries Prediction ####
    for b in 1:nBoundaries
        bFunctionIndex = 2*b-1
        boundaryConditions[bFunctionIndex](mesh, solutionState, b, boundaryConditions[bFunctionIndex+1])
    end

    # Centrally differenced fluxes
    linInterp_3D(mesh, solutionState)

    #### Add JST artificial Diffusion ####
    fDeltas = faceDeltas(mesh, solutionState)
    fDGrads = greenGaussGrad_matrix(mesh, fDeltas, false)
    eps2, eps4 = unstructured_JSTEps(mesh, solutionState, 0.5, (0.5/32), 1)
    # nCells = nFaces - 1
    for f in 1:nFaces-nBdryFaces
        ownerCell = faces[f][1]
        neighbourCell = faces[f][2]
        d = cCenters[neighbourCell] .- cCenters[ownerCell]

        fD = fDeltas[f,:]
        farOwnerfD = fD .- dot(d, fDGrads[ownerCell,:,:])
        farNeighbourfD = fD .+ dot(d, fDGrads[ownerCell,:,:])

        diffusionFlux = eps2[f]*fD - eps4[f]*(farNeighbourfD - 2*fD + farOwnerfD)
        # Add diffusion flux in component form
        unitFA = normalize(fAVecs[f])
        for v in 1:nVars
            i1 = (v-1)*3+1
            i2 = i1+2
            faceFluxes[f,i1:i2] .-= (diffusionFlux[v] .* unitFA)
        end
    end

    return integrateFluxes_unstructured3D(mesh, solutionState, boundaryConditions)
end

######################### TimeStepping #######################

function decodeSolution(solutionState)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState
    nVars = size(cellState, 2)
    if nVars == 3
        return decodeSolution_1D(solutionState)
    else
        return decodeSolution_3D(solutionState)
    end
end

function decodeSolution_1D(solutionState)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState
    nCells = size(cellState, 1)
    for c in 1:nCells
        cellPrimitives[c,:] = decodePrimitives(cellState[c,1], cellState[c,2], cellState[c,3])
        # mass, xMom, eV2 x-direction fluxes
        cellFluxes[c, :] = calculateFluxes1D(cellPrimitives[c,1], cellPrimitives[c,3], cellState[c,2], cellState[c,3])
    end
end

function decodeSolution_3D(solutionState)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState
    nCells = size(cellState, 1)
    for c in 1:nCells
        cellPrimitives[c,:] = decodePrimitives3D(cellState[c,:]...)
        # mass, xMom, eV2 x,y,z-direction fluxes
        cellFluxes[c, :] = calculateFluxes3D(cellPrimitives[c,:]..., cellState[c,:]...)
    end
end

function integrateFluxes_structured1D(dx, solutionState)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState
    nCells = size(cellState, 1)

    # Recomputing flux balances, so wipe existing values
    fill!(fluxResiduals, 0)

    # Face-based integration of fluxes
    #TODO: Multiply by face area here
    for f in 2:nCells
        for v in 1:3
            flux = faceFluxes[f,v]
            # Subtract from owner cell
            fluxResiduals[f-1, v] -= flux
            # Add to neighbour cell
            fluxResiduals[f, v] += flux
        end
    end

    # Divide by cell volume
    for c in 2:nCells-1
        fluxResiduals[c,:] ./= dx[c]
    end

    fluxResiduals[1,:] = [ 0, 0, 0 ]
    fluxResiduals[nCells,:] = [ 0, 0, 0 ]
    fluxResiduals[2,:] = [ 0, 0, 0 ]
    fluxResiduals[nCells-1,:] = [ 0, 0, 0 ]

    return fluxResiduals
end

function integrateFluxes_unstructured3D(mesh, solutionState, boundaryConditions)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

    # Recomputing flux balances, so wipe existing values
    fill!(fluxResiduals, 0)
    nVars = size(fluxResiduals, 2)

    #### Boundaries Correction ####
    for b in 1:nBoundaries
        bFunctionIndex = 2*b-1
        boundaryConditions[bFunctionIndex](mesh, solutionState, b, boundaryConditions[bFunctionIndex+1])
    end

    #### Flux Integration ####
    for f in 1:nFaces
        ownerCell = faces[f][1]
        neighbourCell = faces[f][2]

        for v in 1:nVars
            i1 = (v-1)*3 + 1
            i2 = i1+2
            flow = dot(faceFluxes[f, i1:i2], fAVecs[f])

            if ownerCell > -1
                # Subtract from owner cell
                fluxResiduals[ownerCell, v] -= flow
            end
            if neighbourCell > -1
                # Add to neighbour cell
                fluxResiduals[neighbourCell, v] += flow
            end
        end
    end

    # Divide by cell volume
    for c in 1:nCells
        fluxResiduals[c,:] ./= cVols[c]
    end

    return fluxResiduals
end

######################### Boundary Conditions #######################
# Can work as a supersonic inlet if initial conditions are set to the inlet conditions
function supersonicInletBoundary(mesh, solutionState, boundaryNumber, inletConditions)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState

    Cp = 1005
    P, T, Ux, Uy, Uz = inletConditions
    rho = idealGasRho(T, P)
    xMom, yMom, zMom = [Ux, Uy, Uz] .* rho
    e = calPerfectEnergy(T, Cp)
    eV2 = rho*(e + (mag(U)^2)/2)

    boundaryFluxes = calculateFluxes3D(P, T, Ux, Uy, Uz, rho, xMom, yMom, zMom, eV2)
    currentBoundary = boundaryFaces[boundaryNumber]
    for face in currentBoundary
        faceFluxes[face,:] = boundaryFluxes
    end
end

function zeroGradientBoundary(mesh, solutionState, boundaryNumber, _)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState

    # Directly extrapolate cell center flux to boundary (zero gradient between the cell center and the boundary)
    currentBoundary = boundaryFaces[boundaryNumber]
    for face in currentBoundary
        ownerCell = max(faces[face][1], faces[face][2]) #One of these will be -1 (no cell), the other is the boundary cell we want
        faceFluxes[face, :] = cellFluxes[ownerCell, :]
    end
end

function wallBoundary(mesh, solutionState, boundaryNumber, _)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState

    currentBoundary = boundaryFaces[boundaryNumber]
    for f in currentBoundary
        ownerCell = max(faces[f][1], faces[f][2]) #One of these will be -1 (no cell), the other is the boundary cell we want

        faceP = cellPrimitives[ownerCell, 1]
        # Momentum flux is Pressure in each of the normal directions (dot product)
        faceFluxes[f, 4] = faceP
        faceFluxes[f, 8] = faceP
        faceFluxes[f, 12] = faceP

        # Mass Flux is zero
        faceFluxes[f, 1:3] .= 0.0
        # Energy Flux is zero
        faceFluxes[f, 13:15] .= 0.0
    end
end
symmetryBoundary = wallBoundary

function emptyBoundary(mesh, solutionState, boundaryNumber, _)
    return
end

# Slip wall boundary condition does not require treatment, flux is simply zero across the boundary
#TODO: Gradient calculation at slip walls

######################### Solvers #######################
function central_UnstructuredADFVM(mesh, P, T, U; initDt=0.001, endTime=0.14267, targetCFL=0.2, gamma=1.4, R=287.05, Cp=1005, Cx=0.3, debug=false, silent=true)
    if !silent
        println("Initializing solver")
    end
    ######### MESH ############
    # Extract mesh into local variables for readability
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells = size(cells, 1)
    nFaces = size(faces, 1)
    nBdryFaces = size(boundaryFaces, 1)
    bdryFaceIndices = Array(nFaces-nBdryFaces:nFaces)

    ########### Variable Arrays #############
    # State variables, values are the averages/cell center values
    rho = Array{Float64, 1}(undef, nCells)
    xMom = Array{Float64, 1}(undef, nCells)
    eV2 = Array{Float64, 1}(undef, nCells)
    stateVars = [ rho, xMom, eV2 ]
    # Flux variables for equations 2 and 3 (xMom is flux variable for the continuity eqn)
    rhoU2p = Array{Float64, 1}(undef, nCells)
    rhoUeV2PU = Array{Float64, 1}(undef, nCells)

    # Calc state and flux variables from primitives
    #TODO: Multi-D
    for i in 1:nCells
        rho[i], xMom[i], eV2[i], rhoU2p[i], rhoUeV2PU[i] = encodePrimitives3D(P[i], T[i], U[i])
    end

    if debug
        println("Rho: $rho")
        println("xMom: $xMom")
        println("eV2: $eV2")
    end

    ########### SOLVER ###########
    if !silent
        println("Starting iterations")
    end
    dt = initDt
    currTime = 0
    timeStepCounter = 0
    while currTime < endTime

        if (endTime - currTime) < dt
            dt = endTime - currTime
        end

        # Calculate fluxes through each face
        # TODO: y and z momemtum-fluxes + equations
        # xMassFlux, xMomFlux, xeV2Flux = upwindInterp(mesh, U, xMom, rhoU2p, rhoUeV2PU)
        xMassFlux, xMomFlux, xeV2Flux, faceP = linInterp(mesh, xMom, rhoU2p, rhoUeV2PU, P)
        fluxVars = [ xMassFlux, xMomFlux, xeV2Flux ]

        if debug
            println("xMass Flux: $xMassFlux")
            println("xMom Flux: $xMomFlux")
            println("xEv2 Flux: $xeV2Flux")
        end

        # Use fluxes to update values in each cell
        for i in 1:(nFaces-nBdryFaces)
            fA = mag(fAVecs[i]) #Face Area

            ownerCell = faces[i][1]
            neighbourCell = faces[i][2]
            for v in 1:3
                stateVars[v][ownerCell] -= fluxVars[v][i]*fA*dt/cVols[ownerCell]
                stateVars[v][neighbourCell] += fluxVars[v][i]*fA*dt/cVols[neighbourCell]
            end
        end

        if Cx > 0
            # Apply McCormack Artificial Diffusivity
            lapl_P, lapl_rho, lapl_xMom, lapl_eV2 = laplacian(mesh, "None", P, rho, xMom, eV2)
            S = macCormackAD_S(mesh, Cx , P, lapl_P)

            for i in 2:nCells-1
                avgD = 1/nCells
                rho[i] += S[i]*lapl_rho[i]*avgD^2
                xMom[i] += S[i]*lapl_xMom[i]*avgD^2
                eV2[i] += S[i]*lapl_eV2[i]*avgD^2
            end
        end

        if debug
            println("Rho2: $rho")
            println("xMom2: $xMom")
            println("eV22: $eV2")
        end

        # Boundaries
        ############### Boundaries ################
        # Waves never reach the boundaries, so boundary treatment doesn't need to be good
        copyValues(3, 2, stateVars)
        copyValues(2, 1, stateVars)
        copyValues(nCells-2, nCells-1, stateVars)
        copyValues(nCells-1, nCells, stateVars)

        # Decode primitive values
        for i in 1:nCells
            P[i], T[i], U[i], rhoU2p[i], rhoUeV2PU[i] = decodePrimitives3D(rho[i], xMom[i], eV2[i])
        end

        currTime += dt

        ############## CFL Calculation, timestep adjustment #############
        maxCFL = 0
        for i in 1:nCells
            dx = 1/nCells
            maxCFL = max(maxCFL, CFL(U[i], T[i], dt, dx, gamma, R))
        end
        # Adjust time step to approach target CFL
        dt *= ((targetCFL/maxCFL - 1)/5+1)

        timeStepCounter += 1
        if !silent
            @printf("Timestep: %5.0f, simTime: %8.4g, Max CFL: %8.4g \n", timeStepCounter, currTime, maxCFL)
        end

    end

    return P, U, T, rho
end

function structured1DFVM(dx::Array{Float64, 1}, cellPrimitives::Array{Float64, 2}, timeIntegrationFn=forwardEuler, fluxFunction=structured_JSTFlux1D; initDt=0.001, endTime=0.14267, targetCFL=0.2, gamma=1.4, R=287.05, Cp=1005, silent=true)
    if !silent
        println("Initializing Simulation")
    end
    nCells = size(dx, 1)
    nFaces = nCells+1
    nDims = 1
    # Each dimension adds one momentum equation
    nVars = 2+nDims
    # Each dimension adds a flux for each conserved quantity
    nFluxes = nVars*nDims

    cellState = Array{Float64, 2}(undef, nCells, nVars)
    cellFluxes = Array{Float64, 2}(undef, nCells, nFluxes)
    fluxResiduals = Array{Float64, 2}(undef, nCells, nVars)
    faceFluxes = Array{Float64, 2}(undef, nFaces, nFluxes)
    solutionState = [ cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]

    # Easy access functions
    function P(cell::Int)           return cellPrimitives[cell, 1]      end
    function T(cell::Int)           return cellPrimitives[cell, 2]      end
    function Ux(cell::Int)          return cellPrimitives[cell, 3]      end
    function rho(cell::Int)         return cellState[cell, 1]           end
    function xMom(cell::Int)        return cellState[cell, 2]           end
    function eV2(cell::Int)         return cellState[cell, 3]           end
    function massXFlux(cell::Int)   return cellFluxes[cell, 1]          end
    function xMomXFlux(cell::Int)   return cellFluxes[cell, 2]          end
    function eV2XFlux(cell::Int)    return cellFluxes[cell, 3]          end

    #Initialize States and Fluxes
    for c in 1:nCells
        # rho, xMom, total energy from P, T, U
        cellState[c, :] = encodePrimitives(P(c), T(c), Ux(c))
        # mass, xMom, eV2 x-direction fluxes
        cellFluxes[c, :] = calculateFluxes1D(P(c), Ux(c), xMom(c), eV2(c))
    end

    if !silent
        println("Starting iterations")
    end

    dt = initDt
    currTime = 0
    timeStepCounter = 0
    while currTime < endTime
        ############## Timestep adjustment #############
        maxCFL = 0.0
        for c in 1:nCells
            maxCFL = max(maxCFL, (abs(Ux(c)) + sqrt(gamma * R * T(c))) * dt / dx[c])
        end
        # Slowly approach target CFL
        dt *= ((targetCFL/maxCFL - 1)/5+1)
        # Adjust timestep to hit endtime if this is the final time step
        if (endTime - currTime) < dt
            dt = endTime - currTime
        end

        ############## Take a timestep #############
        solutionState = timeIntegrationFn(dx, fluxFunction, solutionState, dt)
        currTime += dt
        timeStepCounter += 1

        ############### Apply Boundary conditions ################
        # Waves never reach the boundaries, so boundary treatment doesn't need to be good
        allVars = [ cellState, cellFluxes ]
        copyValues(3, 2, allVars)
        copyValues(2, 1, allVars)
        copyValues(nCells-2, nCells-1, allVars)
        copyValues(nCells-1, nCells, allVars)

        if !silent
            @printf("Timestep: %5.0f, simTime: %8.4g, Max CFL: %8.4g \n", timeStepCounter, currTime, maxCFL)
        end
    end

    # P, U, T, rho
    return cellPrimitives[:,1], cellPrimitives[:,3], cellPrimitives[:,2], cellState[:,1]
end

function unstructured3DFVM(mesh, meshPath, cellPrimitives::Array{Float64, 2}, boundaryConditions, timeIntegrationFn=forwardEuler, fluxFunction=unstructured_JSTFlux; initDt=0.001, endTime=0.14267, outputInterval=0.01, targetCFL=0.2, gamma=1.4, R=287.05, Cp=1005, silent=true, restart=false, createRestartFile=true, restartFile="JuliaCFDRestart.txt")
    if !silent
        println("Initializing Simulation")
    end
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    nDims = 3

    # Each dimension adds one momentum equation
    nVars = 2+nDims
    # Each dimension adds a flux for each conserved quantity
    nFluxes = nVars*nDims

    if restart
        cellPrimitives = readRestartFile(restartFile)
    end

    # rho, xMom, total energy from P, T, Ux, Uy, Uz
    cellState = encodePrimitives3D(cellPrimitives, R, Cp)
    cellFluxes = zeros(nCells, nFluxes)
    fluxResiduals = zeros(nCells, nVars)
    faceFluxes = zeros(nFaces, nFluxes)
    solutionState = [ cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]

    # Calculates cell fluxes, primitives from cell state
    decodeSolution(solutionState)

    if !silent
        println("Starting iterations")
    end

    dt = initDt
    currTime = 0
    timeStepCounter = 0
    nextOutputTime = outputInterval
    writeOutputThisIteration = false
    vtkCounter = 1
    while currTime < endTime
        ############## Timestep adjustment #############
        maxCFL = maxCFL3D(mesh, solutionState, dt)
        # Slowly approach target CFL
        dt *= ((targetCFL/maxCFL - 1)/5+1)
        # Adjust timestep to hit endtime if this is the final time step
        if (endTime - currTime) < dt
            dt = endTime - currTime
        elseif (nextOutputTime - currTime) < dt
            dt = nextOutputTime - currTime
            writeOutputThisIteration = true
        end

        ############## Take a timestep #############
        solutionState = timeIntegrationFn(mesh, fluxFunction, solutionState, boundaryConditions, dt)
        currTime += dt
        timeStepCounter += 1

        if !silent
            @printf("Timestep: %5.0f, simTime: %8.4g, Max CFL: %8.4g \n", timeStepCounter, currTime, maxCFL)
        end

        if writeOutputThisIteration
            updateSolutionOutput(cellPrimitives, restartFile, meshPath, vtkCounter, createRestartFile)
            vtkCounter += 1
            writeOutputThisIteration = false
            nextOutputTime = nextOutputTime + outputInterval
        end
    end

    updateSolutionOutput(cellPrimitives, restartFile, meshPath, vtkCounter, createRestartFile)
end
