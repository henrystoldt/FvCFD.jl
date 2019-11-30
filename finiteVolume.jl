using Printf
include("constitutiveRelations.jl")
include("vectorFunctions.jl")

######################### CFL ########################
# TODO: Generalize cell size
# TODO: 3D definition of CFL: Sum up in all directions
function CFL(U, T, dt, dx, gamma=1.4, R=287.05)
    return (abs(U[1]) + sqrt(gamma * R * T)) * dt / dx
end

######################### Gradient Computation #######################
#TODO: make all these function return a single array if you pass in a single value
function leastSqGrad(mesh, values...)
    #TODO
    result = []
    for vals in values
    end
    return result
end

# Pass in values that have already been interpolated to faces to avoid re-interpolating
# Returns array of 3-D vectors
function greenGaussGrad(mesh, valuesAtFaces=false, values...)
    result = []

    # Interpret mesh
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells = size(cells, 1)
    nFaces = size(faces, 1)
    nBdryFaces = size(boundaryFaces, 1)
    bdryFaceIndices = Array(nFaces-nBdryFaces:nFaces)

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
    nCells = size(cells, 1)
    nFaces = size(faces, 1)
    nBdryFaces = size(boundaryFaces, 1)
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
    nCells = size(cells, 1)
    nFaces = size(faces, 1)
    nBdryFaces = size(boundaryFaces, 1)
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
    nFaces = size(faces, 1)
    nBdryFaces = size(boundaryFaces, 1)

    fVals = Array{Float64, 1}(undef, nFaces)

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

function structured_1DlinInterp(dx, values...)
    result = []
    nFaces = size(dx, 1) + 1

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

        push!(fVals, defaultVal)
        ########### Do Interpolations ###########
        for i in 2:nFaces-1
            c1Dist = dx[i-1]/2
            c2Dist = dx[i]/2
            totalDist = c1Dist + c2Dist
            push!(fVals, vals[i-1].*(c2Dist/totalDist) .+ vals[i].*(c1Dist/totalDist))
        end
        push!(fVals, defaultVal)

        push!(result, fVals)
    end
    return result
end

function structured_1DMaxInterp(dx, values...)
    result = []
    nFaces = size(dx, 1) + 1

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

        push!(fVals, defaultVal)
        ########### Do Interpolations ###########
        for i in 2:nFaces-1
            push!(fVals, max(vals[i-1], vals[i]))
        end
        push!(fVals, defaultVal)

        push!(result, fVals)
    end
    return result
end

function structured_1DFaceDelta(dx, values...)
    result = []
    nFaces = size(dx, 1) + 1

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

        push!(fVals, defaultVal)
        ########### Do Interpolations ###########
        for f in 2:nFaces-1
            push!(fVals, vals[f] - vals[f-1])
        end
        push!(fVals, defaultVal)

        push!(result, fVals)
    end
    return result
end

# Matrix versions
# results are stored in faceValues vectors
function checkInputs_structured1DInterp(dx, faceValues::Array{Float64,2}, cellValues::Array{Float64, 2})
    nCells = size(cellValues, 1)
    nFaces = size(faceValues, 1)
    nVars = size(faceValues, 2)
    nVars2 = size(cellValues, 2)
    if nCells != nFaces-1 || nVars2 != nVars
        throw(ArgumentError("ERROR: nCells ($nCells) != nFaces ($nFaces) or Number of variables stored in cells($nVars2) and faces($nVars) are not equal!"))
    end
    nDxs = size(dx, 1)
    if nDxs != nCells
        throw(ArgumentError("ERROR: nCells ($nCells) != nDxs ($nDxs)"))
    end
end

function structured_1DlinInterp(dx, faceValues::Array{Float64,2}, cellValues::Array{Float64, 2})
    nFaces = size(faceValues, 1)
    nVars = size(faceValues, 2)
    checkInputs_structured1DInterp(dx, faceValues, cellValues)

    # Do interpolation
    for f in 2:nFaces-1
        for v in 1:nVars
            c1Dist = dx[f-1]/2
            c2Dist = dx[f]/2
            totalDist = c1Dist + c2Dist
            faceValues[f, v] = cellValues[f-1, v].*(c2Dist/totalDist) .+ cellValues[f, v].*(c1Dist/totalDist)
        end
    end
end

function structured_1DMaxInterp(dx, faceValues::Array{Float64,2}, cellValues::Array{Float64, 2})
    nFaces = size(faceValues, 1)
    nVars = size(faceValues, 2)
    checkInputs_structured1DInterp(dx, faceValues, cellValues)

    # Do interpolation
    for f in 2:nFaces-1
        for v in 1:nVars
            faceValues[f, v] = max(cellValues[f-1,v],  cellValues[f,v])
        end
    end
end

function structured_1DFaceDelta(dx, faceValues::Array{Float64,2}, cellValues::Array{Float64, 2})
    nFaces = size(faceValues, 1)
    nVars = size(faceValues, 2)
    checkInputs_structured1DInterp(dx, faceValues, cellValues)

    # Do interpolation
    for f in 2:nFaces-1
        for v in 1:nVars
            faceValues[f, v] = cellValues[f, v] - cellValues[f-1, v]
        end
    end
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

function structured_1D_JST_sj(P::Array{Float64, 1})
    nCells = size(P, 1)
    sj = Array{Float64, 1}(undef, nCells)

    # Boundary values unset
    sj[1] = 0.0
    sj[nCells] = 0.0

    for c in 2:nCells-1
        sj[c] = abs( (P[c+1] - 2*P[c] + P[c-1]) / (P[c+1] + 2*P[c] + P[c-1]) )
    end

    return sj
end

# Pass in pressure or density for the pressure-based or density-based version respectively
# Density-based version detects contact discontinuities
function structured_1D_JST_sj2(dvar)
    nCells = size(dvar, 1)
    sj = Array{Float64, 1}(undef, nCells)

    # Boundary values unset
    sj[1] = 0.0
    sj[nCells] = 0.0

    for c in 2:nCells-1
        sj[c] = (abs( dvar[c] - dvar[c-1] )/ max( abs(dvar[c]) + abs(dvar[c-1]), 0.0000000001))^2
    end

    return sj
end

function structured_1D_JST_rj(T::Array{Float64, 1}, U::Array{Float64, 1}, gamma=1.4, R=287.05)
    nCells = size(T, 1)
    rj = Array{Float64, 1}(undef, nCells)

    for c in 1:nCells
        rj[c] = abs(U[c]) + sqrt(gamma * R * T[c])
    end
    return rj
end

function structured_1D_JST_Eps(dx, k2, k4, c4, P::Array{Float64, 1}, T::Array{Float64, 1}, U::Array{Float64, 1}, rho::Array{Float64, 1})
    nFaces = size(P, 1) + 1

    # sj = structured_1D_JST_sj(P)

    # Pressure sensor seems to work best for the shock tube
    dP = structured_1DFaceDelta(dx, P)[1]
    sj = structured_1D_JST_sj2(dP)

    # drho = structured_1DFaceDelta(dx, rho)[1]
    # sj = structured_1D_JST_sj2(drho)

    #TODO: Pass through gamma and R
    rj = structured_1D_JST_rj(T, U)
    sjF, rjF = structured_1DMaxInterp(dx, sj, rj)

    eps2 = Array{Float64, 1}(undef, nFaces)
    eps4 = Array{Float64, 1}(undef, nFaces)
    for f in 2:nFaces-1
        eps2[f] = k2 * sjF[f] * rjF[f]
        eps4[f] = max(0, k4*rjF[f] - c4*eps2[f])
    end

    return eps2, eps4
end

# Pass in pressure or density for the pressure-based or density-based version respectively
# Density-based version detects contact discontinuities
function structured_1D_JST_sj2(dvar::Array{Float64, 1}, sj::Array{Float64, 1})
    nCells = size(dvar, 1)

    for c in 2:nCells-1
        sj[c] = (abs( dvar[c] - dvar[c-1] )/ max( abs(dvar[c]) + abs(dvar[c-1]), 0.0000000001))^2
    end
end

function structured_1D_JST_rj(cellPrimitives::Array{Float64, 2}, rj::Array{Float64, 1}, gamma=1.4, R=287.05)
    nCells = size(T, 1)

    for c in 1:nCells
        rj[c] = abs(U[c]) + sqrt(gamma * R * T[c])
    end
end

function structured_1D_JST_Eps(dx, k2, k4, c4, cellPrimitives::Array{Float64,2}, gamma=1.4, R=287.05)
    nFaces = size(cellPrimitives, 1) + 1
    # Pressure sensor seems to work best for the shock tube
    dP = structured_1DFaceDelta(dx, cellPrimitives[:,1])[1]
    sj = structured_1D_JST_sj2(dP)

    rj = structured_1D_JST_rj(cellPrimitives[:,2], cellPrimitives[:,3], gamma, R)
    sjF, rjF = structured_1DMaxInterp(dx, sj, rj)

    eps2 = Array{Float64, 1}(undef, nFaces)
    eps4 = Array{Float64, 1}(undef, nFaces)
    for f in 2:nFaces-1
        eps2[f] = k2 * sjF[f] * rjF[f]
        eps4[f] = max(0, k4*rjF[f] - c4*eps2[f])
    end

    return eps2, eps4
end

function structured_JSTFlux(dx, solutionState)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState
    nCells = size(cellState, 1)
    nFaces = nCells + 1
    faceDeltas = Array{Float64, 2}(undef, nFaces, 3)

    # Centrally differenced fluxes
    structured_1DlinInterp(dx, faceFluxes, cellFluxes)

    #### Add JST artificial Diffusion ####
    structured_1DFaceDelta(dx, faceDeltas, cellState)
    eps2, eps4 = structured_1D_JST_Eps(dx, 0.4, (1/32), 4, cellPrimitives)
    # nCells = nFaces - 1
    for f in 2:(nCells)
        for v in 1:3
            diffusionFlux = eps2[f]*faceDeltas[f, v] - eps4[f]*(faceDeltas[f+1, v] - 2*faceDeltas[f, v] + faceDeltas[f-1, v])
            faceFluxes[f,v] -= diffusionFlux
        end
    end

    return integrateFluxes_structured(dx, solutionState)
end

######################### TimeStepping #######################
function integrateFluxes_structured(dx, solutionState)
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

function forwardEuler(mesh, fluxResidualFn, solutionState, dt)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState

    fluxResiduals = fluxResidualFn(mesh, solutionState)
    cellState .+= fluxResiduals.*dt

    return solutionState
end

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

#TODO: Runge-Kutta Timestepping
# Structured: relies on cells being ordered sequentially
# Expected and used Data Structures defined in dataStructureDefinitions.md
function JST_Structured1DFVM(dx::Array{Float64, 1}, cellPrimitives::Array{Float64, 2}; initDt=0.001, endTime=0.14267, targetCFL=0.2, gamma=1.4, R=287.05, Cp=1005, silent=true)
    if !silent
        println("Initializing Simulation")
    end
    nCells = size(dx, 1)
    nFaces = nCells+1
    nDims = 1
    # Each dimension adds one momentume equation
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

        # Mass Flux
        cellFluxes[c, 1] = cellState[c, 2]
        # x-Momentum Flux
        cellFluxes[c, 2] = xMom(c)*Ux(c) + P(c)
        # x-Total Energy Flux
        cellFluxes[c, 3] = Ux(c)*eV2(c) + P(c)*Ux(c)
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
        if isnan(maxCFL)
            println(cellPrimitives[:,2])
        end
        # Adjust timestep to hit endtime if this is the final time step
        if (endTime - currTime) < dt
            dt = endTime - currTime
        end


        ############## Take a timestep #############
        forwardEuler(dx, structured_JSTFlux, solutionState, dt)
        currTime += dt
        timeStepCounter += 1

        ############## Decode Primitives #############
        for c in 1:nCells
            cellPrimitives[c,:] = decodePrimitives(rho(c), xMom(c), eV2(c))
            # Mass Flux
            cellFluxes[c, 1] = cellState[c, 2]
            # x-Momentum Flux
            cellFluxes[c, 2] = xMom(c)*Ux(c) + P(c)
            # x-Total Energy Flux
            cellFluxes[c, 3] = Ux(c)*eV2(c) + P(c)*Ux(c)
        end

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

    return cellPrimitives[:,1], cellPrimitives[:,3], cellPrimitives[:,2], cellState[:,1]
end

#TODO: Proper boundary treatment
