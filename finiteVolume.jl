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

function structured_1D_JST_rj(T::Array{Float64, 1}, U::Array{Float64, 1}, gamma=1.4, R=287.05)
    nCells = size(T, 1)
    rj = Array{Float64, 1}(undef, nCells)

    for c in 1:nCells
        rj[c] = abs(U[c]) + sqrt(gamma * R * T[c])
    end
    return rj
end

function structured_1D_JST_Eps(dx, k2, k4, c4, P::Array{Float64, 1}, T::Array{Float64, 1}, U::Array{Float64, 1})
    nFaces = size(P, 1) + 1

    sj = structured_1D_JST_sj(P)
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

######################### Solvers #######################
function central_UnstructuredADFVM(mesh, P, T, U; initDt=0.001, endTime=0.14267, targetCFL=0.2, gamma=1.4, R=287.05, Cp=1005, Cx=0.3, debug=false, silent=true)
    if silent == false
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
    if silent == false
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
        if silent == false
            @printf("Timestep: %5.0f, simTime: %8.4g, Max CFL: %8.4g \n", timeStepCounter, currTime, maxCFL)
        end

    end

    return P, U, T, rho
end

#TODO: Runge-Kutta Timestepping
# Structured: relies on cells being ordered sequentially, uses FDM mesh definition
function JST_Structured1DFVM(dx, P, T, U; initDt=0.001, endTime=0.14267, targetCFL=0.2, gamma=1.4, R=287.05, Cp=1005, silent=true)
    nCells = size(dx, 1)

    rho = Array{Float64, 1}(undef, nCells)
    xMom = Array{Float64, 1}(undef, nCells)
    eV2 = Array{Float64, 1}(undef, nCells)
    rhoU2p = Array{Float64, 1}(undef, nCells)
    rhoUeV2PU = Array{Float64, 1}(undef, nCells)

    for i in 1:nCells
        rho[i], xMom[i], eV2[i] = encodePrimitives(P[i], T[i], U[i])
        rhoU2p[i] = xMom[i]*U[i] + P[i]
        rhoUeV2PU[i] = U[i]*eV2[i] + P[i]*U[i]
    end

    stateVars = [ rho, xMom, eV2 ]

    CFL = Array{Float64, 1}(undef, nCells)

    if silent == false
        println("Starting iterations")
    end
    dt = initDt
    currTime = 0
    timeStepCounter = 0
    while currTime < endTime
        if (endTime - currTime) < dt
            dt = endTime - currTime
        end
        
        ############## Predictor #############
        xMassFlux, xMomFlux, xeV2Flux = structured_1DlinInterp(dx, xMom, rhoU2p, rhoUeV2PU)
        xMassDelta, xMomDelta, xeV2Delta = structured_1DFaceDelta(dx, rho, xMom, eV2)
        eps2, eps4 = structured_1D_JST_Eps(dx, 1, (1/32), 4, P, T, U)

        fluxVars = [ xMassFlux, xMomFlux, xeV2Flux ]
        deltas = [ xMassDelta, xMomDelta, xeV2Delta ]

        # nCells = nFaces - 1
        for f in 2:(nCells)
            for v in 1:3
                diffusionFlux = eps2[f]*deltas[v][f] - eps4[f]*(deltas[v][f+1] - 2*deltas[v][f] + deltas[v][f-1])
                flux = (fluxVars[v][f] - diffusionFlux)*dt
                stateVars[v][f-1] -= flux/dx[f-1]
                stateVars[v][f] += flux/dx[f]
            end
        end

        # Decode
        for c in 1:nCells
            P[c], T[c], U[c] = decodePrimitives(rho[c], xMom[c], eV2[c])
            rhoU2p[c] = xMom[c]*U[c] + P[c]
            rhoUeV2PU[c] = U[c]*eV2[c] + P[c]*U[c]
        end

        ############### Boundaries ################
        # Waves never reach the boundaries, so boundary treatment doesn't need to be good
        allVars = [ rho, rhoU2p, rhoUeV2PU, xMom, eV2, U, P, T ]
        copyValues(3, 2, allVars)
        copyValues(2, 1, allVars)
        copyValues(nCells-2, nCells-1, allVars)
        copyValues(nCells-1, nCells, allVars)
        
        ############## CFL Calculation, timestep adjustment #############
        for i in 1:nCells
            CFL[i] = (abs(U[i]) + sqrt(gamma * R * T[i])) * dt / dx[i]
        end
        maxCFL = maximum(CFL)

        # Adjust time step to slowly approach target CFL
        dt *= ((targetCFL/maxCFL - 1)/5+1)

        currTime += dt

        timeStepCounter += 1
        if silent == false
            @printf("Timestep: %5.0f, simTime: %8.4g, Max CFL: %8.4g \n", timeStepCounter, currTime, maxCFL)
        end
    end

    return P, U, T, rho
end

#TODO: Proper boundary treatment