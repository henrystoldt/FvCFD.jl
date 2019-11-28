# finiteVolumeRoe.jl

include("constitutiveRelations.jl")
include("vectorFunctions.jl")

######################### CFL ########################
# TODO: Generalize cell size
# TODO: 3D definition of CFL: Sum up in all directions
function CFL(U, T, dt, dx, gamma=1.4, R=287.05)
    return (abs(U[1]) + sqrt(gamma * R * T)) * dt / dx
end

######################### Gradient Computation #######################
function leastSqGrad(mesh, values...)
    #TODO
    result = []
    for vals in values
    end
    return result
end

# Pass in values that have already been interpolated to faces to avoid re-interpolating
# Returns array of 3-D vectors
function greenGaussGrad(mesh, faceValues...)
    result = []
    
    # Interpret mesh
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells = size(cells, 1)
    nFaces = size(faces, 1)
    nBdryFaces = size(boundaryFaces, 1)
    bdryFaceIndices = Array(nFaces-nBdryFaces:nFaces)

    for vals in faceValues
        # Initialize gradients to zero
        grad = Array{Array{Float64, 1}, 1}(undef, nCells)
        for c in 1:nCells
            grad[c] = [0, 0, 0]
        end

        # Integrate fluxes from each face
        for f in 1:nFaces-nBdryFaces
            faceIntegral = fAVecs[f] .* vals[f]

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

function laplacian(mesh, faceValues...)
    
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

# TODO: TVD Interp

######################### Convective Term Things #######################



######################### Solvers #######################
function upwindFVMRoe1D(mesh, P, T, U; initDt=0.001, endTime=0.14267, targetCFL=0.2, gamma=1.4, R=287.05, Cp=1005, Cx=0.3, debug=false)
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
    dt = initDt
    currTime = 0
    while currTime < endTime

        if (endTime - currTime) < dt
            dt = endTime - currTime
        end
        
        # Calculate fluxes through each face
        # TODO: y and z momemtum-fluxes + equations
        # xMassFlux, xMomFlux, xeV2Flux = upwindInterp(mesh, U, xMom, rhoU2p, rhoUeV2PU)
        xMassFlux, xMomFlux, xeV2Flux, faceP = linInterp(mesh, xMom, rhoU2p, rhoUeV2PU, P)
        fluxVars = [ xMassFlux, xMomFlux, xeV2Flux ]
        
        # TODO: Update with finite volume version of artificial diffusivity
        # dx = 1/nCells
        # pCentralGrad, rhoCG, xMomCG, eV2CG = central2GradNum(dx, P, rho, xMom, eV2)
        # pDenom = central2GradDenom(dx, P)[1]
        
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
            pGrad = greenGaussGrad(mesh, faceP)[1]
            facePGrad = linInterp(mesh, pGrad)[1]
            facePGradX = Array{Float64, 1}(undef, nFaces)
            facePGradY = Array{Float64, 1}(undef, nFaces)
            facePGradZ = Array{Float64, 1}(undef, nFaces)
            for f in 1:size(facePGrad,1)
                facePGradX[f] = facePGrad[f][1]
                facePGradY[f] = facePGrad[f][2]
                facePGradZ[f] = facePGrad[f][3]
            end
            p2Grad = greenGaussGrad(mesh, facePGradX, facePGradY, facePGradZ)
            S = macCormackAD_S(mesh, [ Cx, Cx, Cx ], P, pGrad, p2Grad)
            for i in 2:nCells-1
                rho[i] += S[i]*rhoCG[i]
                xMom[i] += S[i]*xMomCG[i]
                eV2[i] += S[i]*eV2CG[i]
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

    end

    return P, U, T, rho
end

#TODO: Proper boundary treatment