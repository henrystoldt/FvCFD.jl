using Plots
using Plots.PlotMeasures
using LaTeXStrings
plotly()

# TODO: Add cell center and face center coordinates
# TODO: Make mesh global?
################### Mesh Format Definition ####################
# Cells and faces are numbered according to their storage location in the mesh arrays
# Faces are numbered such that boundary faces come last
# Face area vector point outward from the owner cells, into the neighbour cell
# Cells must be convex, composed of planar faces
# Mesh is defined as a face-based tuple containing:
# (
#     cells (list of lists, each containing the indices of the faces that make up the cell)
#     faces (list of lists, each sublist containing two cell indices: the owner cell and the neighbour cell
#     fAVectors (list of face area vectors)
#     fCenters (list of vectors)
#     boundaryFaces (list of lists, each containing the indices of cells on the ith boundary)
#     cellVolumes (list of cell volumes)
#     cellCenters (list of vectors)
# )

################# Nomenclature #################
# All units SI
# e = internal energy
# P = Pressure
# rho = Density
# T = Temperature
# U = Velocity
# xMom = x-Momentum = rho*U
# eV2 = total energy = rho*(e + U^2/2)
# rhoU2p = flux of x-momentum = rho*U^2 + P
# rhoUeV2PU = flux of energy = U*eV2 + P*U

######################### Utility Functions ########################
# Assumes sizes are equal - use it properly!
function dot(vec1, vec2)
    s1 = size(vec1, 1)
    s2 = size(vec2, 1)
    if s1 != s2
        return 0
    end
    sum = 0
    for i in 1:size(vec1, 1)
        sum += vec1[i]*vec2[i]
    end
    return sum
end

function mag(vec)
    sqrSum = 0
    for i = 1:size(vec,1)
        sqrSum += vec[i]*vec[i]
    end
    return sqrt(sqrSum)
end

# TODO: Generalize cell size
# TODO: 3D definition of CFL: Sum up in all directions
function CFL(U, T, dt, dx, gamma=1.4, R=287.05)
    return (abs(U[1]) + sqrt(gamma * R * T)) * dt / dx
end

######################### Constitutive Relations #######################
function idealGasRho(T, P, R=287.05)
    return P/(R*T)
end

function idealGasP(rho, T, R=287.05)
    # PV = mRT
    # P = rho R T
    return rho*R*T
end

function calPerfectEnergy(T, Cp=1005, R=287.05)
    return T*(Cp-R)
end

function calPerfectT(e, Cp=1005, R=287.05)
    return e/(Cp-R)
end

function decodePrimitives(rho, xMom, eV2, R=287.05, Cp=1005)
    U = xMom / rho
    e = (eV2/rho) - U*U/2
    T = calPerfectT(e, Cp)
    P = idealGasP(rho, T, R)
    return P, T, U
end

#TODO: Multi-D
function decodePrimitives3D(rho::Float64, xMom::Float64, eV2::Float64, R=287.05, Cp=1005)
    U = [ xMom/rho, 0.0, 0.0 ]
    e = (eV2/rho) - (mag(U)^2)/2
    T = calPerfectT(e, Cp)
    P = idealGasP(rho, T, R)
    
    Ux = U[1]
    
    # x-Fluxes
    rhoU2p = xMom*Ux + P
    rhoUeV2PU = Ux*eV2 + P.*Ux

    return P, T, U, rhoU2p, rhoUeV2PU
end

# Returns rho, xMom, and eV2
function encodePrimitives(P, T, U, R=287.05, Cp=1005)
    rho = idealGasRho(T, P)
    xMom = U*rho
    e = calPerfectEnergy(T)
    eV2 = rho*(e + U*U/2)
    return rho, xMom, eV2
end

#TODO: Multi-D
function encodePrimitives3D(P::Float64, T::Float64, U::Array{Float64, 1}, R=287.05, Cp=1005)
    # State Variables
    rho = idealGasRho(T, P)
    xMom = U[1]*rho
    e = calPerfectEnergy(T)
    eV2 = rho*(e + (mag(U)^2)/2)

    Ux = U[1]
    
    # x-Fluxes
    rhoU2p = xMom*Ux + P
    rhoUeV2PU = Ux*eV2 + P*Ux

    return rho, xMom, eV2, rhoU2p, rhoUeV2PU
end

######################### Gradient Computation #######################
function forwardGradient(dx, values...)
    result = []
    for vals in values
        n = size(vals,1)
        grad = Array{Float64, 1}(undef, n)
        grad[n] = 0
        for i in 1:n-1
            grad[i] = (vals[i+1] - vals[i]) / dx[i]
        end
        push!(result, grad)
    end
    return result
end

function backwardGradient(dx, values...)
    result = []
    for vals in values
        n = size(vals,1)
        grad = Array{Float64, 1}(undef, n)
        grad[1] = 0
        for i in 2:n
            grad[i] = (vals[i] - vals[i-1])/dx[i-1]
        end
        push!(result, grad)
    end
    return result
end

# Unused
function centralGradient(dx, values...)
    result = []
    for vals in values
        n = size(vals,1)
        grad = Array{Float64, 1}(undef, n)
        grad[1] = 0
        grad[n] = 0
        for i in 2:n-1
            grad[i] = (vals[i+1] - vals[i-1])/(dx[i] + dx[i-1])
        end
        push!(result, grad)
    end
    return result
end

# Unused
function upwindGradient(dx, U, values...)
    result = []
    for vals in values
        n = size(vals,1)
        grad = Array{Float64, 1}(undef, n)

        if U[1] > 0
            grad[1] = 0
        end
        if U[n] <= 0
            grad[n] = 0
        end
        
        for i in 1:n
            if U[i] > 0 && i > 1
                grad[i] = (vals[i] - vals[i-1])/dx[i-1]
            elseif U[i] == 0 && i > 1 && i < n
                grad[i] = (vals[i+1] - vals[i-1]) / (dx[i-1] + dx[i])
            elseif i < n
                grad[i] = (vals[i+1] - vals[i])/dx[i]
            end
        end
        push!(result, grad)
    end
    return result
end

# Just the numerator of the central second derivative, used for artificial diffusion term
function central2GradNum(dx, values...)
    result = []
    for vals in values
        n = size(vals,1)
        grad = Array{Float64, 1}(undef, n)
        grad[1] = 0
        grad[n] = 0
        for i in 2:n-1
            grad[i] = vals[i+1] - 2*vals[i] + vals[i-1]
        end
        push!(result, grad)
    end
    return result
end

# Denominator of artificial diffusion term
function central2GradDenom(dx, values...)
    result = []
    for vals in values
        n = size(vals,1)
        grad = Array{Float64, 1}(undef, n)
        grad[1] = 0
        grad[n] = 0
        for i in 2:n-1
            grad[i] = vals[i+1] + 2*vals[i] + vals[i-1]
        end
        push!(result, grad)
    end
    return result
end

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
    cells, faces, fAVecs, bdryFaces, cellVols = mesh
    nCells = size(cells, 1)
    nFaces = size(faces, 1)
    nBdryFaces = size(bdryFaces, 1)
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
            grad[c] = grad[c] / cellVols[c]
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

####################### Face value interpolation ####################
# Interpolates to all INTERIOR faces
function upwindInterp(mesh, U, values...)
    faces = mesh[2]
    fAVecs = mesh[3]
    bdryFaces = mesh[4]
    nFaces = size(faces, 1)
    nBdryFaces = size(bdryFaces, 1)
    
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

# TODO: Take into account distances to each cell center
function linInterp(mesh, values...)    
    faces = mesh[2]
    bdryFaces = mesh[4]
    nFaces = size(faces, 1)
    nBdryFaces = size(bdryFaces, 1)
    
    fVals = Array{Float64, 1}(undef, nFaces)

    result = []
    for vals in values
        fVals = []

        defaultVal = []
        if size(vals[1], 1) ==1
            defaultVal = 0
        else
            for a in vals[1]
                push!(defaultVal, 0.0)
            end
        end

        for i in 1:nFaces
            if i > nFaces-nBdryFaces
                # Faces on boundaries not treated, must be set externally
                push!(fVals, defaultVal)
                
            else
                # Find velocity at face using linear interpolation
                c1 = faces[i][1]
                c2 = faces[i][2]
                push!(fVals, vals[c1].*0.5 .+ vals[c2].*0.5)
            end
        end

        push!(result, fVals)
    end
    return result
end

# TODO: TVD Interp

######################### Boundary Conditions #######################
function copyValues(fromIndex, toIndex, varArrays)
    for varArray in varArrays
        varArray[toIndex] = varArray[fromIndex]
    end
end

function setValues(value, indices, varArrays)
    for varArray in varArrays
        for i in indices
            if i <= size(varArray, 1)
                varArray[i] = value
            else
                println("Error, index $i out of bounds for array: $varArray")
            end
        end
    end
end

######################### Initialization #######################
function checkPRatio(Pratio)
    if Pratio < 0
        throw(ArgumentError("Pratio must be positive"))
    end
    if Pratio > 1
        Pratio = 1/Pratio
    end

    return Pratio
end

function initializeShockTubeFDM(nCells=100; domainLength=1, Pratio=10)
    Pratio = checkPRatio(Pratio)
    
    # Create arrays to store data (cell # = position in array)
    dx = Array{Float64, 1}(undef, nCells)
    U = Array{Float64, 1}(undef, nCells)
    P = Array{Float64, 1}(undef, nCells)
    T = Array{Float64, 1}(undef, nCells)

    # Apply initial conditions (Fig. 1 in Henry's paper)
    for i in 1:nCells
        if i <= (nCells/2)
            T[i] = 0.00348432
            P[i] = 1
        else
            T[i] = 0.00278746
            P[i] = Pratio
        end
        U[i] = 0

        # dx = Distance between cell centers i and i+1
        dx[i] = domainLength / nCells
    end

    return dx, P, T, U
end

# Wrapper for FDM initialization function, adding a mesh definition suitable for FVM and vector-format velocity
function initializeShockTubeFVM(nCells=100; domainLength=1, Pratio=10)
    dx, P, T, U = initializeShockTubeFDM(nCells, domainLength=domainLength, Pratio=Pratio)
    U = []

    #Shock tube dimensions
    h = 0.1
    w = 0.1

    cells = []
    faces = []
    fAVecs = []
    cVols = []
    boundaryFaces = [ [nCells,], [nCells+1,] ]

    fAVec = [h*w, 0, 0]
    cV = h*w*domainLength/nCells
    for i in 1:nCells      
        if i == 1
            push!(cells, [nCells, 1])    
            push!(faces, [i, i+1]) 
        elseif i == nCells
            push!(cells, [nCells-1, nCells+1])
        else
            push!(cells, [i-1, i])
            push!(faces, [i, i+1]) 
        end

        push!(U, [0.0, 0.0, 0.0] )
        push!(cVols, cV)
        push!(fAVecs, fAVec)
    end

    # Last face
    push!(fAVecs, fAVec)
    # Boundary faces
    push!(faces, [-1,1])
    push!(faces, [nCells, -1])

    # Returns in mesh format
    mesh = [ cells, faces, fAVecs, boundaryFaces, cVols ]
    return mesh, P, T, U
end

############################ Plotting ############################
function plotShockTubeResults(P, U, T, rho)
    plots = []
    xAxis = Array{Float64, 1}(undef, nCells)
    for i in 1:nCells
        xAxis[i] = i/nCells - 1/(2*nCells)
    end

    pPlot = plot(xAxis, P, label="P (Pa)", title="Pressure", xlabel="x (m)", bottom_margin=15mm, left_margin=10mm)
    rhoPlot = plot(xAxis, rho, label="rho (kg/m3)", title="Density", xlabel="x (m)", bottom_margin=15mm, left_margin=10mm)
    uPlot = plot(xAxis, U, label="Velocity (m/s)", title="Velocity", xlabel="x (m)", bottom_margin=15mm, left_margin=10mm)
    TPlot = plot(xAxis, T, label="T (K)", title="Temperature", xlabel="x (m)", bottom_margin=15mm, left_margin=10mm)
    plots = [pPlot, rhoPlot, uPlot, TPlot]
    plot(plots..., layout=(2, 2), size=(1720, 880), window_title="Euler1D_Draft_Henry", legend=false)
    gui()
end

######################### Solvers #######################
function updateCellVal(initVal, flux, faceArea, dt, cellVolume)
    return (initVal*cellVolume + flux*faceArea*dt) / cellVolume
end

# Pass in initial values for each variable
# Shock Tube (undisturbed zero gradient) boundary conditions assumed
#TODO: Issue with shock position
function macCormack1DFDM(dx, P, T, U; initDt=0.001, endTime=0.14267, targetCFL=0.2, gamma=1.4, R=287.05, Cp=1005, Cx=0.3)
    nCells = size(dx, 1)

    rho = Array{Float64, 1}(undef, nCells)
    e = Array{Float64, 1}(undef, nCells)

    for i in 1:nCells
        rho[i] = idealGasRho(T[i], P[i])
        e[i] = calPerfectEnergy(T[i])
    end

    drhoPred = Array{Float64, 1}(undef, nCells)
    duPred = Array{Float64, 1}(undef, nCells)
    dePred = Array{Float64, 1}(undef, nCells)

    rhoPred = Array{Float64, 1}(undef, nCells)
    TPred = Array{Float64, 1}(undef, nCells)
    UPred = Array{Float64, 1}(undef, nCells)
    PPred = Array{Float64, 1}(undef, nCells)
    ePred = Array{Float64, 1}(undef, nCells)

    CFL = Array{Float64, 1}(undef, nCells)

    dt = initDt
    currTime = 0
    while currTime < endTime
        if (endTime - currTime) < dt
            dt = endTime - currTime
        end

        ############## Predictor #############
        drhodx, dudx, dpdx, dedx = backwardGradient(dx, rho, U, P, e)
        pCentralGrad, rhoCG, uCG, eCG = central2GradNum(dx, P, rho, U, e)
        pDenom = central2GradDenom(dx, P)[1]

        for i in 2:(nCells-1)
            # Eq. 6.1, 6.2, 6.4
            drhoPred[i] = -(rho[i]*dudx[i] + U[i]*drhodx[i])
            duPred[i] = -(U[i]*dudx[i] + dpdx[i]/rho[i])
            dePred[i] = -(U[i]*dedx[i] + P[i]*dudx[i]/rho[i])

            S = Cx * abs(pCentralGrad[i])/pDenom[i]
            rhoPred[i] = rho[i] + drhoPred[i]*dt + S*rhoCG[i]
            UPred[i] = U[i] + duPred[i]*dt + S*uCG[i]
            ePred[i] = e[i] + dePred[i]*dt + S*eCG[i]
            TPred[i] = calPerfectT(ePred[i])
            PPred[i] = idealGasP(rhoPred[i], TPred[i])
        end

        ############### Corrector ################
        drhodx, dudx, dpdx, dedx = forwardGradient(dx, rhoPred, UPred, PPred, ePred)
        pCentralGrad, rhoCG, uCG, eCG = central2GradNum(dx, PPred, rhoPred, UPred, ePred)
        pDenom = central2GradDenom(dx, PPred)[1]

        for i in 2:(nCells-1)
            drhoPred2 = -(rhoPred[i]*dudx[i] + UPred[i]*drhodx[i])
            duPred2 = -(UPred[i]*dudx[i] + dpdx[i]/rhoPred[i])
            dePred2 = -(UPred[i]*dedx[i] + PPred[i]*dudx[i]/rhoPred[i])

            # Perform timestep using average gradients
            S = Cx * abs(pCentralGrad[i])/pDenom[i]
            rho[i] += (drhoPred2 + drhoPred[i])*dt/2 + S*rhoCG[i]
            U[i] += (duPred2 + duPred[i])*dt/2 + S*uCG[i]
            e[i] += (dePred2 + dePred[i])*dt/2 + S*eCG[i]
            T[i] = calPerfectT(e[i])
            P[i] = idealGasP(rho[i], T[i])
        end

        ############### Boundaries ################
        # Waves never reach the boundaries, so boundary treatemnt doesn't need to be good
        allVars = [ rho, U, e, T, P ]
        copyValues(3, 2, allVars)
        copyValues(2, 1, allVars)
        copyValues(nCells-2, nCells-1, allVars)
        copyValues(nCells, nCells, allVars)
        
        ############## CFL Calculation, timestep adjustment #############
        for i in 1:nCells
            CFL[i] = (abs(U[i]) + sqrt(gamma * R * T[i])) * dt / dx[i]
        end
        maxCFL = maximum(CFL)

        # Adjust time step to slowly approach target CFL
        dt *= ((targetCFL/maxCFL - 1)/5+1)
        
        currTime += dt
    end

    return P, U, T, rho
end

function macCormack1DConservativeFDM(dx, P, T, U; initDt=0.001, endTime=0.14267, targetCFL=0.2, gamma=1.4, R=287.05, Cp=1005, Cx=0.3)
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
    
    drhoP = Array{Float64, 1}(undef, nCells)
    dxMP = Array{Float64, 1}(undef, nCells)
    deV2P = Array{Float64, 1}(undef, nCells)

    xMomP = Array{Float64, 1}(undef, nCells)
    eV2P = Array{Float64, 1}(undef, nCells)
    rhoP = Array{Float64, 1}(undef, nCells)
    PP = Array{Float64, 1}(undef, nCells)
    rhoU2pP = Array{Float64, 1}(undef, nCells)
    rhoUeV2PUP = Array{Float64, 1}(undef, nCells)

    CFL = Array{Float64, 1}(undef, nCells)

    dt = initDt
    currTime = 0
    while currTime < endTime
        if (endTime - currTime) < dt
            dt = endTime - currTime
        end
        
        ############## Predictor #############
        dxMomdx, drhoU2pdx, drhoUeV2PU = forwardGradient(dx, xMom, rhoU2p, rhoUeV2PU)
        pCentralGrad, rhoCG, xMomCG, eV2CG = central2GradNum(dx, P, rho, xMom, eV2)
        pDenom = central2GradDenom(dx, P)[1]

        for i in 2:(nCells-1)
            # Eq. 2.99, 2.105, 2.106
            drhoP[i] = -dxMomdx[i]
            dxMP[i] = -drhoU2pdx[i]
            deV2P[i] = -drhoUeV2PU[i]

            # Predict
            S = Cx * abs(pCentralGrad[i]) / pDenom[i]
            rhoP[i] = rho[i] + drhoP[i]*dt + S*rhoCG[i]
            xMomP[i] = xMom[i] + dxMP[i]*dt + S*xMomCG[i]
            eV2P[i] = eV2[i] + deV2P[i]*dt + S*eV2CG[i]

            # Decode
            PP[i], TP, UP = decodePrimitives(rhoP[i], xMomP[i], eV2P[i])
            rhoU2pP[i] = xMomP[i]*UP + PP[i]
            rhoUeV2PUP[i] = UP*eV2P[i] + PP[i]*UP
        end

        ############### Corrector ################
        # Rearward differences to compute gradients
        dxMomdxP, drhoU2pdxP, drhoUeV2PUP = backwardGradient(dx, xMomP, rhoU2pP, rhoUeV2PUP)
        pCentralGradP, rhoCGP, xMomCGP, eV2CGP = central2GradNum(dx, PP, rhoP, xMomP, eV2P)
        pDenomP = central2GradDenom(dx, PP)[1]

        for i in 2:(nCells-1)
            drhoP2 = -dxMomdxP[i]
            dxMP2 = -drhoU2pdxP[i]
            deV2P2 = -drhoUeV2PUP[i]

            # Perform timestep using average gradients
            S = Cx * abs(pCentralGradP[i]) / pDenomP[i]
            rho[i] += (drhoP2 + drhoP[i])*dt/2 + S*rhoCGP[i]
            xMom[i] += (dxMP2 + dxMP[i])*dt/2 + S*xMomCGP[i]
            eV2[i] += (deV2P2 + deV2P[i])*dt/2 + S*eV2CGP[i]

            # Decode
            P[i], T[i], U[i] = decodePrimitives(rho[i], xMom[i], eV2[i])
            rhoU2p[i] = xMom[i]*U[i] + P[i]
            rhoUeV2PU[i] = U[i]*eV2[i] + P[i]*U[i]
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
    end

    return P, U, T, rho
end

function upwind1DConservativeFDM(dx, P, T, U; initDt=0.001, endTime=0.14267, targetCFL=0.1, gamma=1.4, R=287.05, Cp=1005, Cx=0.3)
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

    CFL = Array{Float64, 1}(undef, nCells)

    dt = initDt
    currTime = 0
    while currTime < endTime
        if (endTime - currTime) < dt
            dt = endTime - currTime
        end
        
        ############## Predictor #############
        dxMomdx, drhoU2pdx, drhoUeV2PU = upwindGradient(dx, U, xMom, rhoU2p, rhoUeV2PU)
        pCentralGrad, rhoCG, xMomCG, eV2CG = central2GradNum(dx, P, rho, xMom, eV2)
        pDenom = central2GradDenom(dx, P)[1]

        for i in 2:(nCells-1)
            S = Cx * abs(pCentralGrad[i]) / pDenom[i]
            rho[i] = rho[i] -dxMomdx[i]*dt + S*rhoCG[i]
            xMom[i] = xMom[i] -drhoU2pdx[i]*dt + S*xMomCG[i]
            eV2[i] = eV2[i] -drhoUeV2PU[i]*dt + S*eV2CG[i]

            # Decode
            P[i], T[i], U[i] = decodePrimitives(rho[i], xMom[i], eV2[i])
            rhoU2p[i] = xMom[i]*U[i] + P[i]
            rhoUeV2PU[i] = U[i]*eV2[i] + P[i]*U[i]
        end

        ############### Boundaries ################
        # Waves never reach the boundaries, so boundary treatment doesn't need to be good
        allVars = [ rho, rhoU2p, rhoUeV2PU, xMom, eV2, U, P, T ]
        copyValues(3, 2, allVars)
        copyValues(2, 1, allVars)
        copyValues(nCells-2, nCells-1, allVars)
        copyValues(nCells, nCells, allVars)
        
        ############## CFL Calculation, timestep adjustment #############
        for i in 1:nCells
            CFL[i] = (abs(U[i]) + sqrt(gamma * R * T[i])) * dt / dx[i]
        end
        maxCFL = maximum(CFL)

        # Adjust time step to slowly approach target CFL
        dt *= ((targetCFL/maxCFL - 1)/5+1)

        currTime += dt
    end

    return P, U, T, rho
end

# U is a vector for FVM
function upwindFVM(mesh, P, T, U; initDt=0.001, endTime=0.14267, targetCFL=0.2, gamma=1.4, R=287.05, Cp=1005, Cx=0.3, debug=false)
    ######### MESH ############
    # Extract mesh into local variables for readability
    cells, faces, fAVecs, bdryFaces, cellVols = mesh
    nCells = size(cells, 1)
    nFaces = size(faces, 1)
    nBdryFaces = size(bdryFaces, 1)
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
        xMassFlux, xMomFlux, xeV2Flux = linInterp(mesh, xMom, rhoU2p, rhoUeV2PU)
        fluxVars = [ xMassFlux, xMomFlux, xeV2Flux ]
        
        # TODO: Update with finite volume version of artificial diffusivity
        dx = 1/nCells
        pCentralGrad, rhoCG, xMomCG, eV2CG = central2GradNum(dx, P, rho, xMom, eV2)
        pDenom = central2GradDenom(dx, P)[1]
        
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
                stateVars[v][ownerCell] -= fluxVars[v][i]*fA*dt/cellVols[ownerCell]
                stateVars[v][neighbourCell] += fluxVars[v][i]*fA*dt/cellVols[neighbourCell]
            end
        end

        for i in 2:nCells-1
            S = Cx * abs(pCentralGrad[i]) / pDenom[i]
            rho[i] += S*rhoCG[i]
            xMom[i] += S*xMomCG[i]
            eV2[i] += S*eV2CG[i]
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