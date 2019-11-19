using Plots
using Plots.PlotMeasures
using LaTeXStrings
plotly()

######################### Constitutive Relations #######################
function idealGasRho(T, P, R=287.05)
    # PV = nRT
    # PV = mRT
    # m/V = P/RT = rho
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

function calPerfectT(e, Cp=1005)
    return e/Cp
end

function decodePrimitives(rho, xMom, eV2, R=287.05, Cp=1005)
    U = xMom / rho
    e = (eV2/rho) - U*U/2
    T = calPerfectT(e, Cp)
    P = idealGasP(rho, T, R)
    return P, T, U
end

function encodePrimitives(P, T, U, R=287.05, Cp=1005)
    rho = idealGasRho(T, P)
    xMom = U * rho
    e = calPerfectEnergy(T)
    eV2 = rho * (e + U*U/2)
    return rho, xMom, eV2
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

#TODO: TVD gradients

######################### Boundary Conditions #######################
function copyValues(fromIndex, toIndex, varArrays)
    for varArray in varArrays
        varArray[toIndex] = varArray[fromIndex]
    end
end

######################### Initialization #######################
function initializeShockTube(nCells=100, domainLength=1)
    # Create arrays to store data (cell # = position in array)
    dx = Array{Float64, 1}(undef, nCells)
    U = Array{Float64, 1}(undef, nCells)
    P = Array{Float64, 1}(undef, nCells)
    T = Array{Float64, 1}(undef, nCells)

    # Apply initial conditions (Fig. 1 in Henry's paper)
    LHSRho = idealGasRho(0.00348432, 1)
    LHSe = calPerfectEnergy(0.00348432)
    RHSRho = idealGasRho(0.00278746, 0.95)
    RHSe = calPerfectEnergy(0.00278746)
    for i in 1:nCells

        if i <= (nCells/2)
            rho = LHSRho
            e = LHSe
        else
            rho = RHSRho
            e = RHSe
        end
        U[i] = 0
        T[i] = calPerfectT(e)
        P[i] = idealGasP(rho, T[i])

        # dx = Distance between cell centers i and i+1
        dx[i] = domainLength / nCells
    end

    return dx, P, T, U
end

######################### Solution #######################
# Pass in initial values for each variable
# Shock Tube (undisturbed zero gradient) boundary conditions assumed
function macCormack1D(dx, P, T, U; dt=0.001, endTime=0.14267)
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

    currTime = 0
    while currTime < endTime
        if (endTime - currTime) < dt
            dt = endTime - currTime
        end

        ############## Predictor #############
        # Forward Differences to compute gradients
        drhodx, dudx, dpdx, dedx = backwardGradient(dx, rho, U, P, e)
        # drhodx, dudx, dpdx, dedx = upwindGradient(dx, U, rho, U, P, e)
        # drhodx, dudx, dpdx, dedx = centralGradient(dx, rho, U, P, e)
        for i in 2:(nCells-1)

            # Eq. 6.1
            drhoPred[i] = -(rho[i]*dudx[i] + U[i]*drhodx[i])
            # Eq. 6.2
            duPred[i] = -(U[i]*dudx[i] + dpdx[i]/rho[i])
            # Eq. 6.4
            dePred[i] = -(U[i]*dedx[i] + P[i]*dudx[i]/rho[i])

            rhoPred[i] = rho[i] + drhoPred[i]*dt
            UPred[i] = U[i] + duPred[i]*dt
            ePred[i] = e[i] + dePred[i]*dt
            TPred[i] = calPerfectT(ePred[i])
            PPred[i] = idealGasP(rhoPred[i], TPred[i])
        end

        ############### Corrector ################
        # Rearward differences to compute gradients
        drhodx, dudx, dpdx, dedx = forwardGradient(dx, rhoPred, UPred, PPred, ePred)
        # drhodx, dudx, dpdx, dedx = upwindGradient(dx, UPred, rhoPred, UPred, PPred, ePred)
        # drhodx, dudx, dpdx, dedx = centralGradient(dx, rhoPred, UPred, PPred, ePred)
        for i in 2:(nCells-1)

            # Eq. 6.1
            drhoPred2 = -(rhoPred[i]*dudx[i] + UPred[i]*drhodx[i])
            # Eq. 6.2
            duPred2 = -(UPred[i]*dudx[i] + dpdx[i]/rhoPred[i])
            # Eq. 6.4
            dePred2 = -(UPred[i]*dedx[i] + PPred[i]*dudx[i]/rhoPred[i])

            # Perform timestep using average gradients
            rho[i] += (drhoPred2 + drhoPred[i])*dt/2
            U[i] += (duPred2 + duPred[i])*dt/2
            e[i] += (dePred2 + dePred[i])*dt/2
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
        
        currTime += dt
    end

    return P, U, T, rho
end

#TODO: CFL calculation
function macCormack1DConservative(dx, P, T, U; initDt=0.001, endTime=0.14267, targetCFL=0.5, gamma=1.4, R=287.05, Cp=1005)
    nCells = size(dx, 1)

    rho = Array{Float64, 1}(undef, nCells)
    xMom = Array{Float64, 1}(undef, nCells)
    eV2 = Array{Float64, 1}(undef, nCells)

    for i in 1:nCells
        rho[i], xMom[i], eV2[i] = encodePrimitives(P[i], T[i], U[i])
    end
    
    drhoP = Array{Float64, 1}(undef, nCells)
    dxMP = Array{Float64, 1}(undef, nCells)
    deV2P = Array{Float64, 1}(undef, nCells)

    xMomP = Array{Float64, 1}(undef, nCells)
    eV2P = Array{Float64, 1}(undef, nCells)
    rhoP = Array{Float64, 1}(undef, nCells)
    UP = Array{Float64, 1}(undef, nCells)
    PP = Array{Float64, 1}(undef, nCells)
    TP = Array{Float64, 1}(undef, nCells)

    CFL = Array{Float64, 1}(undef, nCells)

    dt = initDt
    currTime = 0
    while currTime < endTime
        if (endTime - currTime) < dt
            dt = endTime - currTime
        end
        
        ############## Predictor #############
        # Forward Differences to compute gradients
        drhodx, dudx, dpdx, deVdx = forwardGradient(dx, rho, U, P, eV2)
        # drhodx, dudx, dpdx, deVdx = upwindGradient(dx, U, rho, U, P, eV2)
        # drhodx, dudx, dpdx, deVdx = centralGradient(dx, rho, U, P, eV2)

        for i in 2:(nCells-1)
            # Eq. 2.82b
            drhoP[i] = -(rho[i]*dudx[i] + U[i]*drhodx[i])
            # Eq. 2.84a
            dxMP[i] = -(dpdx[i] +  2*rho[i]*U[i]*dudx[i] + U[i]*U[i]*drhodx[i])
            # Eq. 2.86
            deV2P[i] = -(U[i]*dpdx[i] + P[i]*dudx[i] + eV2[i]*drhodx[i]*U[i] + rho[i]*deVdx[i]*U[i] + rho[i]*eV2[i]*dudx[i])

            # Predict
            rhoP[i] = rho[i] + drhoP[i]*dt
            xMomP[i] = xMom[i] + dxMP[i]*dt
            eV2P[i] = eV2[i] + deV2P[i]*dt

            # Decode
            PP[i], TP[i], UP[i] = decodePrimitives(rhoP[i], xMomP[i], eV2P[i])
        end

        ############### Corrector ################
        # Rearward differences to compute gradients
        drhodxP, dudxP, dpdxP, deVdxP = backwardGradient(dx, rhoP, UP, PP, eV2P)
        # drhodxP, dudxP, dpdxP, deVdxP = upwindGradient(dx, UP, rhoP, UP, PP, eV2P)
        # drhodxP, dudxP, dpdxP, deVdxP = centralGradient(dx, rhoP, UP, PP, eV2P)

        for i in 2:(nCells-1)
            drhoP2 = -(rhoP[i]*dudxP[i] + UP[i]*drhodxP[i])
            dxMP2 = -(dpdxP[i] +  2*rhoP[i]*UP[i]*dudxP[i] + UP[i]*UP[i]*drhodxP[i])
            deV2P2 = -(UP[i]*dpdxP[i] + PP[i]*dudxP[i] + eV2P[i]*drhodxP[i]*UP[i] + rhoP[i]*deVdxP[i]*UP[i] + rhoP[i]*eV2P[i]*dudxP[i])

            # Perform timestep using average gradients
            rho[i] += (drhoP2 + drhoP[i])*dt/2
            xMom[i] += (dxMP2 + dxMP[i])*dt/2
            eV2[i] += (deV2P2 + deV2P[i])*dt/2

            # Decode
            P[i], T[i], U[i] = decodePrimitives(rho[i], xMom[i], eV2[i])
        end

        ############### Boundaries ################
        # Waves never reach the boundaries, so boundary treatment doesn't need to be good
        allVars = [ rho, xMom, eV2, U, T, P ]
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

################## Output ##################
nCells = 500
P, U, T, rho = macCormack1DConservative(initializeShockTube(nCells)..., initDt=0.0001, endTime=0.2)

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