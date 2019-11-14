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

function calPerfectEnergy(T, Cp=1005)
    return T*Cp
end

function calPerfectT(e, Cp=1005)
    return e/Cp
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
    rho = Array{Float64, 1}(undef, nCells)
    U = Array{Float64, 1}(undef, nCells)
    e = Array{Float64, 1}(undef, nCells)
    P = Array{Float64, 1}(undef, nCells)
    T = Array{Float64, 1}(undef, nCells)


    # Apply initial conditions (Fig. 1 in Henry's paper)
    for i in 1:nCells
        LHSRho = idealGasRho(0.00348432, 1)
        LHSe = calPerfectEnergy(0.00348432)
        RHSRho = idealGasRho(0.00278746, 0.1)
        RHSe = calPerfectEnergy(0.00278746)

        if i <= (nCells/2)
            rho[i] = LHSRho
            U[i] = 0
            e[i] = LHSe
        else
            rho[i] = RHSRho
            U[i] = 0
            e[i] = RHSe
        end
        T[i] = calPerfectT(e[i])
        P[i] = idealGasP(rho[i], T[i])

        # dx = Distance between cell centers i and i+1
        dx[i] = domainLength / nCells
    end

    return dx, rho, U, e, P, T
end

######################### Solution #######################
# Pass in initial values for each variable
# Shock Tube (undisturbed zero gradient) boundary conditions assumed
function macCormack1D(dx, rho, U, e, P, T; dt=0.001, endTime=0.14267)
    nCells = size(rho, 1)
    
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

#TODO: Min initial conditions
function macCormack1DConservative(dx, rho, U, e, P, T; dt=0.001, endTime=0.14267)
    nCells = size(rho, 1)

    xMom = Array{Float64, 1}(undef, nCells)
    eV2 = Array{Float64, 1}(undef, nCells)

    #TODO: calc xMom and eV2 from initCond
    
    drhoPred = Array{Float64, 1}(undef, nCells)
    dxMPred = Array{Float64, 1}(undef, nCells)
    deVPred = Array{Float64, 1}(undef, nCells)

    rhoPred = Array{Float64, 1}(undef, nCells)
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
        drhodx, dudx, dpdx, dedx, deVdx = backwardGradient(dx, rho, U, P, e, eV2)
        # drhodx, dudx, dpdx, dedx = upwindGradient(dx, U, rho, U, P, e)
        # drhodx, dudx, dpdx, dedx = centralGradient(dx, rho, U, P, e)

        for i in 2:(nCells-1)

            # Eq. 2.82b
            drhoPred[i] = -(rho[i]*dudx[i] + U[i]*drhodx[i])
            # Eq. 2.84a
            dxMPred[i] = -(dpdx[i] +  2*rho[i]*U[i]*dudx[i] + U[i]*U[i]*drhodx[i])
            # Eq. 2.86
            deVPred[i] = -(U[i]*dpdx[i] + P[i]*dudx[i] + eV2[i]*drhodx[i]*U[i] + rho[i]*deVdx[i]*U[i] + rho[i]*eV2[i]*dudx[i])

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

        #TODO: Update block
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
        #TODO: Update
        allVars = [ rho, U, e, T, P ]
        copyValues(3, 2, allVars)
        copyValues(2, 1, allVars)
        copyValues(nCells-2, nCells-1, allVars)
        copyValues(nCells, nCells, allVars)
        
        currTime += dt
    end

    return P, U, T, rho
end

# Pass in initial values for each variable
# Shock Tube (undisturbed zero gradient) boundary conditions assumed
function upwind1D(dx, rho, U, e, P, T; dt=0.001, endTime=0.14267)
    nCells = size(rho, 1)
    
    drhoPred = Array{Float64, 1}(undef, nCells)
    duPred = Array{Float64, 1}(undef, nCells)
    dePred = Array{Float64, 1}(undef, nCells)

    currTime = 0
    while currTime < endTime
        if (endTime - currTime) < dt
            dt = endTime - currTime
        end

        ############## Predictor #############
        # Forward Differences to compute gradients
        # drhodx, dudx, dpdx, dedx = backwardGradient(dx, rho, U, P, e)
        drhodx, dudx, dpdx, dedx = upwindGradient(dx, U, rho, U, P, e)
        # drhodx, dudx, dpdx, dedx = centralGradient(dx, rho, U, P, e)
        for i in 2:(nCells-1)
            drhoPred[i] = -(rho[i]*dudx[i] + U[i]*drhodx[i])
            duPred[i] = -(U[i]*dudx[i] + dpdx[i]/rho[i])
            dePred[i] = -(U[i]*dedx[i] + P[i]*dudx[i]/rho[i])

            rho[i] = rho[i] + drhoPred[i]*dt
            U[i] = U[i] + duPred[i]*dt
            e[i] = e[i] + dePred[i]*dt
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

P, U, T, rho = upwind1D(initializeShockTube()..., endTime=0.01)

print(P)

################## Output ##################
plots = []
pPlot = plot(P, label="P (Pa)", title="Pressure", xlabel="x (m)", bottom_margin=15mm, left_margin=10mm)
rhoPlot = plot(rho, label="rho (kg/m3)", title="Density", xlabel="x (m)", bottom_margin=15mm, left_margin=10mm)
uPlot = plot(U, label="Velocity (m/s)", title="Velocity", xlabel="x (m)", bottom_margin=15mm, left_margin=10mm)
TPlot = plot(T, label="T (K)", title="Temperature", xlabel="x (m)", bottom_margin=15mm, left_margin=10mm)
plots = [pPlot, rhoPlot, uPlot, TPlot]
plot(plots..., layout=(2, 2), size=(1720, 880), window_title="Shit Euler = Sheuler", legend=false)