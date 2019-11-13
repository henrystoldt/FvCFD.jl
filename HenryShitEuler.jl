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
#TODO

######################### Boundary Conditions #######################
function copyValues(fromIndex, toIndex, varArrays)
    for varArray in varArrays
        varArray[toIndex] = varArray[fromIndex]
    end
end

######################### Initialization #######################
function doCFD()
    # Create arrays to store data (cell # = position in array)
    nCells = 100
    domainLength = 1 #m
    dt = 0.0001
    # endTime = 0.1
    endTime = 0.14267 #sec
    cellSize = Array{Float64, 1}(undef, nCells)
    rho = Array{Float64, 1}(undef, nCells)
    e = Array{Float64, 1}(undef, nCells)
    P = Array{Float64, 1}(undef, nCells)
    U = Array{Float64, 1}(undef, nCells)
    T = Array{Float64, 1}(undef, nCells)


    # Apply initial conditions (Fig. 1 in Henry's paper)
    for i in 1:nCells
        LHSRho = idealGasRho(0.00348432, 1)
        LHSe = calPerfectEnergy(0.00348432)
        RHSRho = idealGasRho(0.00278746, 0.5)
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

        cellSize[i] = domainLength / nCells
    end


    ######################### Solution #######################
    drhoPred = Array{Float64, 1}(undef, nCells)
    duPred = Array{Float64, 1}(undef, nCells)
    dePred = Array{Float64, 1}(undef, nCells)

    rhoPred = Array{Float64, 1}(undef, nCells)
    TPred = Array{Float64, 1}(undef, nCells)
    uPred = Array{Float64, 1}(undef, nCells)
    PPred = Array{Float64, 1}(undef, nCells)
    ePred = Array{Float64, 1}(undef, nCells)

    currTime = 0
    while currTime < endTime
        if (endTime - currTime) < dt
            dt = endTime - currTime
        end

        for i in 2:(nCells-1)
            # Forward Differences to compute gradients
            dx = (cellSize[i] + cellSize[i+1]) / 2
            drhodx = (rho[i+1] - rho[i]) / dx
            dudx = (U[i+1] - U[i]) / dx
            dpdx = (P[i+1] - P[i]) / dx
            dedx = (e[i+1] - e[i]) / dx

            # Eq. 6.1
            drhoPred[i] = -(rho[i]*dudx + U[i]*drhodx)
            # Eq. 6.2
            duPred[i] = -(U[i]*dudx + dpdx/rho[i])
            # Eq. 6.4
            dePred[i] = -(U[i]*dedx + P[i]*dudx/rho[i])

            rhoPred[i] = rho[i] + drhoPred[i]*dt
            uPred[i] = U[i] + duPred[i]*dt
            ePred[i] = e[i] + dePred[i]*dt
            TPred[i] = calPerfectT(ePred[i])
            PPred[i] = idealGasP(rhoPred[i], TPred[i])
        end

        for i in 2:(nCells-1)
            # Rearward differences to compute gradients
            dx = (cellSize[i] + cellSize[i-1]) / 2
            dudx = (uPred[i] - uPred[i-1]) / dx
            dpdx = (PPred[i] - PPred[i-1]) / dx
            drhodx = (rhoPred[i] - rhoPred[i-1]) / dx
            dedx = (ePred[i] - ePred[i-1]) / dx

            # Eq. 6.1
            drhoPred2 = -(rhoPred[i]*dudx + uPred[i]*drhodx)
            # Eq. 6.2
            duPred2 = -(uPred[i]*dudx + dpdx/rhoPred[i])
            # Eq. 6.4
            dePred2 = -(uPred[i]*dedx + PPred[i]*dudx/rhoPred[i])

            rho[i] += (drhoPred2 + drhoPred[i])*dt/2
            U[i] += (duPred2 + duPred[i])*dt/2
            e[i] += (dePred2 + dePred[i])*dt/2
            T[i] = calPerfectT(e[i])
            P[i] = idealGasP(rho[i], T[i])
        end

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

P, U, T, rho= doCFD()

plots = []
pPlot = plot(P, label="P (Pa)", title="Pressure", xlabel="x (m)", bottom_margin=15mm, left_margin=10mm)
rhoPlot = plot(rho, label="rho (kg/m3)", title="Density", xlabel="x (m)", bottom_margin=15mm, left_margin=10mm)
uPlot = plot(U, label="Velocity (m/s)", title="Velocity", xlabel="x (m)", bottom_margin=15mm, left_margin=10mm)
TPlot = plot(T, label="T (K)", title="Temperature", xlabel="x (m)", bottom_margin=15mm, left_margin=10mm)
plots = [pPlot, rhoPlot, uPlot, TPlot]
plot(plots..., layout=(2, 2), size=(1720, 880), window_title="Shit Euler = Sheuler", legend=false)