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
    T = calPerfectT(e, Cp, R)
    P = idealGasP(rho, T, R)
    return [ P, T, U ]
end

#Actually 1D
function decodePrimitives3D(rho::Float64, xMom::Float64, eV2::Float64, R=287.05, Cp=1005)
    U = [ xMom/rho, 0.0, 0.0 ]
    e = (eV2/rho) - (mag(U)^2)/2
    T = calPerfectT(e, Cp, R)
    P = idealGasP(rho, T, R)

    Ux = U[1]

    # x-Fluxes
    rhoU2p = xMom*Ux + P
    rhoUeV2PU = Ux*eV2 + P.*Ux

    return P, T, U, rhoU2p, rhoUeV2PU
end

function decodePrimitives3D(rho::Float64, xMom::Float64, yMom::Float64, zMom::Float64, eV2::Float64, R=287.05, Cp=1005)
    Ux = xMom/rho
    Uy = yMom/rho
    Uz = zMom/rho
    e = (eV2/rho) - (mag((Ux, Uy, Uz))^2)/2
    T = calPerfectT(e, Cp, R)
    P = idealGasP(rho, T, R)

    return [ P, T, Ux, Uy, Uz ]
end

# Returns rho, xMom, and eV2
function encodePrimitives(P, T, U, R=287.05, Cp=1005)
    rho = idealGasRho(T, P, R)
    xMom = U*rho
    e = calPerfectEnergy(T, Cp, R)
    eV2 = rho*(e + U*U/2)
    return [rho, xMom, eV2]
end

function encodePrimitives3D(P::Float64, T::Float64, U::Array{Float64, 1}, R=287.05, Cp=1005)
    # State Variables
    rho = idealGasRho(T, P, R)
    xMom = U[1]*rho
    e = calPerfectEnergy(T, Cp, R)
    eV2 = rho*(e + (mag(U)^2)/2)

    Ux = U[1]

    # x-Fluxes
    rhoU2p = xMom*Ux + P
    rhoUeV2PU = Ux*eV2 + P*Ux

    return rho, xMom, eV2, rhoU2p, rhoUeV2PU
end

function encodePrimitives3D(cellPrimitives::Array{Float64, 2}, R=287.05, Cp=1005)
    nCells = size(cellPrimitives, 1)

    cellState = zeros(nCells, 5)
    for c in 1:nCells
        cellState[c,1] = idealGasRho(cellPrimitives[c,2], cellPrimitives[c,1], R)
        cellState[c,2:4] .= cellPrimitives[c,3:5] .* cellState[c,1]
        e = calPerfectEnergy(cellPrimitives[c,2], Cp, R)
        cellState[c,5] = cellState[c,1]*(e + (mag(cellPrimitives[c,3:5])^2)/2 )
    end
    return cellState
end

function calculateFluxes1D(P, Ux, xMom, eV2)
    massFlux = xMom
    xMomxFlux = xMom*Ux + P
    eV2xFlux = Ux*eV2 + P*Ux
    return [ massFlux, xMomxFlux, eV2xFlux ]
end

function calculateFluxes3D(P, T, Ux, Uy, Uz, rho, xMom, yMom, zMom, eV2)
    # x-direction fluxes
    xMomxFlux = xMom*Ux + P
    yMomxFlux = xMom*Uy
    zMomxFlux = xMom*Uz

    #y-direction fluxes
    xMomyFlux = yMomxFlux
    yMomyFlux = yMom*Uy + P
    zMomyFlux = yMom*Uz

    #z-direction fluxes
    xMomzFlux = zMomxFlux
    yMomzFlux = zMomyFlux
    zMomzFlux = zMom*Uz + P

    eV2xFlux = Ux*eV2 + P*Ux
    eV2yFlux = Uy*eV2 + P*Uy
    eV2zFlux = Uz*eV2 + P*Uz

    return [ xMom, yMom, zMom, xMomxFlux, xMomyFlux, xMomzFlux, yMomxFlux, yMomyFlux, yMomzFlux, zMomxFlux, zMomyFlux, zMomzFlux, eV2xFlux, eV2yFlux, eV2zFlux ]
end
