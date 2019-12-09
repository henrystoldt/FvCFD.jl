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

function decodePrimitives3D!(primitives, cellState, R=287.05, Cp=1005)
    # Ux = xMom/rho
    primitives[3] = cellState[2]/cellState[1]
    # Uy = yMom/rho
    primitives[4] = cellState[3]/cellState[1]
    # Uz = zMom/rho
    primitives[5] = cellState[4]/cellState[1]
    # e = (eV2/rho) - (mag((Ux, Uy, Uz))^2)/2
    @views e = (cellState[5]/cellState[1]) - (mag(primitives[3:5])^2)/2
    #T
    primitives[2] = calPerfectT(e, Cp, R)
    #P = idealGasP(rho, P, R)
    primitives[1] = idealGasP(cellState[1], primitives[2], R)
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

#Fluxes, prim, state are each assumed to be in the same order as cellFluxes, cellPrimitives and cellState
#Now directly modified the flux vector instead of returning values
function calculateFluxes3D!(fluxes, prim, state)
    # Mass Fluxes
    fluxes[1] = state[2]
    fluxes[2] = state[3]
    fluxes[3] = state[4]

    # x-direction fluxes
    # xMomxFlux = xMom*Ux + P
    fluxes[4] = state[2]*prim[3] + prim[1]
    # yMomxFlux = xMom*Uy
    fluxes[7] = state[2]*prim[4]
    # zMomxFlux = xMom*Uz
    fluxes[10] = state[2]*prim[5]

    #y-direction fluxes
    # xMomyFlux = yMomxFlux
    fluxes[5] = fluxes[7]
    # yMomyFlux = yMom*Uy + P
    fluxes[8] = state[3]*prim[4] + prim[1]
    # zMomyFlux = yMom*Uz
    fluxes[11] = state[3]*prim[5]

    #z-direction fluxes
    # xMomzFlux = zMomxFlux
    fluxes[6] = fluxes[10]
    # yMomzFlux = zMomyFlux
    fluxes[9] = fluxes[11]
    # zMomzFlux = zMom*Uz + P
    fluxes[12] = state[4]*prim[5] + prim[1]

    # eV2xFlux = Ux*eV2 + P*Ux
    fluxes[13] = prim[3]*state[5] + prim[1]*prim[3]
    # eV2yFlux = Uy*eV2 + P*Uy
    fluxes[14] = prim[4]*state[5] + prim[1]*prim[4]
    # eV2zFlux = Uz*eV2 + P*Uz
    fluxes[15] = prim[5]*state[5] + prim[1]*prim[5]

    # return [ xMom, yMom, zMom, xMomxFlux, xMomyFlux, xMomzFlux, yMomxFlux, yMomyFlux, yMomzFlux, zMomxFlux, zMomyFlux, zMomzFlux, eV2xFlux, eV2yFlux, eV2zFlux ]
end
