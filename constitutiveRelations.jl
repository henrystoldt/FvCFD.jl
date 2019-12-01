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
    return [ P, T, U ]
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
    return [rho, xMom, eV2]
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

function calculateFluxes1D(P, Ux, xMom, eV2)
    massFlux = xMom
    xMomxFlux = xMom*Ux + P
    eV2xFlux = Ux*eV2 + P*Ux
    return [ massFlux, xMomxFlux, eV2xFlux ]
end
