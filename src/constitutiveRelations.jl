######################### Internal/Private functions #######################
function idealGasRho(T, P, R=287.05)
    @fastmath return P/(R*T)
end

function idealGasP(rho, T, R=287.05)
    # PV = mRT
    # P = rho R T
    @fastmath return rho*R*T
end

function calPerfectEnergy(T, Cp=1005, R=287.05)
    @fastmath return T*(Cp-R)
end

function calPerfectT(e, Cp=1005, R=287.05)
    @fastmath return e/(Cp-R)
end

######################### Public functions (called by finiteVolume.jl) #######################
#=
    Calculate primitives from cellState
=#
function decodePrimitives3D!(primitives, cellState, R=287.05, Cp=1005)
    ## Velocity ##
    # Ux = xMom/rho
    primitives[3] = cellState[2]/cellState[1]
    # Uy = yMom/rho
    primitives[4] = cellState[3]/cellState[1]
    # Uz = zMom/rho
    primitives[5] = cellState[4]/cellState[1]

    ## Energy ##
    # e = (eV2/rho) - (mag((Ux, Uy, Uz))^2)/2
    @views e = (cellState[5]/cellState[1]) - (mag(primitives[3:5])^2)/2

    ## Temperature ##
    primitives[2] = calPerfectT(e, Cp, R)

    ## Pressure ##
    # P = idealGasP(rho, P, R)
    primitives[1] = idealGasP(cellState[1], primitives[2], R)
end

#=
    Calculate cellState from cellPrimitives

    Arguments:
        cellPrimitives: matrix of primitives: see dataStructuresDefinitions.md
        R: R-value of the fluid (J/(kg*K))
        Cp: Constant pressure specific heat capacity (kJ/(kgK))

    Returns:
        cellState: 2D vector of conserved variable values: see dataStructuresDefinitions.md
=#
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

#=
    Calculates fluxes of transported variables at cell center, from
    Arguments:
        fluxes: (output) vector of cell center fluxes (to be calculated/populated)
        prim:   (input) vector of cell center primitives
        state:  (input) vector of cell center state

    Returns:
        None. Calculation results store in fluxes variable

    Notes:
        See dataStructuresDefinitions.md for definitions of state, primitives, etc...

=#
function calculateFluxes3D!(fluxes, prim, state)
    #### Mass Fluxes ####
    fluxes[1] = state[2]
    fluxes[2] = state[3]
    fluxes[3] = state[4]

    #### Momentum Fluxes ####
    ## x-direction momentum fluxes ##
    # xMomxFlux = xMom*Ux + P
    fluxes[4] = state[2]*prim[3] + prim[1]
    # yMomxFlux = xMom*Uy
    fluxes[7] = state[2]*prim[4]
    # zMomxFlux = xMom*Uz
    fluxes[10] = state[2]*prim[5]

    ## y-direction momentum fluxes ##
    # xMomyFlux = yMomxFlux
    fluxes[5] = fluxes[7]
    # yMomyFlux = yMom*Uy + P
    fluxes[8] = state[3]*prim[4] + prim[1]
    # zMomyFlux = yMom*Uz
    fluxes[11] = state[3]*prim[5]

    ## z-direction momentum fluxes ##
    # xMomzFlux = zMomxFlux
    fluxes[6] = fluxes[10]
    # yMomzFlux = zMomyFlux
    fluxes[9] = fluxes[11]
    # zMomzFlux = zMom*Uz + P
    fluxes[12] = state[4]*prim[5] + prim[1]

    #### Energy Fluxes ###
    # eV2xFlux = Ux*eV2 + P*Ux
    fluxes[13] = prim[3]*state[5] + prim[1]*prim[3]
    # eV2yFlux = Uy*eV2 + P*Uy
    fluxes[14] = prim[4]*state[5] + prim[1]*prim[4]
    # eV2zFlux = Uz*eV2 + P*Uz
    fluxes[15] = prim[5]*state[5] + prim[1]*prim[5]
end
