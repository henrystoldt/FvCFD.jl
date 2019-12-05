function checkInputs_structured1DInterp(dx, faceValues::Array{Float64,2}, cellValues::Array{Float64, 2})
    nCells = size(cellValues, 1)
    nFaces = size(faceValues, 1)
    nVars = size(faceValues, 2)
    nVars2 = size(cellValues, 2)
    if nCells != nFaces-1 || nVars2 != nVars
        throw(ArgumentError("ERROR: nCells ($nCells) != nFaces ($nFaces) or Number of variables stored in cells($nVars2) and faces($nVars) are not equal!"))
    end
    nDxs = size(dx, 1)
    if nDxs != nCells
        throw(ArgumentError("ERROR: nCells ($nCells) != nDxs ($nDxs)"))
    end
end

### Unused ###
function structured_1DlinInterp(dx, faceValues::Array{Float64,2}, cellValues::Array{Float64, 2})
    nFaces = size(faceValues, 1)
    nVars = size(faceValues, 2)
    checkInputs_structured1DInterp(dx, faceValues, cellValues)

    # Do interpolation
    for f in 2:nFaces-1
        for v in 1:nVars
            c1Dist = dx[f-1]/2
            c2Dist = dx[f]/2
            totalDist = c1Dist + c2Dist
            faceValues[f, v] = cellValues[f-1, v].*(c2Dist/totalDist) .+ cellValues[f, v].*(c1Dist/totalDist)
        end
    end
end

function structured_1DMaxInterp(dx, faceValues::Array{Float64,2}, cellValues::Array{Float64, 2})
    nFaces = size(faceValues, 1)
    nVars = size(faceValues, 2)
    checkInputs_structured1DInterp(dx, faceValues, cellValues)

    # Do interpolation
    for f in 2:nFaces-1
        for v in 1:nVars
            faceValues[f, v] = max(cellValues[f-1,v],  cellValues[f,v])
        end
    end
end

function structured_1DFaceDelta(dx, faceValues::Array{Float64,2}, cellValues::Array{Float64, 2})
    nFaces = size(faceValues, 1)
    nVars = size(faceValues, 2)
    checkInputs_structured1DInterp(dx, faceValues, cellValues)

    # Do interpolation
    for f in 2:nFaces-1
        for v in 1:nVars
            faceValues[f, v] = cellValues[f, v] - cellValues[f-1, v]
        end
    end
end

### Used ###
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

# Pass in pressure or density for the pressure-based or density-based version respectively
# Density-based version detects contact discontinuities
function structured_1D_JST_sj2(dvar)
    nCells = size(dvar, 1)
    sj = Array{Float64, 1}(undef, nCells)

    # Boundary values unset
    sj[1] = 0.0
    sj[nCells] = 0.0

    for c in 2:nCells-1
        sj[c] = (abs( dvar[c] - dvar[c-1] )/ max( abs(dvar[c]) + abs(dvar[c-1]), 0.0000000001))^2
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

function structured_1D_JST_Eps(dx, k2, k4, c4, cellPrimitives::Array{Float64,2}, gamma=1.4, R=287.05)
    nFaces = size(cellPrimitives, 1) + 1

    # Pressure sensor seems to work best for the shock tube
    dP = structured_1DFaceDelta(dx, cellPrimitives[:,1])[1]
    sj = structured_1D_JST_sj2(dP)

    rj = structured_1D_JST_rj(cellPrimitives[:,2], cellPrimitives[:,3], gamma, R)
    sjF, rjF = structured_1DMaxInterp(dx, sj, rj)

    eps2 = Array{Float64, 1}(undef, nFaces)
    eps4 = Array{Float64, 1}(undef, nFaces)
    for f in 2:nFaces-1
        eps2[f] = k2 * sjF[f] * rjF[f]
        eps4[f] = max(0, k4*rjF[f] - c4*eps2[f])
    end

    return eps2, eps4
end

# Requires correct cellState and cellPrimitives as input
function structured_JSTFlux1D(dx, solutionState, boundaryConditions, gamma, R)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState
    nCells = size(cellState, 1)
    nFaces = nCells + 1
    faceDeltas = Array{Float64, 2}(undef, nFaces, 3)

    # Centrally differenced fluxes
    structured_1DlinInterp(dx, faceFluxes, cellFluxes)

    #### Add JST artificial Diffusion ####
    structured_1DFaceDelta(dx, faceDeltas, cellState)
    eps2, eps4 = structured_1D_JST_Eps(dx, 0.5, (1.2/32), 0.9, cellPrimitives)
    # nCells = nFaces - 1
    for f in 2:(nCells)
        for v in 1:3
            diffusionFlux = eps2[f]*faceDeltas[f, v] - eps4[f]*(faceDeltas[f+1, v] - 2*faceDeltas[f, v] + faceDeltas[f-1, v])
            faceFluxes[f,v] -= diffusionFlux
        end
    end

    return integrateFluxes_structured1D(dx, solutionState)
end
