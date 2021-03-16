######################### Boundary Conditions #########################
# The job of the BC functions is to calculate/enforce face fluxes at boundary faces.
# The JST method only calculates fluxes at internal faces, then these conditions are applied to calculate them at the boundaries

# InletConditions: [ Static Pressure, Static Temperture, Ux, Uy, Uz, Cp ]
function supersonicInletBoundary(mesh::Mesh, sln::SolutionState, boundaryNumber, inletConditions)
    P, T, Ux, Uy, Uz, Cp = inletConditions
    rho = idealGasRho(T, P)
    xMom, yMom, zMom = [Ux, Uy, Uz] .* rho
    e = calPerfectEnergy(T, Cp)
    eV2 = rho*(e + (mag([Ux, Uy, Uz])^2)/2)

    boundaryFluxes = Vector{Float64}(undef, 15)
    calculateFluxes3D!(boundaryFluxes, [P, T, Ux, Uy, Uz],  [rho, xMom, yMom, zMom, eV2])
    currentBoundary = mesh.boundaryFaces[boundaryNumber]
    @inbounds for face in currentBoundary
        sln.faceFluxes[face,:] = boundaryFluxes
    end
end

# InletConditions: [ totalPressure, totalTemp, nx, ny, nz, gamma, R, Cp ]
# Where n is the unit vector representing the direction of inlet velocity
# Using method from FUN3D solver
function subsonicInletBoundary(mesh, sln::SolutionState, boundaryNumber, inletConditions)
    Pt, Tt, nx, ny, nz, gamma, R, Cp = inletConditions

    boundaryFluxes = Vector{Float64}(undef, 15)
    primitives = Vector{Float64}(undef, 5)
    state = Vector{Float64}(undef, 5)
    velUnitVector = [ nx, ny, nz ]
    currentBoundary = mesh.boundaryFaces[boundaryNumber]

    @inbounds @fastmath for face in currentBoundary
        ownerCell = mesh.faces[face][1]
        @views adjustedVelocity = dot(velUnitVector, sln.cellPrimitives[ownerCell, 3:5])
        primitives[2] = Tt - (gamma-1)/2 * (adjustedVelocity^2)/(gamma*R)
        machNum = abs(adjustedVelocity) / sqrt(gamma * R * primitives[2])
        primitives[1] = Pt*(1 + (gamma-1)/2 * machNum^2)^(-gamma/(gamma-1))
        primitives[3:5] = adjustedVelocity .* velUnitVector

        # Calculate state variables
        state[1] = idealGasRho(primitives[2], primitives[1], R)
        state[2:4] .= primitives[3:5] .* state[1]
        e = calPerfectEnergy(primitives[2], Cp, R)
        state[5] = state[1]*(e + (mag(primitives[3:5])^2)/2 )

        calculateFluxes3D!(boundaryFluxes, primitives, state)
        sln.faceFluxes[face, :] = boundaryFluxes
    end
end

function pressureOutletBoundary(mesh, sln::SolutionState, boundaryNumber, outletPressure)
    nFluxes = size(sln.cellFluxes, 2)

    # Directly extrapolate cell center flux to boundary (zero gradient between the cell center and the boundary)
    currentBoundary = mesh.boundaryFaces[boundaryNumber]
    @inbounds @fastmath for face in currentBoundary
        ownerCell = mesh.faces[face][1]
        # Mass fluxes are unaffected by pressure boundary
        for flux in 1:nFluxes
            sln.faceFluxes[face, flux] = sln.cellFluxes[ownerCell, flux]
        end

        origP = sln.cellPrimitives[ownerCell, 1]
        # Adjust the momentum fluxes containing pressure
        for flux in [4, 8, 12]
            sln.faceFluxes[face, flux] += outletPressure - origP
        end
        #Adjust the energy fluxes
        for d in 1:3
            sln.faceFluxes[face, 12+d] += sln.cellPrimitives[ownerCell, 2+d]*(outletPressure - origP)
        end
    end
end

function zeroGradientBoundary(mesh::Mesh, sln::SolutionState, boundaryNumber, _)
    nFluxes = size(sln.cellFluxes, 2)

    # Directly extrapolate cell center flux to boundary (zero gradient between the cell center and the boundary)
    currentBoundary = mesh.boundaryFaces[boundaryNumber]
    @inbounds for face in currentBoundary
        ownerCell = mesh.faces[face][1] #One of these will be -1 (no cell), the other is the boundary cell we want
        for flux in 1:nFluxes
            sln.faceFluxes[face, flux] = sln.cellFluxes[ownerCell, flux]
        end
    end
end

function wallBoundary(mesh::Mesh, sln::SolutionState, boundaryNumber, _)
    currentBoundary = mesh.boundaryFaces[boundaryNumber]
    @inbounds for f in currentBoundary
        ownerCell = max(mesh.faces[f][1], mesh.faces[f][2]) #One of these will be -1 (no cell), the other is the boundary cell we want

        faceP = sln.cellPrimitives[ownerCell, 1]
        # Momentum flux is Pressure in each of the normal directions (dot product)
        sln.faceFluxes[f, 4] = faceP
        sln.faceFluxes[f, 8] = faceP
        sln.faceFluxes[f, 12] = faceP

        # Mass Flux is zero
        sln.faceFluxes[f, 1:3] .= 0.0
        # Energy Flux is zero
        sln.faceFluxes[f, 13:15] .= 0.0
    end
end

# Symmetry is identical to a wall in Euler simulations
symmetryBoundary = wallBoundary

function emptyBoundary(mesh::Mesh, sln::SolutionState, boundaryNumber, _)
    return
end

function applyBoundaryConditions(mesh::Mesh, sln::SolutionState, boundaryConditions, nBoundaries)
    for boundaryNumber in 1:nBoundaries
        bFunctionIndex = 2*boundaryNumber-1
        boundaryConditionFunction = boundaryConditions[bFunctionIndex]
        boundaryParameters = boundaryConditions[bFunctionIndex+1]

        boundaryConditionFunction(mesh, sln, boundaryNumber, boundaryParameters)
    end
end