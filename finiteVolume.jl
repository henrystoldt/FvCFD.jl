using Printf
include("constitutiveRelations.jl")
include("vectorFunctions.jl")
include("timeDiscretizations.jl")
include("mesh.jl")
include("output.jl")
include("dataStructures.jl")

__precompile__()

######################### Initialization ###########################
# Returns cellPrimitives matrix for uniform solution
function initializeUniformSolution3D(mesh, P, T, Ux, Uy, Uz)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

    cellPrimitives = zeros(nCells, 5)
    for c in 1:nCells
        cellPrimitives[c, :] = [ P, T, Ux, Uy, Uz ]
    end

    return cellPrimitives
end

######################### CFL #######################################
# function maxCFL3D(mesh::Mesh, sln::SolutionState, dt, gamma=1.4, R=287.05)
#     nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
#
#     maxCFL = 0.0
#     for c in 1:nCells
#         a = sqrt(gamma * R * sln.cellPrimitives[c,2])
#
#         CFL = 0.0
#         for d in 1:3
#             cDeltaT = (abs(sln.cellPrimitives[c,d+2]) + a)*dt
#             CFL += cDeltaT / mesh.cellSizes[c, d]
#         end
#         maxCFL = max(maxCFL, CFL)
#     end
#
#     return maxCFL
# end

# Calculates CFL at each cell. Expects sln.cellState, sln.cellPrimitives and sln.faceFluxes to be up to date
function CFL!(CFL, mesh::Mesh, sln::SolutionState, dt=1, gamma=1.4, R=287.05)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

    fill!(CFL, 0.0)
    faceRhoT = linInterp_3D(mesh, hcat(sln.cellState[:,1], sln.cellPrimitives[:,2]))
    for f in 1:nFaces
        ownerCell = mesh.faces[f][1]
        neighbourCell = mesh.faces[f][2]

        faceRho = faceRhoT[f, 1]
        faceT = faceRhoT[f, 2]

        # Use cell center values on the boundary
        if neighbourCell == -1
            faceRho = sln.cellState[ownerCell, 1]
            faceT = sln.cellPrimitives[ownerCell, 2]
        end

        @views faceVel = sln.faceFluxes[f, 1:3] ./ faceRho
        @views flux = abs(dot(faceVel, mesh.fAVecs[f]))*dt

        if faceT <= 0.0
            println("Found it!")
        end

        a = sqrt(gamma * R * faceT)
        flux += mag(mesh.fAVecs[f])*a*dt

        CFL[ownerCell] += flux
        if neighbourCell > -1
            CFL[neighbourCell] += flux
        end
    end

    CFL ./= (2 .* mesh.cVols)
end

######################### Gradient Computation #######################
#TODO: make all these function return a single array if you pass in a single value
#TODO: LSQ Gradient
# Non-functional
function leastSqGrad(mesh::Mesh, matrix::AbstractArray{Float64, 2}, stencil=zeros(2,2))
    # Stencil should be a list of lists, with each sublist containing the cells contained in the stencil of the main cell


    #cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    #bdryFaceIndices = Array(nFaces-nBdryFaces:nFaces)

    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    nVars = size(matrix, 2) # Returns the length of the second dimension of "matrix"

    grad = zeros(nCells, nVars, 3)

    L11 = zeros(nCells,nVars)
    L12 = zeros(nCells,nVars)
    L13 = zeros(nCells,nVars)
    L22 = zeros(nCells,nVars)
    L23 = zeros(nCells,nVars)
    L33 = zeros(nCells,nVars)
    L1f = zeros(nCells,nVars)
    L2f = zeros(nCells,nVars)
    L3f = zeros(nCells,nVars)

    @fastmath for f in 1:nFaces-nBdryFaces

        c1 = mesh.faces[f][1]
        c2 = mesh.faces[f][2]

        dx = mesh.cCenter[c2][1] - mesh.cCenter[c1][1]
        dy = mesh.cCenter[c2][2] - mesh.cCenter[c1][2]
        dz = mesh.cCenter[c2][3] - mesh.cCenter[c1][3]

        weight = 1 / sqrt(dx*dx + dy*dy + dz*dz)

        wdx = weight * dx
        wdy = weight * dy
        wdz = weight * dz

        for v in 1:nVars
            dv = matrix[c2,v] - matrix[c1,v]
            wdv = weight * dv

            L11[c1,v] += wdx^2
            L12[c1,v] += wdx * wdy
            L13[c1,v] += wdx * wdz
            L22[c1,v] += wdy^2
            L23[c1,v] += wdy * wdz
            L33[c1,v] += wdz^2

            L1f[c1,v] += wdx * wdv
            L2f[c1,v] += wdy * wdv
            L3f[c1,v] += wdz * wdv

            L11[c2,v] += wdx^2
            L12[c2,v] += wdx * wdy
            L13[c2,v] += wdx * wdz
            L22[c2,v] += wdy^2
            L23[c2,v] += wdy * wdz
            L33[c2,v] += wdz^2

            L1f[c2,v] += wdx * wdv
            L2f[c2,v] += wdy * wdv
            L3f[c2,v] += wdz * wdv

        end

    end

    # Deal with boundary faces, and add them to the matrix vectors



    return grad
end

#=
    Takes the gradient of (scalar) data provided in matrix form (passed into arugment 'matrix'):
    Cell      x1      x2      x3
    Cell 1    x1_1    x2_1    x3_1
    Cell 2    x1_2    x2_2    x3_2
    ...
    (where x1, x2, and x3 are arbitrary scalars)
    and output a THREE-DIMENSIONAL gradient arrays of the following form
    Cell      x1          x2          x3
    Cell 1    grad(x1)_1  grad(x2)_1  grad(x3)_1
    Cell 2    grad(x1)_2  grad(x2)_2  grad(x3)_2
    ...
    Where each grad(xa)_b is made up of THREE elements for the (x,y,z) directions

    Ex. Gradient @ cell 1 of P would be: greenGaussGrad(mesh, P)[1, 1, :]
=#
function greenGaussGrad(mesh::Mesh, matrix::AbstractArray{Float64, 2}, valuesAtFaces=false)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    nVars = size(matrix, 2)

    #Interpolate values to faces if necessary
    if valuesAtFaces != true
        faceVals = linInterp_3D(mesh, matrix)
    else
        faceVals = matrix
    end

    # Create matrix to hold gradients
    grad = zeros(nCells, nVars, 3)

    # Integrate fluxes from each face
    @inbounds @fastmath for f in eachindex(mesh.faces)
        ownerCell = mesh.faces[f][1]
        neighbourCell = mesh.faces[f][2]

        for v in 1:nVars
            for d in 1:3
                # Every face has an owner
                grad[ownerCell, v, d] += mesh.fAVecs[f][d] * faceVals[f, v]

                # Boundary faces don't - could split into two loops
                if neighbourCell > -1
                    grad[neighbourCell, v, d] -= mesh.fAVecs[f][d] * faceVals[f, v]
                end
            end
        end
    end

    # Divide integral by cell volume to obtain gradients
    @inbounds @fastmath for c in 1:nCells
        for v in 1:nVars
            for d in 1:3
                grad[c,v,d] /= mesh.cVols[c]
            end
        end
    end

    return grad
end

####################### Face value interpolation ######################
#=
    Interpolates to all INTERIOR faces
    Arbitrary value matrix interpolation
    Example Input Matrix:
    Cell      x1       x2       x3
    Cell 1    x1_c1    x2_c1    x3_c1
    Cell 2    x1_c2    x2_c2    x3_c2
    ...
    Outputs a matrix of the following form
    Cell      x1       x2       x3
    Face 1    x1_f1    x2_f1    x3_f1
    Face 2    x1_f2    x2_f2    x3_f2
    ...
=#
function linInterp_3D(mesh::Mesh, matrix::AbstractArray{Float64, 2}, faceVals=zeros(2,2))
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    nVars = size(matrix, 2)

    # Make a new matrix if one is not passed in
    if faceVals == zeros(2,2)
        faceVals = zeros(nFaces, nVars)
    end

    # Boundary face fluxes must be set separately
    @inbounds @fastmath for f in 1:nFaces-nBdryFaces
        # Find value at face using linear interpolation
        c1 = mesh.faces[f][1]
        c2 = mesh.faces[f][2]

        #TODO: Precompute these distances?
        c1Dist = mag(mesh.cCenters[c1] .- mesh.fCenters[f])
        c2Dist = mag(mesh.cCenters[c2] .- mesh.fCenters[f])
        totalDist = c1Dist + c2Dist

        for v in 1:nVars
            faceVals[f, v] = matrix[c1, v]*(c2Dist/totalDist) + matrix[c2, v]*(c1Dist/totalDist)
        end
    end

    return faceVals
end

# Similar to linInterp_3D. Instead of linearly interpolating, selects the maximum value of the two adjacent cell centers as the face value
function maxInterp(mesh::Mesh, vars...)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

    nVars = size(vars, 1)
    faceVals = zeros(nFaces, nVars)

    # Boundary face fluxes must be set separately
    @inbounds @fastmath for f in 1:nFaces-nBdryFaces
        c1 = mesh.faces[f][1]
        c2 = mesh.faces[f][2]

        for v in 1:nVars
            faceVals[f, v] = max(vars[v][c1], vars[v][c2])
        end
    end

    return faceVals
end

# Similar to linInterp_3D. Instead of linearly interpolating, calculates the change (delta) of each variable across each face (required for JST method)
function faceDeltas(mesh::Mesh, sln::SolutionState)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    nVars = size(sln.cellState, 2)

    faceDeltas = zeros(nFaces, nVars)

    # Boundary face fluxes must be set separately (faceDelta is zero at all possible boundary conditions right now)
    @inbounds @fastmath for f in 1:nFaces-nBdryFaces
        ownerCell = mesh.faces[f][1]
        neighbourCell = mesh.faces[f][2]

        for v in 1:nVars
            faceDeltas[f, v] = sln.cellState[neighbourCell, v] - sln.cellState[ownerCell, v]
        end
    end

    return faceDeltas
end

# TODO: MUSCL Interp

######################### Convective Term Things #######################
# Calculates eps2 and eps4, the second and fourth-order artificial diffusion coefficients used in the JST method
function unstructured_JSTEps(mesh::Mesh, sln::SolutionState, k2=0.5, k4=(1/32), c4=1, gamma=1.4, R=287.05)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

    # Calc Pressure Gradient
    @views P = sln.cellPrimitives[:,1]
    P = reshape(P, nCells, :)
    gradP = greenGaussGrad(mesh, P, false)
    gradP = reshape(gradP, nCells, 3)

    sj = zeros(nCells) # 'Sensor' used to detect shock waves and apply second-order artificial diffusion to stabilize solution in their vicinity
    sjCount = zeros(nCells) # Store the number of sj's calculated for each cell, cell-center value will be the average of all of them
    @inbounds @fastmath for f in 1:nFaces-nBdryFaces
        # At each internal face, calculate sj, rj, eps2, eps4
        ownerCell = mesh.faces[f][1]
        neighbourCell = mesh.faces[f][2]
        d = mesh.cCenters[neighbourCell] .- mesh.cCenters[ownerCell]

        # 1. Find pressure at owner/neighbour cells
        oP = P[ownerCell]
        nP = P[neighbourCell]

        # 2. Calculate pressures at 'virtual' far-owner and far-neighbour cells using the pressure gradient (2nd-order)
        @views farOwnerP = nP - 2*dot(d, gradP[ownerCell, :])
        @views farNeighbourP = oP + 2*dot(d, gradP[neighbourCell, :])

        # 3. With the known and virtual values, can calculate sj at each cell center.
        sj[ownerCell] += (abs( nP - 2*oP + farOwnerP )/ max( abs(nP - oP) + abs(oP - farOwnerP), 0.0000000001))^2
        sjCount[ownerCell] += 1
        sj[neighbourCell] += (abs( oP - 2*nP + farNeighbourP )/ max( abs(farNeighbourP - nP) + abs(nP - oP), 0.0000000001))^2
        sjCount[neighbourCell] += 1
    end

    rj = zeros(nCells) # 'Spectral radius' -> maximum possible speed of wave propagation relative to mesh
    @inbounds @fastmath for c in 1:nCells
        @views rj[c] = mag(sln.cellPrimitives[c,3:5]) +  sqrt(gamma * R * sln.cellPrimitives[c,2]) # Velocity magnitude + speed of sound
        sj[c] /= sjCount[c] # Average the sj's computed at each face for each cell
    end

    # Values of rj, sj at faces is the maximum of their values at the two adjacent cell centers
    rjsjF = maxInterp(mesh, rj, sj) # column one is rj, column two is sj, both at face centers

    # Calculate eps2 and eps4
    eps2 = zeros(nFaces)
    eps4 = zeros(nFaces)
    for f in 1:nFaces-nBdryFaces
        eps2[f] = k2 * rjsjF[f,2] * rjsjF[f,1]
        eps4[f] = max(0, k4*rjsjF[f,1] - c4*eps2[f])
    end

    return eps2, eps4
end

#=
    Inputs: Expects that sln.cellState, sln.cellPrimitives and sln.cellFluxes are up-to-date
    Outputs: Updates sln.faceFluxes and sln.cellResiduals
    Returns: Updated sln.cellResiduals

    Applies classical JST method: central differencing + JST artificial diffusion. Each face treated as a 1D problem
    http://aero-comlab.stanford.edu/Papers/jst_2015_updated_07_03_2015.pdf -  see especially pg.5-6
=#
function unstructured_JSTFlux(mesh::Mesh, sln::SolutionState, boundaryConditions, gamma, R)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    nVars = size(sln.cellState, 2)

    #### 1. Centrally differenced fluxes ####
    linInterp_3D(mesh, sln.cellFluxes, sln.faceFluxes)

    #### 2. Add JST artificial Diffusion ####
    fDeltas = faceDeltas(mesh, sln)
    fDGrads = greenGaussGrad(mesh, fDeltas, false)
    eps2, eps4 = unstructured_JSTEps(mesh, sln, 0.5, (1.2/32), 1, gamma, R)

    diffusionFlux = zeros(nVars)
    @inbounds @fastmath for f in 1:nFaces-nBdryFaces
        ownerCell = mesh.faces[f][1]
        neighbourCell = mesh.faces[f][2]
        d = mesh.cCenters[neighbourCell] .- mesh.cCenters[ownerCell]

        @views fD = fDeltas[f,:]
        @views farOwnerfD = fD .- dot(d, fDGrads[ownerCell,:,:])
        @views farNeighbourfD = fD .+ dot(d, fDGrads[ownerCell,:,:])

        diffusionFlux = eps2[f]*fD - eps4[f]*(farNeighbourfD - 2*fD + farOwnerfD)

        # Add diffusion flux in component form
        unitFA = normalize(mesh.fAVecs[f])
        for v in 1:nVars
            i1 = (v-1)*3+1
            i2 = i1+2
            sln.faceFluxes[f,i1:i2] .-= (diffusionFlux[v] .* unitFA)
        end
    end

    #### 3. Apply boundary conditions ####
    for b in 1:nBoundaries
        bFunctionIndex = 2*b-1
        boundaryConditions[bFunctionIndex](mesh, sln, b, boundaryConditions[bFunctionIndex+1])
    end

    #### 4. Integrate fluxes at in/out of each cell (sln.faceFluxes) to get change in cell center values (sln.fluxResiduals) ####
    return integrateFluxes_unstructured3D(mesh, sln, boundaryConditions)
end

######################### TimeStepping #################################
# Calculate sln.cellPrimitives and sln.cellFluxes from sln.cellState
function decodeSolution_3D(sln::SolutionState, R=287.05, Cp=1005)
    nCells = size(sln.cellState, 1)
    for c in 1:nCells
        # Updates cell primitives
        @views decodePrimitives3D!(sln.cellPrimitives[c,:], sln.cellState[c,:], R, Cp)
        # Updates mass, xMom, eV2 x,y,z-direction fluxes
        @views calculateFluxes3D!(sln.cellFluxes[c, :], sln.cellPrimitives[c,:], sln.cellState[c,:])
    end
end

#=
    Calculates sln.fluxResiduals from sln.faceFluxes
=#
function integrateFluxes_unstructured3D(mesh::Mesh, sln::SolutionState, boundaryConditions)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

    # Recomputing flux balances, so wipe existing values
    fill!(sln.fluxResiduals, 0)
    nVars = size(sln.fluxResiduals, 2)

    #### Flux Integration ####
    @inbounds @fastmath for f in eachindex(mesh.faces)
        ownerCell = mesh.faces[f][1]
        neighbourCell = mesh.faces[f][2]

        for v in 1:nVars
            i1 = (v-1)*3 + 1
            i2 = i1+2
            @views flow = dot(sln.faceFluxes[f, i1:i2], mesh.fAVecs[f])

            # Subtract from owner cell
            sln.fluxResiduals[ownerCell, v] -= flow
            # Add to neighbour cell
            if neighbourCell > -1
                sln.fluxResiduals[neighbourCell, v] += flow
            end
        end
    end

    # Divide by cell volume
    for c in 1:nCells
        for v in 1:nVars
            sln.fluxResiduals[c,v] /= mesh.cVols[c]
        end
    end

    return sln.fluxResiduals
end

######################### Boundary Conditions #########################
# The job of the BC functions is to calculate/enforce face fluxes at boundary faces.
# The JST method only calculates fluxes at internal faces, then these conditions are applied to calculate them at the boundaries

# InletConditions: [ Static Pressure, Static Temperture, Ux, Uy, Uz, Cp ]
function supersonicInletBoundary!(mesh::Mesh, sln::SolutionState, boundaryNumber, inletConditions)
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
# TODO: Allow for variation of R and gamma
function subsonicInletBoundary!(mesh, sln::SolutionState, boundaryNumber, inletConditions)
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

function pressureOutletBoundary!(mesh, sln::SolutionState, boundaryNumber, outletPressure)
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

function zeroGradientBoundary!(mesh::Mesh, sln::SolutionState, boundaryNumber, _)
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

function wallBoundary!(mesh::Mesh, sln::SolutionState, boundaryNumber, _)
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
symmetryBoundary! = wallBoundary!

function emptyBoundary!(mesh::Mesh, sln::SolutionState, boundaryNumber, _)
    return
end

######################### Solvers #######################
#=
    This is where the magic happens!
    Do CFD!

    Arguments:
        mesh: see dataStructuresDefinitions.md / dataStructures.jl
        meshPath:           (string) path to folder where mesh is stored
        cellPrimitives:     initial cell-center conditions, see dataStructuresDefinitions.md / dataStructures.jl
        boundaryConditions: Array of alternating boundary condition function references and associated boundary condition parameters:
            Ex: [ emptyBoundary, [], supersonicInletBoundary, [P, T, Ux, Uy, Uz, Cp], zeroGradientBoundary, [], symmetryBoundary, [], zeroGradientBoundary, [], wallBoundary, [] ]
            Order of BC's must match the order of boundaries defined in the mesh's 'boundaries' file
        timeIntegrationFn:  The desired time integration function from timeDiscretizations.jl (should update sln.cellState, sln.cellPrimitives, and sln.cellFluxes based on their old values and sln.cellResiduals (obtained from fluxFunction))
        fluxFunction:       The desired function to calculate sln.cellResiduals from sln.cellState, sln.cellPrimitives, and sln.cellFluxes
        initDt:             Initial time step (s)
        endTime:            Simulation End Time (s). Start Time = 0
        outputInterval:     Writes solution/restart files every outputInterval seconds of simulation time
        targetCFL:          Target maximum CFL in the computation domain (time step will be adjusted based on this value)
        gamma/R/Cp:         Fluid properties
        silent:             Controls whether progress is written to console (can slow down simulation slightly for very simple cases)
        restart:            Controls whether restarting from a restart file
        createRestartFile:  Controls whether to write restart files. These are overwritten every time they are outputted.
        createVTKOutput:    Controls whether to write .vtk output. These are not overwritten every time they are outputted.
        restartFiles:       (string) path to restart file to read/write from

    Returns:
        sln.cellPrimitives at end of simulation

    Outputs:
        restart file (if specified)
        .vtk files (is specified)
=#
function unstructured3DFVM(mesh::Mesh, meshPath, cellPrimitives::Matrix{Float64}, boundaryConditions, timeIntegrationFn=forwardEuler,
        fluxFunction=unstructured_JSTFlux; initDt=0.001, endTime=0.14267, outputInterval=0.01, targetCFL=0.2, gamma=1.4, R=287.05, Cp=1005,
        silent=true, restart=false, createRestartFile=true, createVTKOutput=true, restartFile="JuliaCFDRestart.txt")

    if !silent
        println("Initializing Simulation")
    end
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    nDims = 3

    # Each dimension adds one momentum equation
    nVars = 2+nDims
    # Each dimension adds a flux for each conserved quantity
    nFluxes = nVars*nDims

    if restart
        if !silent
            println("Reading restart file: $restartFile")
        end

        cellPrimitives = readRestartFile(restartFile)

        # Check that restart file has the same nnumber of cells as the mesh
        nCellsRestart = size(cellPrimitives, 1)
        if nCellsRestart != nCells
            throw(ArgumentError("Number of cells in restart file ($nCellsRestart) does not match the present mesh ($nCells). Please provide a matching mesh and restart file."))
        end
    end

    # rho, xMom, total energy from P, T, Ux, Uy, Uz
    cellState = encodePrimitives3D(cellPrimitives, R, Cp)
    cellFluxes = zeros(nCells, nFluxes)
    fluxResiduals = zeros(nCells, nVars)
    faceFluxes = zeros(nFaces, nFluxes)
    # Initialize solution state
    sln = SolutionState(cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes)

    # Calculates cell fluxes, primitives from cell state
    decodeSolution_3D(sln, R, Cp)

    if !silent
        println("Starting iterations")
    end

    dt = initDt
    if timeIntegrationFn==LTSEuler
        dt = zeros(nCells)
    end

    currTime = 0
    timeStepCounter = 0
    nextOutputTime = outputInterval
    writeOutputThisIteration = false
    CFLvec = zeros(nCells)
    while currTime < endTime
        ############## Timestep adjustment #############
        if timeIntegrationFn == LTSEuler
            CFL = 1.0
            if timeStepCounter < 10
                # Ramp up CFL linearly in the first ten time steps to reach the target CFL
                CFL = (timeStepCounter+1) * targetCFL / 10
            else
                CFL = targetCFL
            end

            # Store the target CFL for the present time step in the first element of the dt vector.
            # The time discretization function will determine the actual local time step based on this target CFL
            dt[1] = CFL
        else
            CFL!(CFLvec, mesh, sln, dt, gamma, R)
            CFL = maximum(CFLvec)

            # If CFL too high, attempt to preserve stability by cutting timestep size in half
            if CFL > targetCFL*1.01
                dt *= targetCFL/(2*CFL)
            # Otherwise slowly approach target CFL
            else
                dt *= ((targetCFL/CFL - 1)/10+1)
            end

            # Cut timestep to hit endtime/writeTimes as necessary
            if (endTime - currTime) < dt
                dt = endTime - currTime
            elseif (nextOutputTime - currTime) < dt
                dt = nextOutputTime - currTime
                writeOutputThisIteration = true
            end
        end

        ############## Take a timestep #############
        sln = timeIntegrationFn(mesh, fluxFunction, sln, boundaryConditions, gamma, R, Cp, dt)

        if timeIntegrationFn == LTSEuler
            currTime += CFL
            if (nextOutputTime - currTime) < CFL
                writeOutputThisIteration = true
            end
        else
            currTime += dt
        end
        timeStepCounter += 1

        if !silent
            @printf("Timestep: %5.0f, simTime: %9.4g, Max CFL: %9.4g \n", timeStepCounter, currTime, CFL)
        end

        if writeOutputThisIteration
            updateSolutionOutput(sln.cellPrimitives, restartFile, meshPath, createRestartFile, createVTKOutput)
            writeOutputThisIteration = false
            nextOutputTime += outputInterval
        end
    end

    updateSolutionOutput(sln.cellPrimitives, restartFile, meshPath, createRestartFile, createVTKOutput)
    return sln.cellPrimitives
end
