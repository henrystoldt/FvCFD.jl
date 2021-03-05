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
        parameters = boundaryConditions[bFunctionIndex+1]
        boundaryConditions[bFunctionIndex](mesh, sln, b, parameters)
    end

    #### 4. Integrate fluxes at in/out of each cell (sln.faceFluxes) to get change in cell center values (sln.fluxResiduals) ####
    return integrateFluxes_unstructured3D(mesh, sln, boundaryConditions)
end

