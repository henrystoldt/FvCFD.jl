function forwardEuler(mesh, fluxResidualFn, solutionState, boundaryConditions, gamma, R, Cp, dt)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState

    fluxResiduals = fluxResidualFn(mesh, solutionState, boundaryConditions, gamma, R)
    @fastmath cellState .+= fluxResiduals.*dt
    @fastmath decodeSolution_3D(solutionState, R, Cp)

    return solutionState
end

function LTSEuler(mesh, fluxResidualFn, solutionState, boundaryConditions, gamma, R, Cp, dt)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState
    targetCFL = dt[1]

    fluxResiduals = fluxResidualFn(mesh, solutionState, boundaryConditions, gamma, R)

    OpenFOAMCFL!(dt, mesh, solutionState, 1, gamma, R)
    dt .= targetCFL ./ dt
    smoothTimeStep!(dt, mesh, 0.1)
    smoothTimeStep!(dt, mesh, 0.1)
    cellState .+= fluxResiduals .* dt
    decodeSolution_3D(solutionState, R, Cp)

    return solutionState
end

function smoothTimeStep!(dt, mesh::Mesh, diffusionCoefficient=0.2)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

    timeFluxes = zeros(nCells)
    surfaceAreas = zeros(nCells)
    for f in 1:nFaces-nBdryFaces
        ownerCell = mesh.faces[f][1]
        neighbourCell = mesh.faces[f][2]
        timeFlux = (dt[ownerCell] - dt[neighbourCell]) * mag(mesh.fAVecs[f])
        surfaceAreas[ownerCell] += mag(mesh.fAVecs[f])
        surfaceAreas[neighbourCell] += mag(mesh.fAVecs[f])
        timeFluxes[ownerCell] -= timeFlux
        timeFluxes[neighbourCell] += timeFlux
    end

    timeFluxes .*= (diffusionCoefficient ./ surfaceAreas)

    for i in eachindex(timeFluxes)
        timeFluxes[i] = min(0, timeFluxes[i])
    end

    dt .+= timeFluxes
end

function RK2_Mid(mesh, fluxResidualFn, solutionState, boundaryConditions, gamma, R, Cp, dt)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState

    fluxResiduals1 = fluxResidualFn(mesh, solutionState, boundaryConditions, gamma, R )
    halfwayEstimate = cellState .+ fluxResiduals1.*dt/2
    solutionState2 = [ halfwayEstimate, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]
    decodeSolution_3D(solutionState2, R, Cp)

    fluxResiduals2 = fluxResidualFn(mesh, solutionState2, boundaryConditions, gamma, R)
    cellState .+= fluxResiduals2.*dt
    decodeSolution_3D(solutionState, R, Cp)

    return solutionState
end

function RK4(mesh, fluxResidualFn, solutionState, boundaryConditions, gamma, R, Cp, dt)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState

    fluxResiduals1 = fluxResidualFn(mesh, solutionState, boundaryConditions, gamma, R)
    halfwayEstimate = cellState .+ fluxResiduals1*dt/2
    lastSolutionState = [ halfwayEstimate, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]
    decodeSolution_3D(lastSolutionState, R, Cp)

    fluxResiduals2 = fluxResidualFn(mesh, lastSolutionState, boundaryConditions, gamma, R)
    halfwayEstimate2 = cellState .+ fluxResiduals2*dt/2
    lastSolutionState = [ halfwayEstimate2, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]
    decodeSolution_3D(lastSolutionState, R, Cp)

    fluxResiduals3 = fluxResidualFn(mesh, lastSolutionState, boundaryConditions, gamma, R)
    finalEstimate1 = cellState .+ fluxResiduals3*dt
    lastSolutionState = [ finalEstimate1, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]
    decodeSolution_3D(lastSolutionState, R, Cp)

    fluxResiduals4 = fluxResidualFn(mesh, lastSolutionState, boundaryConditions, gamma, R)
    cellState .+= (fluxResiduals1 .+ 2*fluxResiduals2 .+ 2*fluxResiduals3 .+ fluxResiduals4 )*(dt/6)
    decodeSolution_3D(solutionState, R, Cp)

    return solutionState
end

function ShuOsher(mesh, fluxResidualFn, solutionState, boundaryConditions, gamma, R, Cp, dt)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState

    fluxResiduals1 = fluxResidualFn(mesh, solutionState, boundaryConditions, gamma, R)
    endEstimate = cellState .+ fluxResiduals1.*dt
    lastSolutionState = [ endEstimate, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]
    decodeSolution_3D(lastSolutionState, R, Cp)

    fluxResiduals2 = fluxResidualFn(mesh, lastSolutionState, boundaryConditions, gamma, R)
    estimate2 = (3/4).*cellState .+ (1/4).*(endEstimate .+ fluxResiduals2.*dt)
    lastSolutionState = [ estimate2, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]
    decodeSolution_3D(lastSolutionState, R, Cp)

    fluxResiduals3 = fluxResidualFn(mesh, lastSolutionState, boundaryConditions, gamma, R)
    cellState .= (1/3).*cellState .+ (2/3).*(estimate2 .+ dt.*fluxResiduals3)
    solutionState = [ cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]
    decodeSolution_3D(solutionState, R, Cp)

    return solutionState
end

#TODO: For implicit methods, need to compute the flux Jacobians at each edge, instead of just the fluxes
    # Use Jacobians as coefficients in matrix representing timestepping equations
    # Then solve with GMRES or some other matrix solver
