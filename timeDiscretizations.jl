function forwardEuler(mesh, fluxResidualFn, solutionState, boundaryConditions, dt)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState

    fluxResiduals = fluxResidualFn(mesh, solutionState, boundaryConditions)
    cellState .+= fluxResiduals.*dt
    decodeSolution(solutionState)

    return solutionState
end

function RK2_Mid(mesh, fluxResidualFn, solutionState, boundaryConditions, dt)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState

    fluxResiduals1 = fluxResidualFn(mesh, solutionState, boundaryConditions)
    halfwayEstimate = cellState .+ fluxResiduals1.*dt/2
    solutionState2 = [ halfwayEstimate, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]
    decodeSolution(solutionState2)

    fluxResiduals2 = fluxResidualFn(mesh, solutionState2, boundaryConditions)
    cellState .+= fluxResiduals2.*dt
    decodeSolution(solutionState)

    return solutionState
end

function RK4(mesh, fluxResidualFn, solutionState, boundaryConditions, dt)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState

    fluxResiduals1 = fluxResidualFn(mesh, solutionState, boundaryConditions)
    halfwayEstimate = cellState .+ fluxResiduals1*dt/2
    lastSolutionState = [ halfwayEstimate, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]
    decodeSolution(lastSolutionState)

    fluxResiduals2 = fluxResidualFn(mesh, lastSolutionState, boundaryConditions)
    halfwayEstimate2 = cellState .+ fluxResiduals2*dt/2
    lastSolutionState = [ halfwayEstimate2, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]
    decodeSolution(lastSolutionState)

    fluxResiduals3 = fluxResidualFn(mesh, lastSolutionState, boundaryConditions)
    finalEstimate1 = cellState .+ fluxResiduals3*dt
    lastSolutionState = [ finalEstimate1, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]
    decodeSolution(lastSolutionState)

    fluxResiduals4 = fluxResidualFn(mesh, lastSolutionState, boundaryConditions)
    cellState .+= (fluxResiduals1 .+ 2*fluxResiduals2 .+ 2*fluxResiduals3 .+ fluxResiduals4 )*(dt/6)
    decodeSolution(solutionState)

    return solutionState
end

function ShuOsher(mesh, fluxResidualFn, solutionState, boundaryConditions, dt)
    cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes = solutionState

    fluxResiduals1 = fluxResidualFn(mesh, solutionState, boundaryConditions)
    endEstimate = cellState .+ fluxResiduals1.*dt
    lastSolutionState = [ endEstimate, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]
    decodeSolution(lastSolutionState)

    fluxResiduals2 = fluxResidualFn(mesh, lastSolutionState, boundaryConditions)
    estimate2 = (3/4).*cellState .+ (1/4).*(endEstimate .+ fluxResiduals2.*dt)
    lastSolutionState = [ estimate2, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]
    decodeSolution(lastSolutionState)

    fluxResiduals3 = fluxResidualFn(mesh, lastSolutionState, boundaryConditions)
    cellState .= (1/3).*cellState .+ (2/3).*(estimate2 .+ dt.*fluxResiduals3)
    solutionState = [ cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes ]
    decodeSolution(solutionState)

    return solutionState
end
