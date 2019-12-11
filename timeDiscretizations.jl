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

    CFL!(dt, mesh, solutionState, 1, gamma, R)
    dt .= targetCFL ./ dt
    fluxResiduals = fluxResidualFn(mesh, solutionState, boundaryConditions, gamma, R)
    cellState .+= fluxResiduals .* dt
    decodeSolution_3D(solutionState, R, Cp)

    return solutionState
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
