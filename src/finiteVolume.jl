using Printf
include("constitutiveRelations.jl")
include("vectorFunctions.jl")
include("timeDiscretizations.jl")
include("mesh.jl")
include("output.jl")
include("dataStructures.jl")
include("boundaryConditions.jl")
include("numerics.jl")
include("JST.jl")

__precompile__()

mutable struct SolverStatus
    currentTime::Float64
    nTimeSteps::Int64
    nextOutputTime::Float64
    endTime::Float64
end

######################### Initialization ###########################
# Returns cellPrimitives matrix for uniform solution
function initializeUniformSolution3D(mesh, P, T, Ux, Uy=0, Uz=0, nDims=3)
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

    nVars = 2 + nDims
    cellPrimitives = Matrix{Float64}(undef, nCells, nVars)
    initialValues = [ P, T, Ux, Uy, Uz ][1:nVars]
    
    for c in 1:nCells0
        cellPrimitives[c, :] = initialValues
    end

    return cellPrimitives
end

######################### Taking and adjusting time step sizes #######################
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
            println("Warning: Negative temperature at face $f")
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

function populateSolution(cellPrimitives, nCells, nFaces, R, Cp, nDims=3)
    # Each dimension adds one momentum equation
    nConservedVars = 2+nDims
    # Each dimension adds a flux for each conserved quantity
    nFluxes = nConservedVars*nDims

    # rho, xMom, total energy from P, T, Ux, Uy, Uz
    cellState = encodePrimitives3D(cellPrimitives, R, Cp)
    cellFluxes = zeros(nCells, nFluxes)
    fluxResiduals = zeros(nCells, nConservedVars)
    faceFluxes = zeros(nFaces, nFluxes)

    # Initialize solution state
    sln = SolutionState(cellState, cellFluxes, cellPrimitives, fluxResiduals, faceFluxes)

    # Calculates cell fluxes, primitives from cell state
    decodeSolution_3D(sln, R, Cp)

    return sln
end

function restrictTimeStep(status, desiredDt)
    maxStep = min(status.endTime-status.currentTime, status.nextOutputTime-status.currentTime)

    if desiredDt > maxStep
        return true, maxStep
    else
        return false, desiredDt
    end
end

function adjustTimeStep_LTS(targetCFL, dt, status::SolverStatus)
    CFL = 1.0
    if status.nTimeSteps < 10
        # Ramp up CFL linearly in the first ten time steps to reach the target CFL
        CFL = (status.nTimeSteps+1) * targetCFL / 10
    else
        CFL = targetCFL
    end

    writeOutputThisIteration, CFL = restrictTimeStep(status, CFL)

    # Store the target CFL for the present time step in the first element of the dt vector.
    # The time discretization function will determine the actual local time step based on this target CFL
    dt[1] = CFL # TODO: Cleaner way to pass this information

    return writeOutputThisIteration, dt, CFL
end
        
function adjustTimeStep(maxCFL, targetCFL, dt, status)
    # If CFL too high, attempt to preserve stability by cutting timestep size in half
    if maxCFL > targetCFL*1.01
        dt *= targetCFL/(2*maxCFL)
    # Otherwise slowly approach target CFL
    else
        dt *= ((targetCFL/maxCFL - 1)/10+1)
    end

    writeOutputThisIteration, dt = restrictTimeStep(status, dt)
    return writeOutputThisIteration, dt, maxCFL
end

function advanceStatus!(status::SolverStatus, dt, CFL, timeIntegrationFn, silent)
    status.nTimeSteps += 1

    if timeIntegrationFn == LTSEuler
        status.currentTime += CFL
    else
        status.currentTime += dt
    end        

    if !silent
        @printf("Timestep: %5.0f, simTime: %9.4g, Max CFL: %9.4g \n", status.nTimeSteps, status.currentTime, CFL)
    end
end

######################### Main #######################
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
function solve(mesh::Mesh, meshPath, cellPrimitives::Matrix{Float64}, boundaryConditions, timeIntegrationFn=forwardEuler,
        fluxFunction=unstructured_JSTFlux; initDt=0.001, endTime=0.14267, outputInterval=0.01, targetCFL=0.2, gamma=1.4, R=287.05, Cp=1005,
        silent=true, restart=false, createRestartFile=true, createVTKOutput=true, restartFile="JuliaCFDRestart.txt")

    if !silent
        println("Initializing Simulation")
    end

    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)

    if !silent
        println("Mesh: $meshPath")
        println("Cells: $nCells")
        println("Faces: $nFaces")
        println("Boundaries: $nBoundaries")
    end

    if restart
        if !silent
            println("Reading restart file: $restartFile")
        end

        cellPrimitives = readRestartFile(restartFile)

        # Check that restart file has the same number of cells as the mesh
        nCellsRestart = size(cellPrimitives, 1)
        if nCellsRestart != nCells
            throw(ArgumentError("Number of cells in restart file ($nCellsRestart) does not match the present mesh ($nCells). Please provide a matching mesh and restart file."))
        end
    end

    sln = populateSolution(cellPrimitives, nCells, nFaces, R, Cp, 3)

    dt = initDt
    if timeIntegrationFn==LTSEuler
        # If using local time stepping, the time step will be different for each cell
        dt = zeros(nCells)
    end

    status = SolverStatus(0, 0, outputInterval, endTime)
    writeOutputThisIteration = false
    CFLvec = zeros(nCells)

    if !silent
        println("Starting iterations")
    end

    ### Main Loop ###
    while status.currentTime < status.endTime
        ############## Adjust timestep size #############        
        if timeIntegrationFn == LTSEuler # LTS = Local time stepping
            writeOutputThisIteration, dt, CFL = adjustTimeStep_LTS(targetCFL, dt, status)
        else
            CFL!(CFLvec, mesh, sln, dt, gamma, R)
            writeOutputThisIteration, dt, CFL = adjustTimeStep(maximum(CFLvec), targetCFL, dt, status)
        end

        ############## Take a timestep #############
        sln = timeIntegrationFn(mesh, fluxFunction, sln, boundaryConditions, gamma, R, Cp, dt)
        advanceStatus!(status, dt, CFL, timeIntegrationFn, silent)
        
        ############## Write Output #############
        if writeOutputThisIteration
            writeOutput(sln.cellPrimitives, restartFile, meshPath, createRestartFile, createVTKOutput)
            status.nextOutputTime += outputInterval
        end
    end
    
    # Always create output upon exit
    writeOutput(sln.cellPrimitives, restartFile, meshPath, createRestartFile, createVTKOutput)

    # Return current cell-center properties
    return sln.cellPrimitives
end
