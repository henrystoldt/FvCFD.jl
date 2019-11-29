# finiteVolumeRoe.jl

include("constitutiveRelations.jl")
include("vectorFunctions.jl")

######################### CFL ########################
# TODO: Generalize cell size
# TODO: 3D definition of CFL: Sum up in all directions
function CFL(U, T, dt, dx, gamma=1.4, R=287.05)
    return (abs(U[1]) + sqrt(gamma * R * T)) * dt / dx
end


# Handles scalar or vector-valued variables
function linInterp(mesh, values...)    
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nFaces = size(faces, 1)
    nBdryFaces = size(boundaryFaces, 1)
    
    fVals = Array{Float64, 1}(undef, nFaces)

    result = []
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

        ########### Do Interpolations ###########
        for i in 1:nFaces
            if i > nFaces-nBdryFaces
                # Faces on boundaries not treated, must be set externally
                push!(fVals, defaultVal)
                
            else
                # Find value at face using linear interpolation
                c1 = faces[i][1]
                c2 = faces[i][2]

                #TODO: Precompute these distances
                c1Dist = mag(cCenters[c1] .- fCenters[i])
                c2Dist = mag(cCenters[c2] .- fCenters[i])
                totalDist = c1Dist + c2Dist
                push!(fVals, vals[c1].*(c2Dist/totalDist) .+ vals[c2].*(c1Dist/totalDist))
            end
        end

        push!(result, fVals)
    end
    return result
end


function musclDifference(mesh, rho, xMom, eV2, dx)
    #= Uses MUSCL differencing to get the conserved values at the left and right side of each face
    Then determines the primitive values at the left and right side. BC's are treated as zero grad for now
    =#
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nFaces = size(faces, 1)
    nBdryFaces = size(boundaryFaces, 1)
    
    consLeft = Array{Float64, 2}(undef, nFaces, 3)
    consRight = Array{Float64, 2}(undef, nFaces, 3)

    for i in 1:(nFaces-nBdryFaces)
        #L/R Conserved quantities on non-boundary faces
        #This treatment doesn't deal with connectivity
        #How do we deal with the 2nd and 2ndlast faces?? Neither are on boundaries, but still can't be differenced
        

        #rho
        consLeft[i][1] = leftFaceVals(dx)
        consRight[i][1] = rightFaceVals(dx)

        #xMom
        consLeft[i][2] = leftFaceVals(dx)
        consRight[i][2] = rightFaceVals(dx)

        #eV2
        consLeft[i][3] = leftFaceVals(dx)
        consRight[i][3] = rightFaceVals(dx)



    end

    for i in 1:nBdryFaces
        #Treatment of boundary faces

    end

    primsLeft = Array{Float64, 2}(undef, nFaces, 3)
    primsRight = Array{Float64, 2}(undef, nFaces, 3)

    for i in 1:nFaces
        #Find L/R values of (rho, u, H, P) for each face 


    end


    return consLeft, primsLeft, consRight, primsRight
end

function leftFaceVals()


end

function rightFaceVals()
end

function vanAlbeda()

end


######################### Convective Term Things #######################



######################### Solvers #######################
function upwindFVMRoe1D(mesh, P, T, U; initDt=0.001, endTime=0.14267, targetCFL=0.2, gamma=1.4, R=287.05, Cp=1005, debug=false)
    ######### MESH ############
    # Extract mesh into local variables for readability
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells = size(cells, 1)
    nFaces = size(faces, 1)
    nBdryFaces = size(boundaryFaces, 1)
    bdryFaceIndices = Array(nFaces-nBdryFaces:nFaces)

    ########### Variable Arrays #############
    # State variables, values are the averages/cell center values
    rho = Array{Float64, 1}(undef, nCells)
    xMom = Array{Float64, 1}(undef, nCells)
    eV2 = Array{Float64, 1}(undef, nCells)
    stateVars = [ rho, xMom, eV2 ]
    # Flux variables for equations 2 and 3 (xMom is flux variable for the continuity eqn)
    rhoU2p = Array{Float64, 1}(undef, nCells)
    rhoUeV2PU = Array{Float64, 1}(undef, nCells)

    # Calc state and flux variables from primitives
    #TODO: Multi-D
    for i in 1:nCells
        rho[i], xMom[i], eV2[i], rhoU2p[i], rhoUeV2PU[i] = encodePrimitives3D(P[i], T[i], U[i])
    end
    
    if debug
        println("Rho: $rho")
        println("xMom: $xMom")
        println("eV2: $eV2")
    end

    ########### SOLVER ###########
    dt = initDt
    currTime = 0
    while currTime < endTime

        if (endTime - currTime) < dt
            dt = endTime - currTime
        end
        
        # Calculate fluxes through each face
        # TODO: y and z momemtum-fluxes + equations
        # xMassFlux, xMomFlux, xeV2Flux = upwindInterp(mesh, U, xMom, rhoU2p, rhoUeV2PU)
        #xMassFlux, xMomFlux, xeV2Flux, faceP = linInterp(mesh, xMom, rhoU2p, rhoUeV2PU, P)

        #Use MUSCL differencing to get values at left/right of each face. Boundaries are treated inside this function
        consLeft, primsLeft, consRight, primsRight = musclDifference(mesh, rho, xMom, eV2, dx)

        #Find ROE averaged quantities at each face
        rhoRoe, uRoe, hRoe, aRoe = roeAveraged(rhoLeft,rhoRight,uLeft,uRight,hLeft,hRight)

        #Use left, right and roe averaged quantities to create the flux matrices at every face. Entropy fix happens inside this step
        xMassFlux, xMomFlux, xeV2Flux = fluxVectors(left, right, Roe, faceVector )

        
  
        fluxVars = [ xMassFlux, xMomFlux, xeV2Flux ]
        

        if debug
            println("xMass Flux: $xMassFlux")
            println("xMom Flux: $xMomFlux")
            println("xEv2 Flux: $xeV2Flux")
        end

        # Use fluxes to update values in each cell
        for i in 1:(nFaces-nBdryFaces)
            fA = mag(fAVecs[i]) #Face Area

            ownerCell = faces[i][1]
            neighbourCell = faces[i][2]
            for v in 1:3
                #Not sure if I need the +/-?? I think that's accounted for in the dot area vector? Or maybe not?
                #Also don't think you need to divide each face by cell vol, only the final result?
                stateVars[v][ownerCell] -= fluxVars[v][i]*fA*dt/cVols[ownerCell]
                stateVars[v][neighbourCell] += fluxVars[v][i]*fA*dt/cVols[neighbourCell]
            end
        end

        
        if debug
            println("Rho2: $rho")
            println("xMom2: $xMom")
            println("eV22: $eV2")
        end

        #Unpack state variables
        #Not sure if this step is required??
        rho = stateVars[1][:]
        xMom = stateVars[2][:]
        eV2 = stateVars[3][:]
        
        # Boundaries
        ############### Boundaries ################
        # Waves never reach the boundaries, so boundary treatment doesn't need to be good
        
        #= Farfield neumann-BC is dealt with during the MUSCL-difference I hope
        copyValues(3, 2, stateVars)
        copyValues(2, 1, stateVars)
        copyValues(nCells-2, nCells-1, stateVars)
        copyValues(nCells-1, nCells, stateVars)    
        =#

        # Decode primitive values
        for i in 1:nCells
            P[i], T[i], U[i], rhoU2p[i], rhoUeV2PU[i] = decodePrimitives3D(rho[i], xMom[i], eV2[i])
        end

        currTime += dt
        
        ############## CFL Calculation, timestep adjustment #############
        maxCFL = 0
        for i in 1:nCells
            dx = 1/nCells
            maxCFL = max(maxCFL, CFL(U[i], T[i], dt, dx, gamma, R))
        end
        # Adjust time step to approach target CFL
        dt *= ((targetCFL/maxCFL - 1)/5+1)

    end

    return P, U, T, rho
end

#TODO: Proper boundary treatment