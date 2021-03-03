# finiteVolumeRoe.jl

include("constitutiveRelations.jl")
include("vectorFunctions.jl")

######################### CFL ########################
function CFL(U, T, dt, dx, gamma=1.4, R=287.05)
    #println(T)
    return (abs(U[1]) + sqrt(gamma * R * T)) * dt / dx
end


# Handles scalar or vector-valued variables

function musclDifference(mesh, rho, xMom, eV2, dx, gamma=1.4; debug=false)
    #= Uses MUSCL differencing to get the conserved values at the left and right side of each face
    Then determines the primitive values at the left and right side. BC's are treated as zero grad for now
    =#
    #TODO: Rewrite all of these with parent/neighbor/parent+-1 indexing
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nFaces = size(faces, 1)
    nBdryFaces = size(boundaryFaces, 1)

    consLeft = Array{Float64, 2}(undef, nFaces, 3)
    consRight = Array{Float64, 2}(undef, nFaces, 3)

    #consLeft = []

    #consLeft[1,1] = 52

    #println("RHO: ", consLeft[1,1])
    #println(typeof(consLeft), typeof(xMom), typeof(eV2), typeof(dx))

    for i in 1:(nFaces-nBdryFaces)
        #L/R Conserved quantities on non-boundary faces
        #This treatment doesn't deal with connectivity
        #How do we deal with the 2nd and 2ndlast faces?? Neither are on boundaries, but still can't be differenced
        if i == 1
            #Deal with it as if it were next face for the backwards differences
            consLeft[i,1] = leftFaceVals(rho[3], rho[2], rho[1], dx) #These are fake news numbers until we treat the boundaries properly
            consRight[i,1] = rightFaceVals(rho[i+2], rho[i+1], rho[i], dx)

            consLeft[i,2] = leftFaceVals(xMom[3], xMom[2], xMom[1],dx)
            consRight[i,2] = rightFaceVals(xMom[i+2], xMom[i+1], xMom[i],dx)

            consLeft[i,3] = leftFaceVals(eV2[3], eV2[2], eV2[1], dx)
            consRight[i,3] = rightFaceVals(eV2[i+2], eV2[i+1], eV2[i], dx)

        elseif i==(nFaces-nBdryFaces)
            #TODO: Verify it actually trips this block for the last difference
            #Deal with it as if it were previous face
            consLeft[i,1] = leftFaceVals(rho[i+1], rho[i], rho[i-1], dx)
            consRight[i,1] = rightFaceVals(rho[nFaces-nBdryFaces+1], rho[nFaces-(nBdryFaces)], rho[nFaces-(nBdryFaces-1)], dx)

            consLeft[i,2] = leftFaceVals(xMom[i+1], xMom[i], xMom[i-1],dx)
            consRight[i,2] = rightFaceVals(xMom[nFaces-nBdryFaces+1], xMom[nFaces-(nBdryFaces)], xMom[nFaces-(nBdryFaces-1)],dx)

            consLeft[i,3] = leftFaceVals(eV2[i+1], eV2[i], eV2[i-1], dx)
            consRight[i,3] = rightFaceVals(eV2[nFaces-nBdryFaces+1], eV2[nFaces-(nBdryFaces)], eV2[nFaces-(nBdryFaces-1)], dx)

        else
            #rho
            consLeft[i,1] = leftFaceVals(rho[i+1], rho[i], rho[i-1], dx)
            consRight[i,1] = rightFaceVals(rho[i+2], rho[i+1], rho[i], dx)

            #xMom
            consLeft[i,2] = leftFaceVals(xMom[i+1], xMom[i], xMom[i-1],dx)
            consRight[i,2] = rightFaceVals(xMom[i+2], xMom[i+1], xMom[i],dx)

            #eV2
            consLeft[i,3] = leftFaceVals(eV2[i+1], eV2[i], eV2[i-1], dx)
            consRight[i,3] = rightFaceVals(eV2[i+2], eV2[i+1], eV2[i], dx)

        end

    end

    if debug
        println("Rho Left: ", consLeft[:,1])
        println("xMom Left: ", consLeft[:,2])
        println("E Left: ", consLeft[:,3])
    end


    for i in 1:nBdryFaces
        #Treatment of boundary faces
        #Determine what cell it's connected to, then copy values from that
        neighbourCell = faces[i+nFaces-nBdryFaces][2]
        if neighbourCell==-1 #TODO: Fix mesh treatment to not toss -1
            neighbourCell = nFaces-nBdryFaces
        end
        for j in 1:3
            consLeft[i+nFaces-nBdryFaces,j] = consLeft[neighbourCell,j]
            consRight[i+nFaces-nBdryFaces,j] = consRight[neighbourCell,j]
        end
        #And that's why if the world was zero gradient we would definitely be living in a simulation

    end

    tracker = consLeft - consRight

    primsLeft = Array{Float64, 2}(undef, nFaces, 3) #u, P, H
    primsRight = Array{Float64, 2}(undef, nFaces, 3)

    for i in 1:nFaces
        #Find L/R values of (rho, u, H, P) for each face

        primsLeft[i,1] = consLeft[i,2] / consLeft[i,1] #u =rhou/rho
        primsLeft[i,2] = (gamma-1)*( consLeft[i,3] - 0.5*consLeft[i,1]*primsLeft[i,1]^2 ) # P =(gamma-1)*(E-0.5rho(v^2))
        primsLeft[i,3] = ((consLeft[i,3]/consLeft[i,1])-0.5*primsLeft[i,1]^2  +  primsLeft[i,2])/consLeft[i,1] #H = (e+P)/rho = ( (E/rho -0.5v^2) + P )/rho

        primsRight[i,1] = consRight[i,2] / consRight[i,1]
        primsRight[i,2] = (gamma-1)*( consRight[i,3] - 0.5*consRight[i,1]*primsRight[i,1]^2 )
        primsRight[i,3] = ((consRight[i,3]/consRight[i,1])-0.5*primsRight[i,1]^2  +  primsRight[i,2])/consRight[i,1]

    end

    return consLeft, primsLeft, consRight, primsRight
end

function leftFaceVals(q_ip, q_i, q_im, dx; sig=-1)
    #=MUSCL difference interpolates values to the left face
    q_i - value at your cell
    q_ip - value at the cell opposite the face you're interpolating TODO
    q_im - value at the cell further away from the current face
    =#
    s = vanAlbeda(q_ip, q_i, q_im, dx)
    #s = vanLeer(q_ip, q_i, q_im, dx=dx)
    Q_L = q_i + 0.25*s*( (1-(sig)*s)*(q_i-q_im) + (1+(sig)*s)*(q_ip-q_i) ) * q_i

    return Q_L
end

function rightFaceVals(q_ip, q_i, q_im, dx; sig=-1)
    #=MUSCL difference to the right face
    Slightly different than left face function
        q_i - value at the cell
        q_im - value at the cell on the other side of your face
        q_ip - value at the cell further away from your current face
    =#
    s = vanAlbeda(q_ip, q_i, q_im, dx, direc=1)
    #s = vanLeer(q_ip, q_i, q_im, direc=1, dx=dx)  #1/vanLeer because it should be fwd/backwards
    Q_R = q_i - 0.25*s*( (1-(sig)*s)*(q_ip-q_i) + (1+(sig)*s)*(q_i-q_im) )*q_i

    return Q_R
end

function vanAlbeda(q_ip, q_i, q_im, dx; direc=0)
    #=Impplements the vanAlbeda flux limiter, with a delta value of dx^3
    =#
    delta = dx.*dx.*dx
    #delta = dx
    #delta = 1.0
    if direc ==0
        r = (q_i - q_im)/(q_ip - q_i + delta)
    else
        #r = (q_ip - q_i)/(q_i - q_im + delta)
        r = (q_i - q_im)/(q_ip - q_i + delta)
    end

    if r > 0
        #s = 2*r / (r^2 + 1)
        s = (r^2 + r)/(r^2 + 1)
    else
        s = 0
    end

    #s = (2 * (q_ip-q_i)*(q_i-q_im) + delta )/( (q_ip-q_i)^2 * (q_i-q_im)^2 + delta )

    return s
end

function vanLeer(q_ip, q_i, q_im; direc=0, dx=0.05)
    #=Van Leer slope limiter
    No Delta value required
    =#
    if direc ==0
        r = (q_i - q_im)/(q_ip - q_i + dx)
    else
        r = (q_ip-q_i)/(q_i-q_im + dx)
        #r = (q_i - q_im)/(q_ip - q_i + dx)
    end

    if r > 5
        r = 0.7
        #println("Did some weird stuff with limiter")
    end

    s = 0

    if r > 0
        s = 2*r/(r+1)
    end

    return s
end

function roeAveraged(n, rhoLeft,rhoRight,uLeft,uRight,HLeft,HRight, gamma=1.4)
    #=Takes inputs of left and right primitive quantities, returns roe averaged quantities for each
    =#

    roeRho = Array{Float64, 1}(undef, n)
    roeU = Array{Float64, 1}(undef, n)
    roeH = Array{Float64, 1}(undef, n)
    roeA = Array{Float64, 1}(undef, n)

    roeRho = sqrt.(rhoLeft .* rhoRight)
    roeU = ( sqrt.(rhoLeft).*uLeft + sqrt.(rhoRight).*uRight ) ./ (sqrt.(rhoLeft)+sqrt.(rhoRight))
    roeH = ( sqrt.(rhoLeft).*HLeft + sqrt.(rhoRight).*HRight ) ./ (sqrt.(rhoLeft)+sqrt.(rhoRight))

    println(minimum(roeH))

    #println("\n")
    println(minimum(roeU))


    roeA = sqrt.( (gamma-1).*(roeH - 0.5 .*roeU.^2) )

    return roeRho, roeU, roeH, roeA
end

function fluxVector(cons, prim, fAVec)
    #=This vector takes inputs of conserved and primitive quantities on the L or R side of a face, and returns the flux vector for that side
    =#
    nFaces = size(fAVec, 1)

    fluxVec = Array{Float64, 2}(undef, nFaces, 3) #Over all faces, storing rho*u, rhou^2+P, u(E+P)

    #In 2D, would need to multiply each element by the x/y unit vector components of the face, but in 1D that should already be dealt with??
    for i in 1:nFaces
        fluxVec[i,1] = (cons[i,2])
        fluxVec[i,2] = (cons[i,2]*prim[i,1] + prim[i,2])
        fluxVec[i,3] = (prim[i,1]*(cons[i,3] + prim[i,2]))
    end


    return fluxVec

end

function eigenFluxVectors(rhoRoe, uRoe, HRoe, aRoe, deltaRho, deltaU, deltaP, deltaA, fAVecs)
    #=This function takes (mostly) roe averaged quantities as inputs, and returns the eigenvalue corresponding vectors
        F_1 -> U~
        F_2 -> U~ + a~
        F_3 -> U~ - a~
    =#
    nFaces = size(fAVecs,1)

    mach = uRoe ./ aRoe

    if maximum(mach) > 1
        println("Passed sonic point!!!!")
    end

    #First step is to find the eigenvalues at each face
    eig1, eig2, eig3 = findEigenvalues(uRoe, aRoe, fAVecs)

    #Check, and if required correct them
    #This is the step I'm least sure about in the whole process
    #newEig1, newEig2, newEig3 = checkEigenvalues(eig1, eig2, eig3, deltaU, deltaA, nFaces)
    newEig1, newEig2, newEig3 = checkEigenvaluesUpdate(eig1, eig2, eig3, aRoe, nFaces)


    #Construct vectors

    f1 = Array{Float64, 2}(undef, nFaces, 3)
    f2 = Array{Float64, 2}(undef, nFaces, 3)
    f3 = Array{Float64, 2}(undef, nFaces, 3)

    for i in 1:nFaces
        #TODO: All the * 1 in this loop are places where the *1 should be replaced by the unit vector of the face area
        #First vector is easiest
        f1[i,1] = abs(newEig1[i]) * (deltaRho[i] - (deltaP[i]/(aRoe[i]^2))  )
        f1[i,2] = abs(newEig1[i]) * ( uRoe[i]*( deltaRho[i] - (deltaP[i]/(aRoe[i]^2)) ) + rhoRoe[i]*( deltaU[i] - 1 *deltaU[i] ) ) #The last term only matters for 2D/3D
        f1[i,3] = abs(newEig1[i]) * ( 0.5*uRoe[i]*( deltaRho[i] - (deltaP[i]/(aRoe[i]^2)) ) + rhoRoe[i]*( uRoe[i]*deltaU[i] - 1 * newEig1[i]*deltaU[i] ) ) #same as above

        f2_pre = abs(newEig2[i]) * ( deltaP[i]/(2*aRoe[i]^2) + rhoRoe[i]* 1 * deltaU[i]/(2*aRoe[i]) )
        f2[i,1] = f2_pre
        f2[i,2] = f2_pre * (uRoe[i] + 1 * aRoe[i])
        f2[i,3] = f2_pre * (HRoe[i] + 1 * uRoe[i] * aRoe[i])

        f3_pre = abs(newEig3[i]) * ( deltaP[i]/(2*aRoe[i]^2) - rhoRoe[i]* 1 * deltaU[i]/(2*aRoe[i]) )
        f3[i,1] = f3_pre
        f3[i,2] = f3_pre * (uRoe[i] - 1 * aRoe[i])
        f3[i,3] = f3_pre * (HRoe[i] - 1 * uRoe[i] * aRoe[i])
    end
    return f1, f2, f3
end

function findEigenvalues(u, a, fA)
    #=Given input of u, a, and area vectors find the eigenvalues of the matrix
    TODO: Extend this to 2D, use fA
    =#
    nFaces = size(fA, 1)
    eig1 = Array{Float64,1}(undef, nFaces)
    eig2 = Array{Float64,1}(undef, nFaces)
    eig3 = Array{Float64,1}(undef, nFaces)

    eig1 = u
    eig2 = u + a
    eig3 = u - a

    return eig1, eig2, eig3
end

function checkEigenvalues(eig1, eig2, eig3, dU, dA, n, K=0.04)
    #=Really not sure if this is the way to do it, but
        taking max of (U_r-U_l or 0)

        And using that to replace the U value???
    =#

    new_eig1 = Array{Float64, 1}(undef, n)
    new_eig2 = Array{Float64, 1}(undef, n)
    new_eig3 = Array{Float64, 1}(undef, n)

    for i in 1:n #Entropy fix I hope??
        if dU[i] > 0
            #replace
            eps = K * dU[i]
            if eps > abs(eig1[i])
                new_eig1[i] = (eig1[i]^2 + eps^2)/(2*eps)

                println("Updated eig1 value!")
            else
                new_eig1[i] = eig1[i]
            end
        else
            new_eig1[i] = eig1[i]
        end

        if (dU[i]+dA[i])>0
            #replace
            eps = K * (dU[i] + dA[i])
            #eps = K * c
            if eps > abs(eig2[i])
                new_eig2[i] = (eig2[i]^2 + eps^2)/(2*eps)
                println("Updated eig2 value!")
            else
                new_eig2[i] = eig2[i]
            end
        else
            new_eig2[i] = eig2[i]
        end

        if (dU[i]-dA[i])>0
            #replace
            eps = K * (dU[i] - dA[i])
            #eps = K * c
            if eps > abs(eig3[i])
                new_eig3[i] = (eig3[i]^2 + eps^2)/(2*eps)
                println("Updated eig3 value!")
            else
                new_eig3[i] = eig3[i]
            end
        else
            new_eig3[i] = eig3[i]
        end

    end

    #println("Done with eigenvalues!")

    return new_eig1, new_eig2, new_eig3

end

function checkEigenvaluesUpdate(eig1, eig2, eig3, sound, n, K=0.1)
    #=Really not sure if this is the way to do it, but
        taking max of (U_r-U_l or 0)

        And using that to replace the U value???
    =#
    #c = maximum(sound)

    new_eig1 = Array{Float64, 1}(undef, n)
    new_eig2 = Array{Float64, 1}(undef, n)
    new_eig3 = Array{Float64, 1}(undef, n)

    for i in 1:n #Entropy fix I hope??
        if K*sound[i] > abs(eig1[i])
            #replace
            #eps = K * sound[i]
            #new_eig1[i] = (eig1[i]^2 + eps^2)/(2*eps)
            new_eig1[i] = eig1[i]
            #println("Updated eig1 value!")
        else
            new_eig1[i] = eig1[i]
        end

        if K*sound[i] > abs(eig2[i])
            eps = K * sound[i]
            new_eig2[i] = (eig2[i]^2 + eps^2)/(2*eps)
            println("Updated eig2 value!")

        else
            new_eig2[i] = eig2[i]
        end

        if (K*sound[i]) > abs(eig3[i])
            #replace
            eps = K * sound[i]
            new_eig3[i] = (eig3[i]^2 + eps^2)/(2*eps)
            println("Updated eig3 value!")

        else
            new_eig3[i] = eig3[i]
        end

    end

    println("Done with eigenvalues!")

    return new_eig1, new_eig2, new_eig3

end

function findFluxes(lF, rF, e1F, e2F, e3F, fA)
    #=This function combines the decomposed vectors into a single flux vector
    =#
    n = size(fA, 1)
    xMassFlux = Array{Float64, 1}(undef, n)
    xMomFlux = Array{Float64, 1}(undef, n)
    xEneFlux = Array{Float64, 1}(undef, n)

    for i in 1:n
        xMassFlux[i] = 0.5 * (lF[i,1] + rF[i,1] - e1F[i,1] - e2F[i,1] - e3F[i,1])
        xMomFlux[i] = 0.5 * (lF[i,2] + rF[i,2] - e1F[i,2] - e2F[i,2] - e3F[i,2])
        xEneFlux[i] = 0.5 * (lF[i,3] + rF[i,3] - e1F[i,3] - e2F[i,3] - e3F[i,3])
    end

    return xMassFlux, xMomFlux, xEneFlux
end


######################### Convective Term Things #######################

function convectFlux(mesh, rho, xMom, eV2, dx; gamma=1.4, verbose=false, debug=false)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells = size(cells, 1)
    nFaces = size(faces, 1)

    #Use MUSCL differencing to get values at left/right of each face. Boundaries are treated inside this function
    consLeft, primsLeft, consRight, primsRight = musclDifference(mesh, rho, xMom, eV2, dx, debug=debug)

    rhoLeft = consLeft[:,1]
    rhoRight = consRight[:,1]

    uLeft = primsLeft[:,1]
    uRight = primsRight[:,1]

    pLeft = primsLeft[:,2]
    pRight = primsRight[:,2]

    HLeft = primsLeft[:,3]
    HRight = primsRight[:,3]

    aLeft = (gamma-1)*(HLeft - 0.5*uLeft.^2)
    aRight = (gamma-1)*(HRight - 0.5*uRight.^2)

    #println("Checkpoint 2!\n")


    #Find ROE averaged quantities at each face
    rhoRoe, uRoe, HRoe, aRoe = roeAveraged(nFaces, rhoLeft, rhoRight, uLeft, uRight, HLeft, HRight)

    #Construct flux partial vectors
    leftFlux = fluxVector(consLeft, primsLeft, fAVecs)
    rightFlux = fluxVector(consRight, primsRight, fAVecs)

    deltaRho = rhoRight - rhoLeft
    deltaU = uRight - uLeft
    deltaP = pRight - pLeft
    deltaA = aRight - aLeft

    #println("Checkpoint 3!\n")

    eigen1Flux, eigen2Flux, eigen3Flux = eigenFluxVectors(rhoRoe, uRoe, HRoe, aRoe, deltaRho, deltaU, deltaP, deltaA, fAVecs) #Eigenvalue check and correction for entropy happens in this step

    #Combines decomposed flux vectors into flux at every step
    xMassFlux, xMomFlux, xeV2Flux = findFluxes(leftFlux, rightFlux, eigen1Flux, eigen2Flux, eigen3Flux, fAVecs )

    #println("Checkpoint 4!\n")



    fluxVars = [ xMassFlux, xMomFlux, xeV2Flux ]


    if debug
        println("xMass Flux: $xMassFlux")
        println("xMom Flux: $xMomFlux")
        println("xEv2 Flux: $xeV2Flux")
    end

    return fluxVars

end



######################### Solvers #######################
function upwindFVMRoe1D(mesh, P, T, U; initDt=0.001, endTime=0.14267, targetCFL=0.2, gamma=1.4, R=287.05, Cp=1005, debug=false, verbose=false, timeStep="EulerExp")
    ######### MESH ############
    # Extract mesh into local variables for readability
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells = size(cells, 1)
    nFaces = size(faces, 1)
    nBdryFaces = size(boundaryFaces, 1)
    bdryFaceIndices = Array(nFaces-nBdryFaces:nFaces)




    #dx = cCenters[2] - cCenters[1]
    dx = 1 / nCells

    ########### Variable Arrays #############
    # State variables, values are the averages/cell center values
    rho = Array{Float64, 1}(undef, nCells)
    xMom = Array{Float64, 1}(undef, nCells)
    eV2 = Array{Float64, 1}(undef, nCells)
    stateVars = [ rho, xMom, eV2 ]
    # Flux variables for equations 2 and 3 (xMom is flux variable for the continuity eqn)
    rhoU2p = Array{Float64, 1}(undef, nCells)
    rhoUeV2PU = Array{Float64, 1}(undef, nCells)


    stateVars_old = Array{Float64,2}(undef, 3, nCells)



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

    if verbose
        #println("-------------------------------------------------------------------------------------------------------")
        println("------------------------------  Intialize complete, commencing time stepping  -------------------------")
        println("-------------------------------------------------------------------------------------------------------")
        println("")
        println("CurrentTime ------ dt ------ max mom res")
    end

    for i in 1:3
        stateVars_old[i,:] = stateVars[i][:]
    end


    ########### SOLVER ###########
    dt = initDt
    currTime = 0

    if timeStep=="EulerExp"
        while currTime < endTime

            if (endTime - currTime) < dt
                dt = endTime - currTime
            end

            # Calculate fluxes through each face
            # TODO: y and z momemtum-fluxes + equations
            # xMassFlux, xMomFlux, xeV2Flux = upwindInterp(mesh, U, xMom, rhoU2p, rhoUeV2PU)
            #xMassFlux, xMomFlux, xeV2Flux, faceP = linInterp(mesh, xMom, rhoU2p, rhoUeV2PU, P)

            fluxVars = convectFlux(mesh, rho, xMom, eV2, dx; gamma=gamma, verbose=verbose, debug=debug)

            #[xMassFlux, xMomFlux, xEnerFlux] = fluxVars

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


            stateVars[2][1] = stateVars[2][2]
            stateVars[2][nCells] = stateVars[2][nCells-1]

            #Unpack state variables
            #Not sure if this step is required??
            rho = stateVars[1][:]
            xMom = stateVars[2][:]
            eV2 = stateVars[3][:]

            # Decode primitive values
            for i in 1:nCells
                P[i], T[i], U[i], rhoU2p[i], rhoUeV2PU[i] = decodePrimitives3D(rho[i], xMom[i], eV2[i])
            end

            if verbose
                println("\n")
                println(" ", currTime, " --- ", dt, " --- ", maximum(fluxVars[1][:]))
                println("\n")
            end

            currTime += dt

            #println("Checkpoint 6!\n")

            #Catch negative temps for now - until can find cause
            #=
            for z in 1:nCells
                if T[z]< 0
                    T[z] = 0.0001222
                end
            end
            =#

            ############## CFL Calculation, timestep adjustment #############
            maxCFL = 0
            #println("Temp is: ", T)
            for i in 1:nCells
                dx = 1/nCells
                maxCFL = max(maxCFL, CFL(U[i], T[i], dt, dx, gamma, R))
            end
            # Adjust time step to approach target CFL
            dt *= ((targetCFL/maxCFL - 1)/5+1)

        end
    elseif timeStep=="RK4"
        while currTime < endTime

            if (endTime - currTime) < dt
                dt = endTime - currTime
            end

            fluxResid_1 = Array{Float64, 2}(undef, 3, nCells)
            fluxResid_2 = Array{Float64, 2}(undef, 3, nCells)
            fluxResid_3 = Array{Float64, 2}(undef, 3, nCells)
            fluxResid_4 = Array{Float64, 2}(undef, 3, nCells)

            stateVars_1 = Array{Float64, 2}(undef, 3, nCells)
            stateVars_2 = Array{Float64, 2}(undef, 3, nCells)
            stateVars_3 = Array{Float64, 2}(undef, 3, nCells)
            stateVars_new = Array{Float64, 2}(undef, 3, nCells)

            ##########################  Step 1 #################################
            # Calculate fluxes through each face
            fluxVars_1 = convectFlux(mesh, rho, xMom, eV2, dx; gamma=gamma, verbose=verbose, debug=debug)

            # Use fluxes to update values in each cell
            for i in 1:(nFaces-nBdryFaces)
                fA = mag(fAVecs[i]) #Face Area

                ownerCell = faces[i][1]
                neighbourCell = faces[i][2]
                for v in 1:3
                    #Not sure if I need the +/-?? I think that's accounted for in the dot area vector? Or maybe not?
                    #Also don't think you need to divide each face by cell vol, only the final result?
                    fluxResid_1[v,ownerCell] -= fluxVars_1[v][i]*fA/cVols[ownerCell]
                    fluxResid_1[v,neighbourCell] += fluxVars_1[v][i]*fA/cVols[neighbourCell]
                end
            end

            stateVars_1 = stateVars_old + fluxResid_1*dt*0.5

            stateVars_1[2,1] = stateVars_1[2,2]
            stateVars_1[2,nCells] = stateVars_1[2,nCells-1]

            ##########################  Step 2 #################################


            #Unpack state variables
            #Not sure if this step is required??
            rho_1 = stateVars_1[1,:]
            xMom_1 = stateVars_1[2,:]
            eV2_1 = stateVars_1[3,:]

            fluxVars_2 = convectFlux(mesh, rho_1, xMom_1, eV2_1, dx; gamma=gamma, verbose=verbose, debug=debug)

            # Use fluxes to update values in each cell
            for i in 1:(nFaces-nBdryFaces)
                fA = mag(fAVecs[i]) #Face Area

                ownerCell = faces[i][1]
                neighbourCell = faces[i][2]
                for v in 1:3
                    #Not sure if I need the +/-?? I think that's accounted for in the dot area vector? Or maybe not?
                    #Also don't think you need to divide each face by cell vol, only the final result?
                    fluxResid_2[v,ownerCell] -= fluxVars_2[v][i]*fA/cVols[ownerCell]
                    fluxResid_2[v,neighbourCell] += fluxVars_2[v][i]*fA/cVols[neighbourCell]
                end
            end

            stateVars_new = stateVars_old + fluxResid_2*dt

            stateVars_new[2,1] = stateVars_new[2,2]
            stateVars_new[2,nCells] = stateVars_new[2,nCells-1]
            #=
            ##########################  Step 3 #################################

            #Unpack state variables
            #Not sure if this step is required??
            rho_2 = stateVars_2[1,:]
            xMom_2 = stateVars_2[2,:]
            eV2_2 = stateVars_2[3,:]

            fluxVars_3 = convectFlux(mesh, rho_2, xMom_2, eV2_2, dx; gamma=gamma, verbose=verbose, debug=debug)

            # Use fluxes to update values in each cell
            for i in 1:(nFaces-nBdryFaces)
                fA = mag(fAVecs[i]) #Face Area

                ownerCell = faces[i][1]
                neighbourCell = faces[i][2]
                for v in 1:3
                    #Not sure if I need the +/-?? I think that's accounted for in the dot area vector? Or maybe not?
                    #Also don't think you need to divide each face by cell vol, only the final result?
                    fluxResid_3[v,ownerCell] -= fluxVars_3[v][i]*fA/cVols[ownerCell]
                    fluxResid_3[v,neighbourCell] += fluxVars_3[v][i]*fA/cVols[neighbourCell]
                end
            end

            stateVars_3 = stateVars_old + fluxResid_3*dt

            stateVars_3[2,1] = stateVars_3[2,2]
            stateVars_3[2,nCells] = stateVars_3[2,nCells-1]

            ##########################  Step 4 #################################

            #Unpack state variables
            #Not sure if this step is required??
            rho_3 = stateVars_3[1,:]
            xMom_3 = stateVars_3[2,:]
            eV2_3 = stateVars_3[3,:]

            fluxVars_4 = convectFlux(mesh, rho_3, xMom_3, eV2_3, dx; gamma=gamma, verbose=verbose, debug=debug)

            # Use fluxes to update values in each cell
            for i in 1:(nFaces-nBdryFaces)
                fA = mag(fAVecs[i]) #Face Area

                ownerCell = faces[i][1]
                neighbourCell = faces[i][2]
                for v in 1:3
                    #Not sure if I need the +/-?? I think that's accounted for in the dot area vector? Or maybe not?
                    #Also don't think you need to divide each face by cell vol, only the final result?
                    fluxResid_4[v,ownerCell] -= fluxVars_4[v][i]*fA/cVols[ownerCell]
                    fluxResid_4[v,neighbourCell] += fluxVars_4[v][i]*fA/cVols[neighbourCell]
                end
            end

            stateVars_new = stateVars_old + (dt/6)*(fluxResid_1 + 2*fluxResid_2 + 2*fluxResid_3 + fluxResid_4)

            stateVars_new[2,1] = stateVars_new[2,2]
            stateVars_new[2,nCells] = stateVars_new[2,nCells-1]
            =#

            ##########################  Decode for next step ############################
            rho = deepcopy(stateVars_new[1,:])
            xMom = deepcopy(stateVars_new[2,:])
            eV2 = deepcopy(stateVars_new[3,:])

            # Decode primitive values
            for i in 1:nCells
                P[i], T[i], U[i], rhoU2p[i], rhoUeV2PU[i] = decodePrimitives3D(rho[i], xMom[i], eV2[i])
            end

            if verbose
                println("\n")
                println(" ", currTime, " --- ", dt, " --- ")
                println("\n")
            end

            currTime += dt

            ############## CFL Calculation, timestep adjustment #############
            maxCFL = 0
            #println("Temp is: ", T)
            for i in 1:nCells
                dx = 1/nCells
                maxCFL = max(maxCFL, CFL(U[i], T[i], dt, dx, gamma, R))
            end
            # Adjust time step to approach target CFL
            dt *= ((targetCFL/maxCFL - 1)/5+1)

        end

    end

    #println("P is: ", P)
    #println("U is: ", U)
    #println("Temp is: ", T)



    return P, U, T, rho
end