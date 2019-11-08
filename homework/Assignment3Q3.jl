include("LinearSolvers.jl")

# zeta will be referred to as z or Z

# Initial value functions for T and eta(Z)
function initZ(eta)
    return eta/20
end

function initT(eta)
    return 1 - eta/20
end

function addCentralDerivative(stencil, multiple, position, matrix)
    stencil2 = stencil * multiple
    stencilLength = size(stencil, 2)
    stencilBandwidth = floor(Int32, (stencilLength - 1) / 2)
    for i in 1:stencilLength
        matrix[position, position+i-stencilBandwidth-1] += stencil2[i]
    end

    return matrix
end

function addStencil(stencil, multiple, row, col, matrix)
    stencil2 = stencil * multiple
    for i in 1:size(stencil,2)
        matrix[row, col + i - 1] += stencil2[i]
    end
    return matrix
end

function solve_Q3(nNodes, z1, t1, t2, epsilon=0.0001, iterLimit=200)
    eta1 = 0
    eta2 = 20
    dx = (eta2 - eta1) / (nNodes + 1)

    #Create initial value vectors
    # 3 and 2 additional equations for the boundary conditions
    additionalZEqns = 3
    additionalTEqns = 2
    z = Array{Float64, 1}(undef, nNodes + additionalZEqns)
    t = Array{Float64, 1}(undef, nNodes + additionalTEqns)
    zNodes = size(z,1)
    tNodes = size(t,1)

    for i in 1:(nNodes+additionalZEqns)
        eta = dx*(i-2)
        
        z[i] = initZ(eta)        
        if i > 1 && i <= nNodes + 3
            t[i-1] = initT(eta)
        end
    end

    # Second order central stencils for each derivative
    cdn3 = (1/(2* dx*dx*dx)) .* [ -1 2 0 -2 1 ]
    cdn2 = (1/(dx*dx)) .* [ 1 -2 1 ]
    cdn = (1/(2*dx)) .* [ -1 0 1 ]

    maxDx = 1
    iterationCounter = 1
    while maxDx > epsilon
        ####################### Build z Matrix ########################
        zMatrix = zeros(Float64, zNodes+tNodes, zNodes+tNodes+1)

        # Apply the appropriate stencils for interior nodes
        for i in 3:(zNodes-2)
            zMatrix = addCentralDerivative(cdn3, 1, i, zMatrix)
            secondDerMultiple = 3*z[i]
            zMatrix = addCentralDerivative(cdn2, secondDerMultiple, i, zMatrix)
            thirdDerMultiple = -2*( (z[i+1] - z[i-1]) / (2*dx) )
            zMatrix = addCentralDerivative(cdn, thirdDerMultiple, i, zMatrix)
            zMatrix[i, zNodes+i-1] = 1
        end

        # Zeta Boundary Conditions: First row: derivative condition
        cdn = [ -1 0 1 ]
        zMatrix = addStencil(cdn, 1, 1, 1, zMatrix)

        # Second row is equation for first boundary node
        zMatrix[2, 2] = 1
        zMatrix[2, zNodes+tNodes+1] = z1

        # Second last row: Modified third derivative term in equation
         cdn3Offset = [ -1 7 -22 34 -25 7 ] / (2 * dx*dx*dx) # f-4 to f+1
         zMatrix = addStencil(cdn3Offset, 1, zNodes-1, zNodes-5, zMatrix)
         
         # Add rest of equation as normal
         i = zNodes-1
         secondDerMultiple = 3*z[i]
         zMatrix = addCentralDerivative(cdn2, secondDerMultiple, i, zMatrix)
         thirdDerMultiple = -2*( (z[i+1] - z[i-1]) / (2*dx) )
         zMatrix = addCentralDerivative(cdn, thirdDerMultiple, i, zMatrix)
         zMatrix[i, zNodes+i-1] = 1

        # Last row: derivative condition
        cdnR = [ 1 -4 3 ]
        zMatrix = addStencil(cdnR, 1, zNodes, zNodes-2, zMatrix)


        ######################## Build t Matrix ##########################
        # Apply appropriate stencils for interior nodes
        for i in 2+zNodes:(zNodes+tNodes-1)
            zMatrix = addCentralDerivative(cdn2, 1, i, zMatrix)
            firstDerMultiple = 3 * 0.71 * z[i+1 - zNodes]
            zMatrix = addCentralDerivative(cdn, firstDerMultiple, i, zMatrix)
        end

        #Boundary Conditions: For the first row, value of T0 is known to be t1
        zMatrix[1+zNodes,1+zNodes] = 1
        zMatrix[1+zNodes, zNodes+tNodes+1] = t1

        # For the last row, value of T1 is known to be t2
        zMatrix[tNodes+zNodes, tNodes+zNodes] = 1
        zMatrix[tNodes+zNodes, zNodes+tNodes+1] = t2

        ######################## Solve ##########################
        # With sufficient equation reordering we could make this use a thomas-type solver, but use Gauss Elim for now for simplicity
        newTZ = Solve_GaussElim!(zMatrix)

        newZ = newTZ[1:zNodes]
        newT = newTZ[zNodes+1:zNodes+tNodes]

        ################### Evaluate dt, dz #####################
        maxDx = 0
        for i in 3:nNodes+1
            dz = abs(newZ[i] - z[i])
            dt = abs(newT[i-1] - t[i-1])
            maxDx = max(max(dz, dt), maxDx)
        end
       
        println("Iteration $iterationCounter, Max dx = $maxDx")
        println("")
        
        # Under-relaxation to stabilize the solution
        t = (newT + t) / 2
        z = (newZ + z) / 2

        iterationCounter += 1

        if iterationCounter == iterLimit
            print("ERROR: solution did not converge in $iterLimit iterations")
            return z[2:nNodes+3], t
        end
    end

    # Cuts out-of-domain values off of z before returning
    return z[2:nNodes+3], t
end