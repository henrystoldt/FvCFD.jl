include("LinearSolvers.jl")
using Plots
using Plots.PlotMeasures
plotly()

# zeta will be referred to as z or Z

function residual1(d3Zdn, d2Zdn, dZdn, T, Z)
    return d3Zdn + 3*Z*d2Zdn - 2 *(dZdn*dZdn) + T
end

function residual2(d2Tdn, dTdn, Z)
    return d2Tdn + 3*0.71*Z*dTdn
end

# Initial value functions for T and eta(Z)
function initZ(eta)
    return eta/20
end

function initT(eta)
    return 1 - eta/20
end

#nNodes = number of interior nodes
function solve_Q3_SecondOrder(nNodes, epsilon=0.0001)    
    eta1 = 0
    eta2 = 20
    dx = (eta2 - eta1) / (nNodes + 1)
    
    println("Nodes = $nNodes, dx=$dx")

    nVars = nNodes * 2

    #Create initial value vectors
    x = Array{Float64, 1}(undef, nVars)
    for i in 1:nNodes
        eta = dx*i
        x[2*i - 1] = initZ(eta)
        x[2*i] = initT(eta)
    end
    println("Initial X-vector: $x")

    # Second order central stencils for each derivative
    cdn3 = (1/(2* dx*dx*dx)) .* [ -1 2 0 -2 1 ]
    cdn2 = (1/(dx*dx)) .* [ 1 -2 1 ]
    cdn = (1/(2*dx)) .* [ -1 0 1 ]

    println("3rd Derivative Stencil: $cdn3")
    println("2nd Derivative Stencil: $cdn2")
    println("1st Derivative Stencil: $cdn")
    
    # Build the linear parts of the matrix (first term of both equations)
    # Make a single stencil, replicate it for each row, adjust the first and last for the boundary conditions
    bandwidth = size(cdn3, 2) * 2
    
    # Order of variables is z1, T1, z2, T2, z3, T3 etc...
    # Create template z eqn
    thomasZRow = zeros(bandwidth)
    # Add third derivative of zeta
    for i in 1:size(cdn3, 2)
        thomasZRow[2*i-1] += cdn3[i]
    end
    # Add Ti
    Zcenter = floor(Int32, bandwidth / 2)
    Tcenter = Zcenter + 1
    thomasZRow[Tcenter] += 1

    println("Bandwidth: $bandwidth")
    println("Center: $Zcenter")

    println("Template Z Eqn Stencil:")
    println(thomasZRow)

    # Create template T eqn
    thomasTRow = zeros(bandwidth)
    #Add second derivative of T
    offset = floor(Int32, (size(cdn3, 2) - size(cdn2, 2))/2)
    for i in 1:size(cdn2,2)
        thomasTRow[2*(i+offset)] += cdn2[i]
    end

    println("Template T Eqn Stencil:")
    println(thomasTRow)

    # Is stored with z and T values interleaved. Z at odd indices, T at even
    function getXVal(i)
        # Boundary conditions for T
        # Given
        if i == nVars + 2
            return 1
        elseif i == 0
            return 0
        # Boundary conditions for z
        # From second order central difference and derivative boundary condition
        elseif i == nVars + 3
            return x[nVars-1]
        elseif i == -3
            return x[1]
        # From second order backward/forward differences and derivative boundary condition
        elseif i == nVars + 1
            return 0 #(4*x[nVars-1] - x[nVars-3]) / 3
        elseif i == -1
            return (4*x[1] - x[3]) / 3
        elseif i == -2 || i == nVars + 2 || i == 0 || i == nVars + 4
            return 0
        else
            return x[i]
        end
    end
    
    println("")

    maxDx = 1
    iterationCounter = 1
    while maxDx > epsilon
        # Create matrix, alternating T and zeta equations
        matrix = zeros(Float64, nVars, nVars+1)
        semiSpan = floor(Int32, (bandwidth - 2) / 2)

        # Calculate values of nonlinear terms for this iteration, lagging the lowest order terms in each nonlinear term
        for r in 1:nNodes
            # Initialize result column to zero
            matrix[2*r - 1, nVars+1] = 0
            matrix[2*r, nVars+1] = 0

            # Copy numbers in to matrix
            for colIndex in 1:bandwidth
                column = 2*r - 1 + colIndex - Zcenter
                # println("Col $column")

                #Evaluate and move elements to RHS
                if column < 1 || column > nVars
                    matrix[2*r - 1, nVars+1] += -1 * thomasZRow[colIndex] * getXVal(column)
                    matrix[2*r, nVars+1] += -1 * thomasTRow[colIndex] * getXVal(column)
                    # line = matrix[3,:]
                    # println("Added to RHS: $line")
                    continue
                end

                # Populate LHS
                matrix[2*r - 1, column] = thomasZRow[colIndex]
                matrix[2*r, column] = thomasTRow[colIndex]
                # line = matrix[3,:]
                # println("Added to LHS: $line")
            end


            # line = matrix[3,:]
            # println(line)

            # Add all the 3-element derivative discretizations
            # Adding second derivative of zeta
            cdn2_Z2 = cdn2 .* 3 * getXVal(2*r - 1)
            # Adding squared first derivative of zeta
            cdn1_Z1 = cdn .* -2 * (getXVal(2*r + 1) - getXVal(2*r - 3)) / (2*dx)
            #Adding first derivative of T1
            cdn_T1 = cdn .* 3 * 0.71 * getXVal(2*r - 1)
            semiSpan = floor(Int32, (size(cdn2, 2) -1) /2)
            for a in 1:size(cdn2, 2)
                offCenter = a - (semiSpan + 1)
                col = 2*r + 2*offCenter - 1
                # println("Col $col")
                if col < 1 || col > nVars
                    matrix[2*r - 1, nVars+1] += -1 * cdn2_Z2[a] * getXVal(col)
                    matrix[2*r - 1, nVars+1] += -1 * cdn1_Z1[a] * getXVal(col)
                    matrix[2*r, nVars+1] += -1 * cdn_T1[a] * getXVal(col+1)
                    # line = matrix[3,:]
                    # println("Added to RHS: $line")
                    continue
                end
                matrix[2*r-1, col] += cdn2_Z2[a]
                matrix[2*r-1, col] += cdn1_Z1[a]
                matrix[2*r, col] += cdn_T1[a]
                # line = matrix[3,:]
                # println("Added to LHS: $line")
            end
        end

        # printMatrix(matrix)

        # Solve matrix
        newX = Solve_GaussElim!(matrix)

        # Calculate residuals/dxs
        maxDx = 0
        totalDx = 0
        for i in 1:size(newX, 1)
            dx = abs(newX[i] - x[i])
            if dx > maxDx
                maxDx = dx
            end
            totalDx += dx
        end
        meanDx = totalDx/nVars

        x = newX
        avgX = sum(x)/nVars

        println("Iteration $iterationCounter, Max dx = $maxDx, Mean dx = $meanDx, Avg x = $avgX")
        println("")
        iterationCounter += 1
    end

    return x
end

function addCentralDerivative(stencil, multiple, position, matrix)
    stencil2 = stencil * multiple
    # Assuming stencil is symmetricl
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

function solve_Q3_SecondOrder2(nNodes, z1, t1, t2, epsilon=0.0001, debug=false)
    eta1 = 0
    eta2 = 20
    dx = (eta2 - eta1) / (nNodes + 1)
    
    if debug
        println("Nodes = $nNodes, dx = $dx")
    end

    #Create initial value vectors
    additionalZEqns = 2
    z = Array{Float64, 1}(undef, nNodes+additionalZEqns)
    t = Array{Float64, 1}(undef, nNodes)
    zNodes = size(z,1)

    for i in 1:(nNodes+additionalZEqns)
        eta = dx*(i-1)
        
        z[i] = initZ(eta)
        
        if i > 1 && i < nNodes + additionalZEqns
            t[i-1] = initT(eta)
        end
    end

    if debug
        println("")
        println("Initial z-vector: $z")
        println("Initial t-vector: $t")
    end

    # Second order central stencils for each derivative
    cdn3 = (1/(2* dx*dx*dx)) .* [ -1 2 0 -2 1 ]
    cdn2 = (1/(dx*dx)) .* [ 1 -2 1 ]
    cdn = (1/(2*dx)) .* [ -1 0 1 ]

    if debug
        println("")
        println("3rd Derivative Stencil: $cdn3")
        println("2nd Derivative Stencil: $cdn2")
        println("1st Derivative Stencil: $cdn")
        println("")
    end

    # Solve in a block-iterative method, where we solve a matrix representing equation 1, then another matrix representing equation 2, then iterate again
    maxDx = 1
    iterationCounter = 1
    while maxDx > epsilon
        ################### Build and solve z Matrix ###################
        # One additional equation required on the RHB, where we only have a derivative boundary condition
        zMatrix = zeros(Float64, zNodes, zNodes+1)

        if debug
            println("Init:")
            printMatrix(zMatrix)
        end

        # Will simply apply the stencils for interior nodes
        for i in 3:(zNodes-2)
            zMatrix = addCentralDerivative(cdn3, 1, i, zMatrix)
            secondDerMultiple = 3*z[i]
            zMatrix = addCentralDerivative(cdn2, secondDerMultiple, i, zMatrix)
            thirdDerMultiple = -2*( (z[i+1] - z[i-1]) / (2*dx) )
            zMatrix = addCentralDerivative(cdn, thirdDerMultiple, i, zMatrix)
            zMatrix[i, zNodes+1] = -t[i-1]
        end

        if debug
            println("Interior Stencils Applied:")
            printMatrix(zMatrix)
        end

        # Now need to apply boundary conditions, make the first two and the last two rows
        # First row is equation for first boundary node
        row = 1
        col = 1
        cdn3LL = [ 0 -4 1 ] / (2 * dx*dx*dx)
        zMatrix = addStencil(cdn3LL, 1, row, col, zMatrix)
        RHS = -3 * z1 / (2* dx*dx*dx)
        zMatrix[1, zNodes+1] += RHS

        if debug
            println("Top Row 1:")
            printMatrix(zMatrix)
        end
        
        cdn2L = [ 0 2 ] / (dx*dx)
        secondDerMultiple = 3*z1
        zMatrix = addStencil(cdn2L, secondDerMultiple, row, col, zMatrix)
        RHS = secondDerMultiple * 2 * z1 / (dx*dx)
        zMatrix[1, zNodes+1] += RHS

        if debug
            println("Top Row 2:")
            printMatrix(zMatrix)
        end

        zMatrix[1, zNodes+1] += -t1

        if debug
            println("Top Row Final Boundary Conditions Applied:")
            printMatrix(zMatrix)
        end

        # Now need to do the second row, only modifying the first term
        row = 2
        col = 1
        cdn3R = [ 2 -1 -2 1 ] / (2 * dx*dx*dx )
        zMatrix = addStencil(cdn3R, 1, row, col, zMatrix)

        if debug
            println("Second Row 1:")
            printMatrix(zMatrix)
        end

        i = row
        secondDerMultiple = 3*z[i]
        zMatrix = addCentralDerivative(cdn2, secondDerMultiple, i, zMatrix)
        thirdDerMultiple = -2*( (z[i+1] - z[i-1]) / (2*dx) )
        zMatrix = addCentralDerivative(cdn, thirdDerMultiple, i, zMatrix)
        zMatrix[i, zNodes+1] += -t[i]
        
        if debug   
            println("Second Row Final:")
            printMatrix(zMatrix)
        end

        # On the right hand boundary, from the boundary conditions, z[n+1] = z[n-1] and z[n-2] = z[n+2]
        row = zNodes
        col = zNodes - 2
        cdn3RR = [ -1 4 -3 ] / (2 * dx*dx*dx)
        zMatrix = addStencil(cdn3RR, 1, row, col, zMatrix)

        if debug
            println("Last Row 1:")
            printMatrix(zMatrix)
        end

        col = zNodes - 1
        cdn2RR = [ 2 -2 ] / ( dx*dx )
        cdn2RRMul = 3*z[row]
        zMatrix = addStencil(cdn2RR, cdn2RRMul, row, col, zMatrix)

        if debug
            println("Last Row 2:")
            printMatrix(zMatrix)
        end
        
        RHS = -t2
        zMatrix[row, zNodes+1] += RHS

        if debug
            println("Last Row Final:")
            printMatrix(zMatrix)
        end

        # One node inside the right hand boundary, only need to adjust the third derivative term
        row = zNodes - 1
        col = zNodes - 3
        cdn3R = [ -1 2 1 -2 ] / (2 * dx*dx*dx )
        zMatrix = addStencil(cdn3R, 1, row, col, zMatrix)

        if debug
            println("Second Last Row 1:")
            printMatrix(zMatrix)
        end

        i = row
        secondDerMultiple = 3*z[i]
        zMatrix = addCentralDerivative(cdn2, secondDerMultiple, i, zMatrix)
        thirdDerMultiple = -2*( (z[i+1] - z[i-1]) / (2*dx) )
        zMatrix = addCentralDerivative(cdn, thirdDerMultiple, i, zMatrix)
        zMatrix[i, zNodes+1] += -t[i-1]

        if debug
            println("Second Last Row Final:")
            printMatrix(zMatrix)
        end

        # Now zMatrix should be complete, solve it
        newZ = Solve_GaussElim!(zMatrix)

        if debug
            println("New Zeta: $newZ")
        end

        #################### Build and solve t Matrix #####################
        # No additional equations required for t, dirichlet on both sides
        tMatrix = zeros(Float64, nNodes, nNodes+1)

        # Apply appropriate stencils for interior nodes
        for i in 2:(nNodes-1)
            tMatrix = addCentralDerivative(cdn2, 1, i, tMatrix)
            firstDerMultiple = 3 * 0.71 * newZ[i]
            tMatrix = addCentralDerivative(cdn, firstDerMultiple, i, tMatrix)
        end

        if debug
            println("Interior Stencils Applied:")
            printMatrix(tMatrix)
        end

        # Now apply boundary conditions
        # For the first row, value of T0 is known to be t1
        row = 1
        col = 1

        cdn2R = [ -2 1 ] / (dx*dx)
        tMatrix = addStencil(cdn2R, 1, row, col, tMatrix)
        RHS = -t1 / (dx*dx)
        tMatrix[1, nNodes+1] += -RHS

        cdnR = [ 0 1 ] / (2*dx)
        cdnRMul = 3*0.71*newZ[row]
        tMatrix = addStencil(cdnR, cdnRMul, row, col, tMatrix)
        RHS = cdnRMul * t1 / (2*dx)
        tMatrix[1, nNodes+1] += RHS
        
        if debug
            println("Second Last Row Final:")
            printMatrix(tMatrix)
        end

        #For the last row, value of Tn is known to be t2
        row = nNodes
        col = nNodes-1
        
        cdn2R = [-2 1] / (dx*dx)
        tMatrix = addStencil(cdn2R, 1, row, col, tMatrix)
        RHS = -t2 / (dx*dx)
        tMatrix[nNodes, nNodes+1] += RHS

        if debug
            println("Last Row 1:")
            printMatrix(tMatrix)
        end

        cdnR = [ 0 -1 ] / (2*dx)
        cdnRMul = 3 * 0.71 * newZ[nNodes]
        tMatrix = addStencil(cdnR, cdnRMul, row, col, tMatrix)
        RHS = -t2 * cdnRMul / (2*dx)
        tMatrix[nNodes, nNodes+1] += RHS

        if debug
            println("Last Row Final:")
            printMatrix(tMatrix)
        end

        # Solve
        newT = Solve_GaussElim!(tMatrix)

        if debug
            println("New T: $newT")
        end

        ############## Evaluate dt, dz ##############
        maxDx = 0
        for i in 1:nNodes
            dz = abs(newZ[i] - z[i])
            dt = abs(newT[i] - t[i])
            maxDx = max(max(dz, dt), maxDx)
        end
       
        println("Iteration $iterationCounter, Max dx = $maxDx")
        
        iterationCounter += 1
        t = newT
        z = newZ

        println("")

        if iterationCounter == 200
            break
        end
    end

    return z, t
end

function solve_Q3_SecondOrder3(nNodes, z1, t1, t2, epsilon=0.0001, debug=false)
    eta1 = 0
    eta2 = 20
    dx = (eta2 - eta1) / (nNodes + 1)
    
    if debug
        println("Nodes = $nNodes, dx = $dx")
    end

    #Create initial value vectors
    additionalZEqns = 3
    additionalTEqns = 2
    z = Array{Float64, 1}(undef, nNodes + additionalZEqns)
    t = Array{Float64, 1}(undef, nNodes + additionalTEqns)
    zNodes = size(z,1)
    tNodes = size(t,1)

    for i in 1:(nNodes+additionalZEqns)
        eta = dx*(i-2)
        
        z[i] = initZ(eta)
        
        if i > 1 && i <= nNodes + 1
            t[i-1] = initT(eta)
        end
    end

    if debug
        println("")
        println("Initial z-vector: $z")
        println("Initial t-vector: $t")
    end

    # Second order central stencils for each derivative
    cdn3 = (1/(2* dx*dx*dx)) .* [ -1 2 0 -2 1 ]
    cdn2 = (1/(dx*dx)) .* [ 1 -2 1 ]
    cdn = (1/(2*dx)) .* [ -1 0 1 ]

    if debug
        println("")
        println("3rd Derivative Stencil: $cdn3")
        println("2nd Derivative Stencil: $cdn2")
        println("1st Derivative Stencil: $cdn")
        println("")
    end

    # Solve in a block-iterative method, where we solve a matrix representing equation 1, then another matrix representing equation 2, then iterate again
    maxDx = 1
    iterationCounter = 1
    while maxDx > epsilon
        ################### Build and solve z Matrix ###################
        # One additional equation required on the RHB, where we only have a derivative boundary condition
        zMatrix = zeros(Float64, zNodes, zNodes+1)

        if debug
            println("Init:")
            printMatrix(zMatrix)
        end

        # Will simply apply the stencils for interior nodes
        for i in 3:(zNodes-2)
            zMatrix = addCentralDerivative(cdn3, 1, i, zMatrix)
            secondDerMultiple = 3*z[i]
            zMatrix = addCentralDerivative(cdn2, secondDerMultiple, i, zMatrix)
            thirdDerMultiple = -2*( (z[i+1] - z[i-1]) / (2*dx) )
            zMatrix = addCentralDerivative(cdn, thirdDerMultiple, i, zMatrix)
            zMatrix[i, zNodes+1] = -t[i-2]
        end

        if debug
            println("Interior Stencils Applied:")
            printMatrix(zMatrix)
        end

        # Now need to apply boundary conditions, make the first two and the last two rows
        # Second row: derivative condition
        cdn = [ -1 0 1 ]
        zMatrix = addStencil(cdn, 1, 1, 1, zMatrix)

        # Second row is equation for first boundary node
        zMatrix[2, 2] = 1
        zMatrix[2, zNodes+1] = z1

        # Second last row: Modified third derivative term in equation
         cdn3Offset = [ 1 -6 12 -10 3 ] # f-3 to f+1
         zMatrix = addStencil(cdn3Offset, 1, zNodes-1, zNodes-4, zMatrix)
         
         i = zNodes-1
         secondDerMultiple = 3*z[i]
         zMatrix = addCentralDerivative(cdn2, secondDerMultiple, i, zMatrix)
         thirdDerMultiple = -2*( (z[i+1] - z[i-1]) / (2*dx) )
         zMatrix = addCentralDerivative(cdn, thirdDerMultiple, i, zMatrix)
         zMatrix[i, zNodes+1] = -t[i-2]

        # Last row: derivative condition
        cdnR = [ -1 0 1 ]
        zMatrix = addStencil(cdnR, 1, zNodes, zNodes-2, zMatrix)
        
        if debug
            println("Complete Z matrix:")
            printMatrix(zMatrix)
        end

        # Now zMatrix should be complete, solve it
        newZ = Solve_GaussElim!(zMatrix)

        if debug
            println("New Zeta: $newZ")
        end

        #################### Build and solve t Matrix #####################
        # No additional equations required for t, dirichlet on both sides
        tMatrix = zeros(Float64, tNodes, tNodes+1)

        # Apply appropriate stencils for interior nodes
        for i in 2:(tNodes-1)
            tMatrix = addCentralDerivative(cdn2, 1, i, tMatrix)
            firstDerMultiple = 3 * 0.71 * newZ[i+1]
            tMatrix = addCentralDerivative(cdn, firstDerMultiple, i, tMatrix)
        end

        if debug
            println("Interior Stencils Applied:")
            printMatrix(tMatrix)
        end

        # Now apply boundary conditions
        # For the first row, value of T0 is known to be t1
        tMatrix[1,1] = 1
        tMatrix[1, tNodes+1] = t1

        # For the last row, value of T1 is known to be t2
        tMatrix[tNodes, tNodes] = 1
        tMatrix[tNodes, tNodes+1] = t2

        if debug
            println("Last Row Final:")
            printMatrix(tMatrix)
        end

        # Solve
        newT = Solve_GaussElim!(tMatrix)

        if debug
            println("New T: $newT")
        end

        ############## Evaluate dt, dz ##############
        maxDx = 0
        for i in 3:nNodes
            dz = abs(newZ[i] - z[i])
            dt = abs(newT[i-1] - t[i-1])
            maxDx = max(max(dz, dt), maxDx)
        end
       
        println("Iteration $iterationCounter, Max dx = $maxDx")
        
        iterationCounter += 1
        t = newT
        z = newZ

        println("")

        if iterationCounter == 20
            break
        end
    end

    return z, t
end

function solve_Q3_SecondOrder4(nNodes, z1, t1, t2, epsilon=0.0001, debug=false)
    eta1 = 0
    eta2 = 20
    dx = (eta2 - eta1) / (nNodes + 1)
    
    if debug
        println("Nodes = $nNodes, dx = $dx")
    end

    #Create initial value vectors
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

    if debug
        println("")
        println("Initial z-vector: $z")
        println("Initial t-vector: $t")
    end

    # Second order central stencils for each derivative
    cdn3 = (1/(2* dx*dx*dx)) .* [ -1 2 0 -2 1 ]
    cdn2 = (1/(dx*dx)) .* [ 1 -2 1 ]
    cdn = (1/(2*dx)) .* [ -1 0 1 ]

    if debug
        println("")
        println("3rd Derivative Stencil: $cdn3")
        println("2nd Derivative Stencil: $cdn2")
        println("1st Derivative Stencil: $cdn")
        println("")
    end

    maxDx = 1
    iterationCounter = 1
    while maxDx > epsilon
        ################### Build and solve z Matrix ###################
        # One additional equation required on the RHB, where we only have a derivative boundary condition
        zMatrix = zeros(Float64, zNodes+tNodes, zNodes+tNodes+1)

        if debug
            println("Init:")
            printMatrix(zMatrix)
        end

        # Will simply apply the stencils for interior nodes
        for i in 3:(zNodes-2)
            zMatrix = addCentralDerivative(cdn3, 1, i, zMatrix)
            secondDerMultiple = 3*z[i]
            zMatrix = addCentralDerivative(cdn2, secondDerMultiple, i, zMatrix)
            thirdDerMultiple = -2*( (z[i+1] - z[i-1]) / (2*dx) )
            zMatrix = addCentralDerivative(cdn, thirdDerMultiple, i, zMatrix)
            zMatrix[i, zNodes+i-1] = 1
        end

        if debug
            println("Interior Stencils Applied:")
            printMatrix(zMatrix)
        end

        # Now need to apply boundary conditions, make the first two and the last two rows
        # Second row: derivative condition
        cdn = [ -1 0 1 ]
        zMatrix = addStencil(cdn, 1, 1, 1, zMatrix)

        # Second row is equation for first boundary node
        zMatrix[2, 2] = 1
        zMatrix[2, zNodes+tNodes+1] = z1

        # Second last row: Modified third derivative term in equation
         cdn3Offset = [ -1 7 -22 34 -25 7 ] / (2 * dx*dx*dx) # f-3 to f+1
         zMatrix = addStencil(cdn3Offset, 1, zNodes-1, zNodes-5, zMatrix)
         
         i = zNodes-1
         secondDerMultiple = 3*z[i]
         zMatrix = addCentralDerivative(cdn2, secondDerMultiple, i, zMatrix)
         thirdDerMultiple = -2*( (z[i+1] - z[i-1]) / (2*dx) )
         zMatrix = addCentralDerivative(cdn, thirdDerMultiple, i, zMatrix)
         zMatrix[i, zNodes+i-1] = 1

        # Last row: derivative condition
        cdnR = [ 1 -4 3 ]
        zMatrix = addStencil(cdnR, 1, zNodes, zNodes-2, zMatrix)
        
        if debug
            println("Complete Z matrix:")
            printMatrix(zMatrix)
        end

        #################### Build and solve t Matrix #####################

        # Apply appropriate stencils for interior nodes
        for i in 2+zNodes:(zNodes+tNodes-1)
            zMatrix = addCentralDerivative(cdn2, 1, i, zMatrix)
            firstDerMultiple = 3 * 0.71 * z[i+1 - zNodes]
            zMatrix = addCentralDerivative(cdn, firstDerMultiple, i, zMatrix)
        end

        if debug
            println("Interior Stencils Applied:")
            printMatrix(zMatrix)
        end

        # Now apply boundary conditions
        # For the first row, value of T0 is known to be t1
        zMatrix[1+zNodes,1+zNodes] = 1
        zMatrix[1+zNodes, zNodes+tNodes+1] = t1

        # For the last row, value of T1 is known to be t2
        zMatrix[tNodes+zNodes, tNodes+zNodes] = 1
        zMatrix[tNodes+zNodes, zNodes+tNodes+1] = t2

        if debug
            println("Last Row Final:")
            printMatrix(zMatrix)
        end

        # Solve
        newTZ = Solve_GaussElim!(zMatrix)

        newZ = newTZ[1:zNodes]
        newT = newTZ[zNodes+1:zNodes+tNodes]
        if debug
            println("New z: $newZ")
            println("New T: $newT")
        end

        ############## Evaluate dt, dz ##############
        maxDx = 0
        for i in 3:nNodes+1
            dz = abs(newZ[i] - z[i])
            dt = abs(newT[i-1] - t[i-1])
            maxDx = max(max(dz, dt), maxDx)
        end
       
        println("Iteration $iterationCounter, Max dx = $maxDx")
        
        iterationCounter += 1
        # Under-relaxation to stabilize the solution
        t = (newT + t) / 2
        z = (newZ + z) / 2

        println("")

        if iterationCounter == 200
            return z, t
        end
    end

    return z, t
end

# Plot Results
nNodes = 50
z, T = solve_Q3_SecondOrder4(nNodes, 0, 1, 0, 0.001, false)

z = z[2:nNodes+3]
eta = Array{Float64, 1}(undef, nNodes+2)
for i in 1:nNodes+2
    eta[i] = ((i-1)/(nNodes+1)) * 20
end

println(z[nNodes])

plot(eta, z, label="zeta", size=(1720, 880))
plot!(eta, T, label="T")
gui()