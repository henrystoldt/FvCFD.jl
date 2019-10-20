include("LinearSolvers.jl")

# zeta will be referred to as z or Z

function residual1(d3Zdn, d2Zdn, dZdn, T, Z)
    return d3Zdn + 3*Z*d2Zdn - 2 *(dZdn*dZdn) + T
end

function residual2(d2Tdn, dTdn, Z)
    return d2Tdn + 3*0.71*Z*dTdn
end

# Initial value functions for T and eta(Z)
function initZ(eta)
    return 0
end

function initT(eta)
    return 1 - eta/20
end

#nNodes = number of interior nodes
function solve_Q3(nNodes, epsilon=0.0001)
    eta1 = 0
    eta2 = 20
    dx = (eta2 - eta1) / (nNodes + 1)

    nVars = nNodes * 2

    #Create initial value vectors
    x = Array{Float64, 1}(undef, nVars)
    for i in 1:nNodes
        eta = dx*i
        x[2*i - 1] == initZ(eta)
        x[2*i] == initT(eta)
    end

    # Second order central stencils for each derivative
    cdn3 = (1/(2* dx*dx*dx)) .* [ -1 2 0 -2 1 ]
    cdn2 = (1/(dx*dx)) .* [ 1 -2 1 ]
    cdn = (1/(2*dx)) .* [ -1 0 1 ]
    
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
    thomasZRow[Zcenter] += 1

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

    # Create matrix, alternating T and zeta equations
    matrix = Array{Float64, 2}(undef, nVars, nVars+1)
    semiSpan = floor(Int32, (bandwidth - 2) / 2)
    for r in 1:nNodes
        # Copy numbers in to matrix
        for colIndex in 1:bandwidth
            column = 2*r - 1 + colIndex - Zcenter
            if column < 1
                continue
            elseif column > nVars
                break
            end
            println("Row: $r, Column: $column")
            matrix[2*r - 1, column] = thomasZRow[colIndex]
            matrix[2*r, column] = thomasTRow[colIndex]
        end
    end

    # Add linear boundary conditions!, matbe set the other augmented terms to zero

    println("Linear terms matrix:")
    printMatrix(matrix)
        
    maxDx = 1
    while maxDx > epsilon
        # Create matrix by substituting in nonlinear values from last iteration
        matrixI = copy(matrix)

        # Calculate values of nonlinear terms for this iteration, lagging the lowest order terms in each nonlinear term
        for r in 1:nNodes
            # Adding second derivative of zeta
            cdn2_Z2 = cdn2 .* 3 * x[2*r - 1]
            # Adding squared first derivative of zeta
            cdn1_Z1 = cdn .* 2 * (x[2*r + 1] - x[2*r - 3]) / (2*dx)
            semiSpan = floor(Int32, size(cdn2, 2) -1 /2)
            for a in 1:size(cdn2, 2)
                offCenter = a - (semiSpan + 1)
                col = r + offCenter
                if col < 1
                    continue
                elseif col > nVars
                    break
                end
                matrix[2*i-1, col] += cdn2_Z2[a]
                matrix[2*i-1, col] += cdn1_Z1[a]
            end

            printMatrix(matrixI)

            #Adding first derivative of T1
            cdn_T1 = cdn .* 3 * 0.71 * x[2*i]
            offset = size(cdn, 2) - size(cdn3, 2)
            for a in 1+offset:size(cdn, 2)+offset
                matrix[2*i, 2*i] += cdn_T1[a]
            end

        end

        printMatrix(matrix)
        # Adjust first and last rows with boundary conditions


        # Solve matrix


        # Calculate residuals/dxs

        break
    end
end

solve_Q3(3)