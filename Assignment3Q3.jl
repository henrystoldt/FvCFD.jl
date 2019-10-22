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
    return eta*eta
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
            return (4*x[nVars-1] - x[nVars-3]) / 3
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

# Plot Results
nNodes = 100
x = solve_Q3_SecondOrder(nNodes, 0.002)
z = Array{Float64, 1}(undef, nNodes)
T = Array{Float64, 1}(undef, nNodes)
eta = Array{Float64, 1}(undef, nNodes)
for i in 1:nNodes
    z[i] = x[2*i-1]
    T[i] = x[2*i]
    eta[i] = (i/nNodes) * 20
end

plot(eta, z, label="zeta")
plot!(eta, T, label="T")
gui()

function solve_Q3_Shooting(nNodes, epsilon=0.001)
    # Need to assume two boundary conditions
    # dz/dn = y
    # dy/dn = x
    # dx/dn = -3zx +2Y^2 - T

    # dT/dn = u
    # du/dn = -3Pr*z*u

    # T init = 1
    # zeta init = ?
    # y init = 0
    # x init = ?
    # u init = ?

end

function verlet(x1, y1, vx1, vy1, timeStep, endTime)
    currTime = 0
    x = x1
    y = y1
    vx = vx1
    vy = vy1
    while currTime < endTime
        

        currTime += timeStep
end