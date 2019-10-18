inlude("LinearSolvers.jl")

# zeta will be referred to as z or Z

function residual1 (d3Zdn, d2Zdn, dZdn, T, Z)
    return d3Zdn + 3*Z*d2Zdn - 2 *(dZdn*dZdn) + T
end

function residual2 (d2Tdn, dTdn, Z)
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
        x[i] == initZ(eta)
        x[i+1] == initT(eta)
    end

    # Second order central stencils for each derivative
    cdn3 = (1/(2* dx*dx*dx)) .* [ -1 2 0 -2 1 ]
    cdn2 = (1/(dx*dx)) .* [ 1 -2 1 ]
    cdn = (1/(2*dx)) .* [ -1 0 1 ]
    
    # Build the linear parts of the matrix (first term of both equations)
    # Make a single stencil, replicate it for each row, adjust the first and last for the boundary conditions
    bandwidth = size(cdn3) * 2
    
    # Order of variables is z1, T1, z2, T2, z3, T3 etc...
    thomasZRow = Array{Float64, 1}(0, bandwidth)
    for i in 1:size(cdn3)
        thomasZRow[2*i-1] += cdn3[i]
    end
    Zcenter = floor((size(cdn3) + 1) / 2)
    thomasZRow[Zcenter] += 1

    thomasTRow = Array{Float64, 1}(0, bandwidth)
    offset = size(cdn3) - size(cdn2)
    for i in 1+offset:size(cdn2)+offset
        thomasTRow[2*i] += cdn2[i]
    end

    matrix = Array{Float64, 2}(0, nVars, nVars+1)
    for r in 1:nNodes
        matrix[r, 1:nVars] += thomasZRow
        matrix[r+1, 1:nVars] += thomasTRow
    end
        

    while maxDx > epsilon
        # Create matrix by substituting in nonlinear values from last iteration
        tZRow = copy(thomasZRow)
        
        offset = size(cdn2) - size(cdn3)
        cdn2_Z2 = 
        for i in 1+offset:size(cdn2)+offset
        
        tTRow = copy(thomasTRow)

    end


end