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
function solve_Q3(nNodes)
    eta1 = 0
    eta2 = 20
    dx = (eta2 - eta1) / (nNodes + 1)

    nVars = nNodes * 2

    # Second order central stencils for each derivative
    cdn3 = (1/(2* dx*dx*dx)) .* [ -1 2 0 -2 1 ]
    cdn2 = (1/(dx*dx)) .* [ 1 -2 1 ]
    cdn = (1/(2*dx)) .* [ -1 0 1 ]
    
    # Augmented matrix
    # Order of variables is z1, T1, z2, T2, z3, T3 etc...
    matrix = Array{Float64, 2}(undef, nVars, nVars+1)

    # Build the linear parts of the matrix (first term of both equations)
    for i in 1:nNodes
        centerPt = i
        stencilSemiSpan = floor(size(cdn3) / 2)
        minIndex = max(1, centerPt - stencilSemiSpan)
        maxIndex = min(nVars, )
        for c in range minIndex:
    end


    


end