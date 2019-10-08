function getFunction()
    println("Enter the following:
    n
    Eqn i (x1, x2, ..., xn)
    Eqn n
    
    Example:
    2
    x1 + x2
    2*x1 - x2
    
    Notes:
    Each equation is implicitly = 0
    Each eqaution must be written as valid Julia code, variables must be called x1, x2, ... xn")

    nEquations = chomp(readline())
    nEquations = parse(Int64, nRows)

    eqns = Array{String, 1}(undef, nEquations)
    for i in 1:nEquations
        eqns[i] = readline()
    end

    #TODO turn strings into equations with Meta.parse

    return eqns
end

# Rows of the nonlinear function
function f1(x1, x2, x3)
    return exp(2*x1) - x2 + 4
end

function f2(x1, x2, x3)
    x2 - x3*x3 - 1
end

function f3(x1, x2, x3)
    x3 - sin(x1)
end

#Function to calculate partial derivatives numerically
#Pass in arrays of functiosn and x-values
function ddx(fns, x::Array{Float64}, epsilon=0.00000001)
    nEqns = size(fns, 1)
    nXs = size(x, 1)

    ddxs = Array{Float64, 2}(undef, nEqns, nXs)
    for i in 1:nEqns
        for a in 1:nXs
            xCall = copy(x)
            xCall[a] = x[a] + epsilon
            ddxs[i,a] = (fns[i](xCall...) - fns[i](x...)) / epsilon
        end
    end
    return ddxs
end

# Pass in array of functions and array of x's
# Fn to calculate residual vector
function calcResiduals(fns, x)
    nFns = size(fns, 1)
    nXs = size(x, 1)
    if nFns != nXs
        throw(ArgumentError("Number of functions should equal number of Xs"))
    end

    res = Array{Float64, 1}(undef, )
    for i in 1:nEqns
        res[i] = fns[i](x...)
    end

    return res
end

#Function to solve matrix, calculate new x-vector

#Function to evaluate arbitrary functions which are inputted?