include("LinearSolvers.jl")

function getFunctions()
    println("Enter the following:
    # Eqns
    Fi(x1, x2, ..., xn)
    
    Ex:
    2
    x1 + x2 + 3
    2*x1 - x2
    
    Notes:
    Each equation is implicitly = 0
    Each equation must be valid Julia code, variables = x1, x2, ... xn")

    nEquations = chomp(readline())
    nEquations = parse(Int64, nEquations)

    Fns = Array{Function, 1}(undef, nEquations)
    for i in 1:nEquations
        eqn = readline()
        highestX = 1
        fnString = "("
        while occursin("x$highestX", eqn)
            fnString = string(fnString, "x$highestX, ")
            highestX += 1
        end
        len = length(fnString)
        fnString = fnString[1:len-2]
        fnString = string(fnString, ") -> $eqn")
        println("Function $i: $fnString")
        Fns[i] = eval(Meta.parse(fnString))
    end

    return Fns
end

function getInitEstimates()
    println("Enter space-separated initial guesses for x1...xn")
    xInit = chomp(readline())
    xInit = split(xInit)
    xI = Array{Float64, 1}(undef, length(xInit))
    for i in 1:length(xInit)
        xI[i] = parse(Float64, xInit[i])
    end
    return xI
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

#Function to calculate partial derivatives numerically using 2nd order central method
#Pass in arrays of functiosn and x-values
function ddx(fns, x::Array{Float64}, epsilon=0.00000001)
    nEqns = size(fns, 1)
    nXs = size(x, 1)

    ddxs = Array{Float64, 2}(undef, nEqns, nXs)
    for i in 1:nEqns
        for a in 1:nXs
            xCall1 = copy(x)
            xCall2 = copy(x)
            xCall1[a] = x[a] + epsilon
            xCall2[a] = x[a] - epsilon
            ddxs[i,a] = (fns[i](xCall1...) - fns[i](xCall2...)) / (2*epsilon)
        end
    end
    return ddxs
end

# Pass in array of functions and array of x's
function calcResiduals(fns, x)
    nFns = size(fns, 1)

    res = Array{Float64, 1}(undef, nFns)
    for i in 1:nFns
        res[i] = fns[i](x...)
    end

    return res
end

#Function to solve matrix, calculate new x-vector
function solve_NonLinear!(fns, xInit, iterLimit=100; printResiduals=false)
    if printResiduals
        println("Nonlinear Newton-Raphson Solver:")
    end
    nFns = size(fns, 1)
    nXs = size(xInit, 1)
    AugmentedMatrix = Array{Float64, 2}(undef, nFns, nXs + 1)
    
    maxResidual = 1000
    iterationCounter = 1
    while maxResidual > 0.00001 && iterationCounter <= iterLimit
        partialDerivatives = ddx(fns, xInit)
        
        # println("Residuals")
        residuals = calcResiduals(fns, xInit)
        # println(residuals)

        maxResidual = maximum(abs.(residuals))
        if printResiduals
            println("Iteration $iterationCounter, maxRes = $maxResidual")
        end

        #Combine partial derivatives and residuals to make augmented matrix
        for i in 1:nFns
            for a in 1:nXs + 1
                if a <= nXs
                    AugmentedMatrix[i,a] = partialDerivatives[i,a]
                else
                    AugmentedMatrix[i,a] = -1 * residuals[i]
                end
            end
        end

        # println("Augmented Matrix:")
        # printMatrix(AugmentedMatrix)
        
        dx = Solve_GaussElim!(AugmentedMatrix)
        # println("Dx:")
        # println(dx)

        xInit += dx
        iterationCounter += 1
    end

    if iterationCounter == iterLimit
        println("Error: Convergence not achieved in 1000 iterations")
    end

    return xInit
end