include("LinearSolvers.jl")
__precompile__()

################# Interactivity ####################
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

    EqnStrings = Array{String, 1}(undef, nEquations)
    Fns = Array{Function, 1}(undef, nEquations)
    maxX = 1
    for i in 1:nEquations
        eqn = readline()
        #eqn2 will be mutilated below searching for x-variables
        eqn2 = eqn

        # Check whether the maxmimum x-number in this equation is higher than the previously stored value
        while occursin(r"x[0-9]+", eqn2)
            #Find x-variable, store number if higher than those previously encountered
            matchLoc = findnext(r"x[0-9]+", eqn2, 1)
            varText = eqn2[matchLoc]
            varText = varText[2:length(varText)]
            varNum = parse(Int32, varText)
            if varNum > maxX
                maxX = varNum
            end

            # Remove the x-variable located
            if matchLoc[1] == 1
                eqn2 = eqn2[matchLoc[2]+1:length(eqn2)]
            else
                eqn2 = eqn2[1:matchLoc[1]-1] * eqn2[matchLoc[2]:length(eqn2)]
            end
                        
        end

        # Store clean equation string for future use
        EqnStrings[i] = eqn
    end

    # Convert all the functions to be functions of the overall highest x
    variableString = "( "
    for i in 1:maxX-1
        variableString = variableString * "x$i, "
    end
    variableString = variableString * "x$maxX ) -> "

    # println("Functions:")
    for i in 1:nEquations
        fnString = variableString * EqnStrings[i]
        # println("Function $i: $fnString")
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

################# Nonlinear functions - used only for Assignment 2 #######################
function f1(x1, x2, x3)
    return exp(2*x1) - x2 + 4
end
function f2(x1, x2, x3)
    x2 - x3*x3 - 1
end
function f3(x1, x2, x3)
    x3 - sin(x1)
end

################## Nonlinear Newton-Raphson Solver ###################
#Function to calculate partial derivatives numerically using 2nd order central method
#Pass in arrays of functions and x-values at which to calculate derivatives
function ddx(fns, x::Array{Float64}, epsilon=0.00000001)
    nEqns = size(fns, 1)
    nXs = size(x, 1)

    ddxs = Array{Float64, 2}(undef, nEqns, nXs)
    for i in 1:nEqns
        for a in 1:nXs
            xCall1 = copy(x)
            xCall2 = copy(x)

            # Modify each copy of the initial x-vector - in one move x[a] slightly forward, in the other slightly backward
            xCall1[a] = x[a] + epsilon
            xCall2[a] = x[a] - epsilon
            # Central Method
            ddxs[i,a] = (fns[i](xCall1...) - fns[i](xCall2...)) / (2*epsilon)
        end
    end
    return ddxs
end

# Pass in array of functions and array of x's, simply evaluates them at x
function calcResiduals(fns, x)
    nFns = size(fns, 1)

    res = Array{Float64, 1}(undef, nFns)
    for i in 1:nFns
        res[i] = fns[i](x...)
    end

    return res
end

# Function to solve set of nonlinear equations using Newton-Raphson method, calculate x-vector
    # Solves [Adx = -residuals] matrix using Gauss Elimination, iterates until convergence or iterLimit is reached
function solve_NonLinear!(fns, xInit, iterLimit=100; printResiduals=true, residualTolerance=0.00001)
    if printResiduals
        println("Nonlinear Newton-Raphson Solver:")
    end
    nFns = size(fns, 1)
    nXs = size(xInit, 1)
    AugmentedMatrix = Array{Float64, 2}(undef, nFns, nXs + 1)
    
    maxResidual = 1000
    iterationCounter = 1
    while maxResidual > residualTolerance && iterationCounter <= iterLimit
        partialDerivatives = ddx(fns, xInit)        
        residuals = calcResiduals(fns, xInit)
        maxResidual = maximum(abs.(residuals))

        if printResiduals
            println("Iteration $iterationCounter: maxRes = $maxResidual")
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
        
        dx = Solve_GaussElim!(AugmentedMatrix)

        xInit += dx
        iterationCounter += 1
    end

    if iterationCounter > iterLimit
        println("ERROR: Convergence not achieved in $iterLimit iterations")
    end

    return xInit
end