include("LinearSolvers.jl")
include("NonLinearSolvers.jl")

# Parameters:
    # vars should provide initial conditions for x1, ..., xn
    # ddts should provide d/dt functions for x1, ... , xn
        # parameters for each ddt function must be x1, ..., xn
    # dt = marching variable step size
    # nTimeSteps = Number of time steps to complete

# Return Value:
    # Functions return a 2D array
    # Columns correspond to marching steps, with the first column containing initial values of each variable x1 through xn    
    # The last row contains the marching variable value for each step

function solve_ExplicitEuler(vars, ddts, dt, nTimeSteps)
    #Setup
    nVars = length(vars)
    timeStepCounter = 2
    results = Array{Float64, 2}(undef, nVars+1, nTimeSteps+1)

    for i in 1:nVars
        results[i,1] = vars[i]
    end
    results[nVars+1,1] = 0

    while timeStepCounter <= nTimeSteps+1
        for i in 1:nVars
            results[i, timeStepCounter] = results[i, timeStepCounter-1] + dt*ddts[i](results[1:nVars, timeStepCounter-1]...)
        end
        results[nVars+1, timeStepCounter] = results[nVars+1, timeStepCounter-1] + dt

        timeStepCounter += 1
    end

    return results
end

# Uses Newton-Raphson iterations, not predictor-corrector method
function solve_ImplicitEuler(vars, ddts, dt, nTimeSteps)
    #Setup
    nVars = length(vars)
    timeStepCounter = 2
    results = Array{Float64, 2}(undef, nVars+1, nTimeSteps+1)

    for i in 1:nVars
        results[i,1] = vars[i]
    end
    results[nVars+1,1] = 0

    augmentedMatrix = Array{Float64, 2}(undef, length(ddts), nVars+1)

    while timeStepCounter <= nTimeSteps+1
        #Predict solution with explicit Euler
        initGuesses = Array{Float64, 1}(undef, nVars)
        for i in 1:nVars
            initGuesses[i] = results[i, timeStepCounter-1] + dt*ddts[i](results[1:nVars, timeStepCounter-1]...)
        end
    
        #Refine solution using Newton-Raphson iterations to get backward Euler solution
        residualTolerance = 0.00001
        maxRes = 1
        maxIter = 100
        iterCount = 0
        while maxRes > residualTolerance && iterCount < maxIter
            augmentedMatrix[:,1:nVars] = ddx(ddts, initGuesses)*dt

            for i in 1:nVars
                augmentedMatrix[i, i] -= 1
            end

            maxRes = 0
            for i in 1:nVars
                residual = -(results[i, timeStepCounter-1] + ddts[i](initGuesses...)*dt - initGuesses[i])
                if abs(residual) > maxRes
                    maxRes = abs(residual)
                end
                augmentedMatrix[i, nVars+1] = residual
            end
            initGuesses += Solve_GaussElim!(augmentedMatrix)
            iterCount += 1
        end

        if iterCount == maxIter
            println("ERROR: timestep $timeStepCounter not converged, maxRes = $maxRes")
        end

        results[1:nVars, timeStepCounter] = initGuesses
        results[nVars+1, timeStepCounter] = results[nVars+1, timeStepCounter-1] + dt

        timeStepCounter += 1
    end

    return results
end

# Uses predictor-corrector method
function solve_ModifiedEuler(vars, ddts, dt, nTimeSteps)
    #Setup
    nVars = length(vars)
    timeStepCounter = 2
    results = Array{Float64, 2}(undef, nVars+1, nTimeSteps+1)

    for i in 1:nVars
        results[i,1] = vars[i]
    end
    results[nVars+1,1] = 0

    predictor1 = Array{Float64, 1}(undef, nVars)
    predictor2 = Array{Float64, 1}(undef, nVars)
    while timeStepCounter <= nTimeSteps+1
        for i in 1:nVars
            predictor1[i] = results[i, timeStepCounter-1] + dt*ddts[i](results[1:nVars, timeStepCounter-1]...)
        end
        
        for i in 1:nVars
            predictor2[i] = results[i, timeStepCounter-1] + dt*ddts[i](predictor1...)
            results[i, timeStepCounter] = (predictor1[i] + predictor2[i]) / 2
        end
        
        results[nVars+1, timeStepCounter] = results[nVars+1, timeStepCounter-1] + dt
        timeStepCounter += 1
    end

    return results
end

function solve_RK4(vars, ddts, dt, nTimeSteps)
    #Setup
    nVars = length(vars)
    timeStepCounter = 2
    results = Array{Float64, 2}(undef, nVars+1, nTimeSteps+1)

    for i in 1:nVars
        results[i, 1] = vars[i]
    end
    results[nVars+1, 1] = 0

    while timeStepCounter <= nTimeSteps+1
        dy1 = Array{Float64, 1}(undef, nVars)
        for i in 1:nVars
            dy1[i] = ddts[i](results[1:nVars, timeStepCounter-1]...) * dt
        end
        
        y2 = results[1:3, timeStepCounter-1] + dy1/2
        dy2 = Array{Float64, 1}(undef, nVars)
        for i in 1:nVars
            dy2[i] = ddts[i](y2...) * dt
        end

        y3 = results[1:3, timeStepCounter-1] + dy2/2
        dy3 = Array{Float64, 1}(undef, nVars)
        for i in 1:nVars
            dy3[i] = ddts[i](y3...) * dt
        end

        y4 = results[1:3, timeStepCounter-1] + dy3
        dy4 = Array{Float64, 1}(undef, nVars)
        for i in 1:nVars
            dy4[i] = ddts[i](y4...) * dt
        end

        for i in 1:nVars
            results[i, timeStepCounter] = results[i, timeStepCounter-1] + (1/6)*(dy1[i] + 2*dy2[i] + 2*dy3[i] + dy4[i])
        end

        results[nVars+1, timeStepCounter] = results[nVars+1, timeStepCounter-1] + dt
        timeStepCounter += 1
    end

    return results
end
