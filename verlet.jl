include("ODESolvers.jl")
using Plots
using Plots.PlotMeasures
plotly()

############## Specific to Extra HW Assignment 1 ###############
# a = 0.1
# m = 1
# w = 2
# g = 9.81
# b = 1

function dx2dt(initXs, initVs)
    # -a/m * sqrt(vx^2 + vy^2)*vx - w*b*vy/m
    vx = initVs[1]
    vy = initVs[2]
    return -0.1*sqrt(vx*vx + vy*vy)*vx - 2*vy
end

function dy2dt(initXs, initVs)
    # -g -a/m*sqrt(vx^2 + vy^2)vy + w*b*vx/m
    vx = initVs[1]
    vy = initVs[2]
    return -9.81 - 0.1*sqrt(vx*vx + vy*vy)*vy + 2*vx
end

########### Verlet Motion Solver #############
# Pass in arrays of initial positions, velocities, and acceleration functions
    # Acceleration functions must accept the positions and velocities as arguments
# Returns results matrix, with rows containing the folowing, with N coordinate directions
    # Time values
    # X1 Position
    # ...
    # XN Position
    # X1 Velocity
    # ...
    # XN Velocity
    # X1 Accel
    # ...
    # XN Accel
function solve_Verlet(initXs, initVs, accelFns, timeStep, endTime)
    nTimeSteps = ceil(Int32, endTime/timeStep)

    nPos = size(initXs, 1)
    nResultsVars = nPos*3 + 1
    results = Array{Float64, 2}(undef, nResultsVars, nTimeSteps+1)
    
    iterationCounter = 1
    results[1,1] = 0

    # Initialize results matrix with initial values
    for i in 1:nPos
        results[i+1, 1] = initXs[i]
        results[i+1+nPos, 1] = initVs[i]
        results[i+1+2*nPos, 1] = accelFns[i](initXs, initVs)
    end

    ############ Utility Functions #############
    function vel(var, iteration) 
        return results[var+1+nPos, iteration]
    end

    function pos(var, iteration)
        return results[var+1, iteration]
    end

    function accel(var, iteration)
        return accelFns[var](getPosVector(iteration), getVelVector(iteration))
    end

    function getPosVector(iteration)
        pstns = []
        for i in 1:nPos
            push!(pstns, results[i+1, iteration])
        end
        return pstns
    end

    function getVelVector(iteration)
        vels = []
        for i in 1:nPos
            push!(vels, results[i+nPos+1, iteration])
        end
        return vels
    end

    function storePos(var, iteration, pos)
        results[var+1, iteration] = pos
    end

    function storeVel(var, iteration, vel)
        results[var+nPos+1, iteration] = vel
    end

    function storeAccel(var, iteration, accel)
        results[var + 2*nPos + 1, iteration] = accel
    end

    ########## Verlet integration loop #############
    currTime = 0
    while currTime + 0.00000001 < endTime
        # Solve for new accelerations and positions, velocity predictor
        predVels = []
        for i in 1:nPos
            ax = accel(i, iterationCounter)
            newX = pos(i, iterationCounter) + vel(i, iterationCounter)*timeStep + 0.5*ax*timeStep*timeStep
            push!(predVels, vel(i, iterationCounter) + ax*timeStep)
    
            storeAccel(i, iterationCounter, ax)
            storePos(i, iterationCounter+1, newX)
        end
        
        # Velocity corrector loop
        for i in 1:nPos
            ax = results[2*nPos + 1 + i, iterationCounter]
            ax2 = accelFns[i](getPosVector(iterationCounter+1), predVels)
            corrVx = vel(i, iterationCounter) + (ax + ax2)*timeStep/2
            
            storeVel(i, iterationCounter+1, corrVx)
        end

        currTime += timeStep
        results[1, iterationCounter] = currTime
        iterationCounter += 1      
    end

    return results
end

########### Call integration function, plot results ##############
initPos = [0, 10]
initVels = [5, 10]
accelFns = [dx2dt, dy2dt]
dt = 0.1
endTime = 10
results = solve_Verlet(initPos, initVels, accelFns, dt, endTime)
nTimeSteps = round(Int32, endTime/dt)
time = results[1,1:nTimeSteps]
x = results[2, 1:nTimeSteps]
y = results[3, 1:nTimeSteps]
Vx = results[4, 1:nTimeSteps]
Vy = results[5, 1:nTimeSteps]

xPlot = plot(time, x, label="X Position", title="Position vs Time", xlabel="t")
xPlot = plot!(time, y, label="Y Position")
vPlot = plot(time, Vx, label="X-Velocity", title="Velocity vs Time", xlabel="t")
vPlot = plot!(time, Vy, label="Y-Velocity")
plot(xPlot, vPlot, aPlot, size=(1720, 880), window_title="Verlet")
gui()
