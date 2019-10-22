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

function dx2dt(vx, vy)
    # -a/m * sqrt(vx^2 + vy^2)*vx - w*b*vy/m
    return -5 #-0.1*sqrt(vx*vx + vy*vy)*vx - 2*vy
end

function dy2dt(vx, vy)
    # -g -a/m*sqrt(vx^2 + vy^2)vy + w*b*vx/m
    return -10 #-9.81 - 0.1*sqrt(vx*vx + vy*vy)*vy + 2*vx
end

function solve_Verlet(x1, y1, vx1, vy1, timeStep, endTime)
    nTimeSteps = ceil(Int32, endTime/timeStep)
    results = Array{Float64, 2}(undef, 7, nTimeSteps+1)
    
    iterationCounter = 1
    results[1,iterationCounter] = 0
    results[2,iterationCounter] = x1
    results[3,iterationCounter] = y1
    results[4,iterationCounter] = vx1
    results[5,iterationCounter] = vy1
    results[6,iterationCounter] = dx2dt(vx1, vy1)
    results[7,iterationCounter] = dy2dt(vx1, vy1)

    # Compute second set of values to start the method
    iterationCounter += 1
    currTime = timeStep
    x2 = x1 + vx1*timeStep + 0.5*results[6,1]*timeStep*timeStep
    y2 = y1 + vy1*timeStep + 0.5*results[7,1]*timeStep*timeStep
    predVx = vx1 + results[6,1]*timeStep
    predVy = vy1 + results[7,1]*timeStep
    ax2 = dx2dt(predVx, predVy)
    ay2 = dy2dt(predVx, predVy)
    corrVx = vx1 + (results[6,1] + ax2)*0.5*timeStep
    corrVy = vy1 + (results[7,1] + ay2)*0.5*timeStep

    results[1,iterationCounter] = timeStep
    results[2,iterationCounter] = x2
    results[3,iterationCounter] = y2
    results[4,iterationCounter] = corrVx
    results[5,iterationCounter] = corrVy
    results[6,iterationCounter] = ax2
    results[7,iterationCounter] = ay2

    while currTime + 0.00000001 < endTime
        # Solve for new accelerations and positions
        ax = dx2dt(results[4, iterationCounter], results[5, iterationCounter])
        newX = results[2, iterationCounter] + results[4, iterationCounter]*timeStep + 0.5*ax*timeStep*timeStep

        ay = dy2dt(results[4, iterationCounter], results[5, iterationCounter])
        newY = results[3, iterationCounter] + results[5, iterationCounter]*timeStep + 0.5*ay*timeStep*timeStep

        # Solve for v @ i+1 using predictor-corrector method
        # Predictor
        lastVx = results[4, iterationCounter]
        predNewVx = lastVx + ax*timeStep
        lastVy = results[5, iterationCounter]
        predNewVy = lastVy + ay*timeStep
        
        # Corrector
        newAx = dx2dt(predNewVx, predNewVy)
        newAy = dy2dt(predNewVx, predNewVy)
        newVx = lastVx + (ax + newAx)*timeStep/2
        newVy = lastVy + (ay + newAy)*timeStep/2

        currTime += timeStep
        iterationCounter += 1
        
        results[1, iterationCounter] = currTime
        results[2, iterationCounter] = newX
        results[3, iterationCounter] = newY
        results[4, iterationCounter] = newVx
        results[5, iterationCounter] = newVy
        results[6, iterationCounter] = ax
        results[7, iterationCounter] = ay        
    end

    return results
end

results = verlet(0, 10, 5, 10, 0.01, 10)
time = results[1,:]
x = results[2, :]
y = results[3, :]
Vx = results[4, :]
Vy = results[5, :]
ax = results[6, :]
ay = results[7, :]

xPlot = plot(time, x, label="X Position", title="Position vs Time", xlabel="t")
xPlot = plot!(time, y, label="Y Position")
vPlot = plot(time, Vx, label="X-Velocity", title="Velocity vs Time", xlabel="t")
vPlot = plot!(time, Vy, label="Y-Velocity")
aPlot = plot(time, ax, label="X-Accel", title="Acceleration vs Time", xlabel="t")
aPlot = plot!(time, ay, label="Y-Accel")
plot(xPlot, vPlot, aPlot, size=(1720, 880), window_title="Verlet R")
gui()

# Solve with explicit Euler to check
# init = [5,10]
# dt = 0.01
# timeSteps = 1000
# ddts = [ dx2dt, dy2dt ]
# results = solve_ExplicitEuler(init, ddts, dt, timeSteps)
# Vx = results[1, :]
# Vy = results[2, :]
# time = results[3,:]
# vPlot = plot(time, Vx, label="X-Velocity", title="Velocity vs Time", xlabel="t", size=(860, 880), window_title="Verlet R")
# vPlot = plot!(time, Vy, label="Y-Velocity")
# gui()
