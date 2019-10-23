# ENME 631 Numerical Methods
Name: Henry Stoldt  
ID: 10127324  
Date: Thursday Oct. 24th, 2019  
Language: Julia 

The file "verlet.jl" runs all of code requested for Additional Homework #1. Explanation and comments below.  
Plots and code are attached.

### The Solver:
I've constructed a verlet position integrator for n spatial coordinates. It works based on the process outlined in the additional notes Dr. Mohamad posted on D2l (section 7.5.4, "The Verlet Method").
The algorithm first calculates the position at time i+1 explicitly, using the position, velocity and accelerations from time i.
Next, it calculates the velocity at time i+1 implicitly, in my case using a predictor-corrector method.
The predictor corrector method I implemented predicts the velocity at time i+1 using the explicit first-order euler method, then corrects it using equation 7.52 from the provided notes, where the acceleration at time i+1 is calculated using the predicted values of the velocity.

### Solving the sample problem:
As requested, the sample problem was solved for the time interval 0-10 seconds. Time steps smaller than the suggested value of 0.1 seconds did not noticeably change the solution.