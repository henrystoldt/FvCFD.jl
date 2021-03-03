# JuliaCFD

![](https://github.com/henrystoldt/JuliaCFD/workflows/Tests/badge.svg)

<img src="https://raw.githubusercontent.com/henrystoldt/JuliaCFD/master/Resources/fStep.gif?raw=true" alt="Forward step simulation"
  title="MAPLEAF" height=225 style="padding-right: 10px;"/>

JuliaCFD is a simple explicit compressible Euler solver for 3D unstructured polyhedral meshes written in the Julia programming language.
JuliaCFD is compatible with OpenFOAM meshes and outputs .vtu files which can be viewed with most post processing software (ex. Paraview).
The code is compact enough to be read fully by individual users/developers (<2000 lines) and makes for an excellent hands-on introduction to CFD.
This code was originally created as a final project for a graduate numerical methods class, and would make for a good starting point for additional projects or low-overhead numerical research.
Ongoing projects include the implementation of adaptive meshing and implicit time stepping.
Contributions are welcome.
In terms of numerics, JuliaCFD uses the JST or Roe/MUSCL second-order convective schemes and explicit local/global Runge-Kutta time stepping of up to fourth order.
Have a look at or run runSolver.jl to get started!

## Install
TODO: Julia Package

## Sample results

## For newcomers to Julia
Julia is unlike most other languages in that it makes extensive use of Just-in-Time Compilation (JIT).
In scientific computing, this approach can often deliver something approximating the speed of C++ combined with the simplicity of Python.
Unfortunately, since compilation and execution are mixed together, Julia's speed is not always immediately obvious.

### Compilation (Julia user experience)
The first time you execute a piece of code Julia performs type-inference, compiles any code for which it is able to infer types, and then runs it.
As a result, it often takes longer than you would expect to run code for the first time (compared to languages like Python).
Fortunately, after you've run code once, instances of Julia will cache the compiled code and will only need to recompile parts that you change, making subsequent runs much faster.
To take advantage of this, you can start an instance of Julia and run your code in it repeatedly, without closing the REPL.

However, this can also cause problems in that any function that has been executed/compiled in an instance of Julia remains cached and available for use until the REPL is closed.
As such, instances of Julia are stateful.
This can lead to unexpected/non-replicable behavior in your code if it uses an old or deleted function which would not be present in a brand new instance of Julia executing your code.
Check for these problems by restarting Julia and running your code/tests from scratch every so often or after significant changes.

### Type-Stability (Julia developer experience)
Delving into the previous topic a bit more, Julia can execute your code in two distinct ways: compiled or not-compiled/interpreted.
As you can imagine, executing compiled code is orders of magnitude faster than interpreting code and Julia will compile any code it can to take advantage of this.
As such, the first and most important step to writing high performance Julia code is to ensure that Julia can compile it.

To be able to compile code, Julia needs to know the types of all variables involved at all times.
This allows Julia to both allocate the appropriate amount of memory for each variable, and (later) to understand/interpret the data stored there.
Julia can infer most variable types and does not need type annotations everywhere,  with type inference.
Whenever possible, make sure functions return variables of a consistent type (a type-stable function).
Similarly, avoid changing the types of variables within functions - keep variables 'type-stable'.
Use Julia's [@code-warntype macro](https://docs.julialang.org/en/v1/manual/performance-tips/#man-code-warntype) to check for type-stability problems.


A common workflow is to start by a creating an prototype without worrying excessively about type stability or performance, resulting in code similar to Python.
Then, as/if the project progresses, the code's speed can be increased dramatically by adding a few type annotations and minor refactoring, without having to rewrite it in another language.

More information about Julia: https://julialang.org/  
General documentation: https://docs.julialang.org/en/v1/  
Julia performance tips: https://docs.julialang.org/en/v1/manual/performance-tips/
# Nomenclature
All units SI  
For the simple variables:

Var | Meaning
--- | ---
e | internal energy  
P | Pressure  
rho | Density  
T | Temperature  
U | Velocity   

And the state & flux variables:  

Var | Meaning | Definition  
--- | --- | ---  
xMom | x-Momentum | rho*U  
eV2 | total energy | rho*(e + U^2/2)  
rhoU2p | flux of x-momentum | rho*U^2 + P   
rhoUeV2PU | x-direction flux of total energy | U*eV2 + P*U  

# References
These are referred to throughout the code:  
**Moukalled et al.**: The Finite Volume Method in Computational Fluid Dynamics: An Advanced Introduction with OpenFOAM and Matlab  
**Hoffman**: Numerical Methods for Engineers and Scientists  
**Versteeg et al.**: An Introduction to Computational Fluid Dynamics: The Finite Volume Method (2nd Edition)  
**Anderson**: Computational Fluid Dynamics: The Basics with Applications  
