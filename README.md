# JuliaCFD

# Plan:
1. Euler 1D Solver (shock tube)
2. 3D Adjoint multigridding implicit coupled runge-kutta


# Required Components
1. Parser for SU2 Mesh files
 - Inputs: mesh file name (path)
 - Outputs: 
	- See meshDataStructureDefinition.md

2. Need a standard input file format and a parser for it, which includes
	- Mesh file path
	- Freestream/initialization conditions
	- Boundary conditions
	- FD/FVM (if we have time to implement both?)
	- Numerical convective scheme (MacCormack, Lax-Wendroff, Lax-Friedrich, JST) 

3. Each numerical convective scheme will need a function to setup the solution matrices, and a method to deal with boundary conditions

4. Initialize solution (apply freestream conditions to every DOF) - primitives of rho, u,v,w (depending on dimensions), and e
	-Shock tube needs better ICs - from analytical??

5. Apply convective numerical scheme to take one time-step. Setup matrix equations with unknowns

6. Matrix solver (SOR or Thomas)

7. Convert conserved quatities back into primitives

8. Add results of matrix solution to current state to create state t+dt. Depends on time-marching scheme (euler-impl, euler-exp, rk4)

9. Calculate residuals (on primitives or conserved??)

9a. Print values to console

10. Check convergence criteria

11. Repeat if required

12. Output file format - primitive values at every cell
	- Make vtk file if possible, csv if not
	- If transient, should output primitives for each cell at each time step

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
