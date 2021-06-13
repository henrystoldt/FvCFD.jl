---
title: 'FvCFD.jl: A Julia package for high-speed fluid dynamics'
tags:
  - Julia
  - CFD
  - Simulation
  - Finite volume
  - Compressible
authors:
  - name: Henry Stoldt
    affiliation: 1
  - name: Ben Dalman
    affiliation: 1
  - name: Craig Johansen
    affiliation: 1
affiliations:
 - name: Department of Mechanical and Manufacturing Engineering, University of Calgary
   index: 1
date: 6 May 2021
bibliography: paper.bib
---

# Summary

Computational fluid dynamics (CFD) solvers are used in many engineering disciplines to estimate how fluids will flow around or interact with objects or each other.
Compared to experimental studies, they are very flexible in that they can quickly adapt to changing conditions or changing designs while also being relatively inexpensive.
CFD simulations can be conducted at various levels of accuracy, with the cost of increased accuracy usually coming in the form of additional computational cost.
When combined with a suitable wall friction correction, solutions to the compressible Euler equations offer a useful compromise between accuracy and computational cost for many high-speed aerospace applications.
This is especially the case in the context of design optimizations, where many simulations must be performed and computational cost becomes very important.

# Statement of need
`FvCFD.jl` is a compact explicit CFD solver suited for the simulation of high-speed compressible flows while neglecting heat conduction and viscosity (solving the compressible Euler equations).
While there are more scalable tools for industrial and production use, this solver is well suited for learning about CFD and experimenting with numerical methods, and technologically well positioned for future expansion.
In this vein, `FvCFD.jl` is compatible with 3D unstructured OpenFOAM meshes and many common post-processing tools (.vtk output), so it can be substituted for one of the many OpenFOAM solvers in simulation workflows when moving to/from a large-scale applications.

The code's main advantage are that it is very compact (<2000 lines) and could easily be fully read by individual students/developers in a matter of days, yet it is fully functional CFD solver.
This compactness is partially a result of using the Julia programming language, which is particularly good at extracting high performance from simple and general code [@Bezanson:2017].
This is in contrast to many existing open-source CFD solvers (e.g., OpenFOAM, SU2) which due to their comprehensive feature set and corresponding degree of code generalization, have very large and complex code bases.
In addition to educational applications, `FvCFD.jl` provides a convenient starting point for low-overhead numerical methods research, particularly for new spatial and temporal discretizations.
Technologically, since the code is written entirely in Julia, there are ample opportunities for future integration with other Julia packages such as `DifferentialEquations.jl` [@Rackauckas:2017] or `ForwardDiff.jl` [@Revels:2016].
As a compressible flow solver, `FvCFD.jl` is also introducing a new capability to the open-source Julia ecosystem, which has thus far been focused on incompressible flow, through the excellent `Oceananigans.jl` package [@Ramadhan:2020].

`FvCFD.jl` simulations are straightforward to run, requiring a mesh, boundary conditions, information about the fluid, and time stepping parameters.
In terms of numerical schemes, `FvCFD.jl` implements the JST \& MUSCL + Roe convective schemes, the Green-Gauss \& weighted Least-Squares gradient schemes, explicit Runge-Kutta time stepping of orders 1-4 for unsteady cases, and first-order local time stepping for steady state cases.
Sample `FvCFD.jl` results for a supersonic delta wing and Woodward \& Collela's Mach 3 forward step case [@Woodward1984] are shown below in \autoref{fig:deltaWing} and \autoref{fig:forwardStep}, respectively.

![Cut through a FvCFD.jl supersonic delta wing simulation. Showing the portion of the flow field in which the pressure is higher than ambient. Wing (shown in green) is sharp-edged with 67.5 degree leading-edge sweep, flying at Mach 2 at a twelve degree angle of attack. Corresponds to wind-tunnel experiments run by Miller \& Wood [@Miller1985].\label{fig:deltaWing}](DeltaWingImage.png)

<!-- TODO: New, labeled forward step figure and performance comparison to rhoCentralFoam (current picture is very old)  -->
![Comparison of density contour results obtained using OpenFOAM's rhoCentralFoam solver (left) [@Greenshields2010] and FvCFD.jl (right) for the Mach 3 forward step benchmark case.\label{fig:forwardStep}](ForwardStepComparison2.png)

Going forward, the solver's explicit nature combined with Julia's parallelization features have the potential to close the performance gap with existing solvers in large-scale simulations and make this a highly scalable code.
Ongoing projects include the implementation of adaptive meshing and implicit time-stepping.
All questions and contributions are welcome.

# References