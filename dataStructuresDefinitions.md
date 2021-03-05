# Unstructured FVM Mesh Format Definition
Immutable struct  
Passed to nearly all FVM functions  
Defined as follows:  
```julia
mesh =
[
   cells,           #(list of lists, each containing the indices of the faces that make up the cell)  
   cVols,           #(list of scalar cell volumes)  
   cCenters,        #(list of position vectors of cell centers)  
   cellSizes,       #(Matrix: rows are cells, columns are x-, y-, and z-direction sizes of the cell)
   faces,           #(list of lists, each sublist containing two cell indices: the owner cell and the neighbour cell)  
   fAVecs,          #(list of face area vectors)  
   fCenters,        #(list of position vectors of face centers)  
   boundaryFaces    #(list of lists, each containing the indices of faces on the ith boundary)  
]
```

## Notes
- Cells and faces are numbered according to their storage location in the mesh arrays
- Faces are numbered such that boundary faces come last
- Face area vectors point outward from the owner cells, into the neighbour cell
- Cells must be composed of planar faces
- "List" above, in the context of Julia, means a 1-D Array  

OpenFOAM meshes can be parsed into the format above using the "OpenFOAMMesh" function in "mesh.jl".
The current parse is slightly less flexible than OpenFOAM in terms of formatting, so if problems occur, try running the OpenFOAM utility "renumberMesh -overwrite" to ensure the mesh format is exactly as expected by mesh.jl.

# Solution State Definition
Mutable (modifiable) struct.  
Passed around FVM functions, intended to contain all universally-applicable info for FVM computations  
Examples offered are all for 1D. In a multidimensional computation, the states/residuals include an additional momentum term for each additional dimension, and the cell/face flux matrices will include one flux per coordinate direction per conserved variable.
For multidimensional computations, fluxes are ordered first by flux type, then by flux direction (ex: x-dirMassFlux, y-dirMassFlux, z-dirMassFlux, x-dirXMomentumFlux, y-dirXMomentumFlux, etc...)
```julia
solutionState =  
[  
   cellState,        # rho, rhoU, eV2 (conserved variables) at cell centers
   cellFluxes,       # massFluxes, momentumFluxes, energyFluxes at cell centers
   cellPrimitives,   # P, T, U at cell centers
   fluxResiduals,    # d/dt of each "cellState" variable (flux balances for each cell)
   faceFluxes        # massFluxes, momentumFluxes, and energyFluxes at face centers
]  
```
where:  
#### Cell State definition
```julia
# CellState =
# Cell      rho      x-Momentum   Total Energy
# Cell 1    rho_1    xMom_1         eV2_1
# Cell 2    rho_2    xMom_2         eV2_2
# ...
```
Ex. to access the x-Momentum in cell 23: `solutionState.cellState[23,2]`
#### Cell Fluxes definition
```julia
# CellFluxes =
# Cell      x_dirMassFlux   x-dir_x-MomentumFlux    x-dir_TotalEnergyFlux
# Cell 1    rhoU_1          (rho*U^2 + P)_1         (eV2*U + P*U)_1
# Cell 2    rhoU_2          (rho*U^2 + P)_2         (eV2*U + P*U)_2
# ...
```
#### Cell Primitives definition
```julia
# CellPrimitives =
# Cell      P       T       Ux
# Cell 1    P_1    T_1      Ux_1
# Cell 2    P_2    T_2      Ux_2
# ...
```
#### Flux Residuals definition
Note that "residual" is taken to mean flux balance, which is only strictly true in a steady-state computations.
In transient simulations, the "residual" represents the time derivative of the conserved quantity in each cell.
```julia
# fluxResiduals =
# Cell      rho            x-Momentum     Total Energy
# Cell 1    d(rho_1)/dt    d(xMom1)/dt      d(eV2_1)/dt
# Cell 2    d(rho_2)/dt    d(xMom2)/dt      d(eV2_2)/dt
# ...
```
#### Face Fluxes definition
```julia
# faceFluxes =
# Face      rho         x-Momentum            Total Energy
# Face 1    rhoU_f1     (rho*U^2 + P)_f1      (eV2*U + P*U)_f1
# Face 2    rhoU_f2     (rho*U^2 + P)_f2      (eV2*U + P*U)_f2
# ...
```