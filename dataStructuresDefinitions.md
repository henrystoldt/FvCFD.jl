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

# Unstructured FVM AdaptMesh Format Definition
Mutable struct
Created by passing in Mesh and mesh path to function  
Same as Mesh, but holds face point locations, and is mutable to allow addition of points and faces
Defined as follows:  
```julia
adaptMesh =
[
   cells,           #(list of lists, each containing the indices of the faces that make up the cell)  
   cVols,           #(list of scalar cell volumes)  
   cCenters,        #(list of position vectors of cell centers)  
   cellSizes,       #(Matrix: rows are cells, columns are x-, y-, and z-direction sizes of the cell)
   faces,           #(list of lists, each sublist containing two cell indices: the owner cell and the neighbour cell)  
   fAVecs,          #(list of face area vectors)  
   fCenters,        #(list of position vectors of face centers)  
   boundaryFaces,    #(list of lists, each containing the indices of faces on the ith boundary)  
   fPoints,         #(list of lists, each sublist containing the indices of points on the ith face Each matrix corresponds to a face)
   fPointLocs,      #(Matrix: rows are points, columns are x/y/z location)
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

# FacesData Format Defintion
Mutable struct
Created by passing info into createFacesDataStruct() function
Designed to be used for tracking changes to mesh during adaption process
```julia
faces =
[
   nFaces,        #(Integer, tracking the total number of faces in the mesh)
   nIntFaces,     #(Integer, tracking the number of internal faces in the mesh)
   nBdry,         #(Integer, tracking the number of boundary types)
   nBdryTypes,    #(list of type Any, holding the types of boundaries)
   nBdryIndices,  #(list of integers, holding the indices at which each boundary type begins in the faces list)
]
```

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
   faceState         # rho, rhoU, eV2 (conserved variables) at face centers. Used primarily for holding linInterp results
   faceFluxes        # massFluxes, momentumFluxes, and energyFluxes at face centers
   facePrimitives    # P, T, U at face centers
]  
```
Where  
#### Cell State definition
```julia
# CellState =
# Cell      rho      x-Momentum   Total Energy
# Cell 1    rho_1    xM_1         eV2_1
# Cell 2    rho_2    xM_2         eV2_2
# ...
```
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
Note that "residual" is taken to mean flux balance, which is only strictly true in a steady-state computations. Term used loosely here.
```julia
# fluxResiduals =
# Cell      rho            x-Momentum     Total Energy
# Cell 1    d(rho_1)/dt    d(xM_1)dt      d(eV2_1)/dt
# Cell 2    d(rho_2)/dt    d(xM_2)dt      d(eV2_2)/dt
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
