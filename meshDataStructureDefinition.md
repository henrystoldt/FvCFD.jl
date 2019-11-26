# Mesh Format Definition
Mesh data structure is defined as follows:  
```julia
mesh =
[
   cells,           #(list of lists, each containing the indices of the faces that make up the cell)  
   faces,           #(list of lists, each sublist containing two cell indices: the owner cell and the neighbour cell)  
   fAVectors,       #(list of face area vectors)  
   fCenters,        #(list of position vectors)  
   boundaryFaces,   #(list of lists, each containing the indices of cells on the ith boundary)  
   cellVolumes,     #(list of scalar cell volumes)  
   cellCenters      #(list of position vectors)  
]
```

## Notes
- Cells and faces are numbered according to their storage location in the mesh arrays
- Faces are numbered such that boundary faces come last
- Face area vector point outward from the owner cells, into the neighbour cell
- Cells must be convex, composed of planar faces
- "List" above, in the context of Julia, means a 1-D Array