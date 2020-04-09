using Plots
using DelimitedFiles
include("mesh.jl")
pyplot()

#=
    Displays surface plot of cellValue for a 2D mesh.
=#
function plot2DResult(mesh, cellValue)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    x = zeros(nCells)
    y = zeros(nCells)
    for c in 1:nCells
        x[c] = cCenters[c][1]
        y[c] = cCenters[c][2]
    end

    plot(x, y, callValue, st=:surface)
    gui()
end

#### Read / Write Restart Files ####
function writeRestartFile(cellPrimitives, path="JuliaCFDRestart.txt")
    writedlm(path, cellPrimitives)
end

function readRestartFile(path="JuliaCFDRestart.txt")
    cellPrimitives = readdlm(path)
end


function outputVTK(meshPath, cellPrimitives, fileName="solution.vtk")
    open(fileName, "w") do f
        write(f, "# vtk DataFile Version 2.0\n")
        write(f, "JuliaCFD\n")
        write(f, "ASCII\n")
        write(f, "DATASET UNSTRUCTURED_GRID\n")

        #### Output all POINTS in the mesh, one per line ####
        points, cellPtIndices = OpenFOAMMesh_findCellPts(meshPath)
        nCells = size(cellPtIndices, 1)
        nPts = size(points, 1)
        write(f, "\nPOINTS $nPts float\n")
        for pt in 1:nPts
            x = points[pt, 1]
            y = points[pt, 2]
            z = points[pt, 3]
            write(f, "$x $y $z\n")
        end

        #### Output CELLS, and CELL_TYPES ####
        # Count and output the total number of points that we will have to output to make up all cells
        totalCellListSize = nCells
        for c in 1:nCells
            totalCellListSize += size(cellPtIndices[c], 1)
        end
        write(f, "\nCELLS $nCells $totalCellListSize\n")

        # Output one cell per line. The cell is represented by the indices of all the points that make it up (referring to the list of points printed earlier). Point indices are decremented by 1 to convert from Julia's 1-based indexing to .vtk's 0-based indexing
        for c in 1:nCells
            ptCount = size(cellPtIndices[c], 1)
            str = "$ptCount"
            for pt in cellPtIndices[c]
                str = string(str, " ", pt-1 )
            end
            str = string(str, "\n")
            write(f, str)
        end

        # Output the .vtk "type" of each cell, one per line
        cellType = [ "1", "3", "5", "10", "14", "13", "ERROR", "12" ] # This array maps from number of points in a cell to the .vtk numeric cell type. Example: 8 pts -> "12", which is .vtk code for "VTK_HEXAHEDRON"
            # Corresponding .vtk cell types: [ "VTK_VERTEX", "VTK_LINE", "VTK_TRIANGLE", "VTK_TETRA", "VTK_PYRAMID", "VTK_WEDGE", "ERROR", "VTK_HEXAHEDRON" ]
        write(f, "CELL_TYPES $nCells\n")
        for c in 1:nCells
            nPts = size(cellPtIndices[c], 1)
            cT = cellType[nPts]
            write(f, "$cT\n")
        end

        #### Output SCALAR, VECTOR data at each cell center ####
        write(f, "\nCELL_DATA $nCells\n")

        # Pressure, Temperature
        dataNames = [ "P", "T" ]
        for d in 1:2
            dataName = dataNames[d]
            write(f, "SCALARS $dataName float 1\n") # Name=P dataType=float, 1 component
            write(f, "LOOKUP_TABLE default\n") # No custom color lookup table
            for c in 1:nCells
                P = cellPrimitives[c, d]
                write(f, "$P\n")
            end
        end

        # Velocity
        write(f, "VECTORS U float\n")
        for c in 1:nCells
            Ux, Uy, Uz = cellPrimitives[c, 3:5]
            write(f, "$Ux $Uy $Uz\n")
        end
    end
end

#=
    Calls above functions to output restart and .vtk files, if desired.

    Inputs:
        cellPrimitives: should come from sln.cellPrimitives
        restartFile: (string) path to which to write restart file
        meshPath: (string) path to OpenFOAM mesh FOLDER
        createRestartFile: (bool)
        createVTKOutput: (bool)

    Will overwrite existing restart files
    Will not overwrite existing .vtk files
=#
function updateSolutionOutput(cellPrimitives, restartFile, meshPath, createRestartFile, createVTKOutput)
    if createRestartFile
        println("Writing Restart File: $restartFile")
        writeRestartFile(cellPrimitives, restartFile)
    end

    if createVTKOutput
        # Check for next available filename
        files = readdir()
        maxNum = 0
        for item in files
            if occursin("solution", item)
                slnNumber = parse(Int, item[10:end-4])
                maxNum = max(maxNum, slnNumber)
            end
        end
        vtkCounter = maxNum + 1

        # Write vtk file
        solnName = "solution.$vtkCounter.vtk"
        println("Writing $solnName")
        outputVTK(meshPath, cellPrimitives, solnName)
    end
end
