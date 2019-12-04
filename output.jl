using Plots
using DelimitedFiles
include("mesh.jl")
pyplot()

function plot2DResult(mesh, cellValue)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells, nFaces, nBoundaries, nBdryFaces = unstructuredMeshInfo(mesh)
    x = zeros(nCells)
    y = zeros(nCells)
    for c in 1:nCells
        x[c] = cCenters[c][1]
        y[c] = cCenters[c][2]
    end

    plot(x, y, P, st=:surface)
    gui()
end

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
        write(f, "DATASET POLYDATA\n")

        #### Output POINTS ####
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
        #cellType = [ "VTK_VERTEX", "VTK_LINE", "VTK_TRIANGLE", "VTK_TETRA", "VTK_PYRAMID", "VTK_WEDGE", "ERROR", "VTK_HEXAHEDRON" ]
        cellType = [ "1", "3", "5", "10", "14", "13", "ERROR", "12" ]
        totalCellListSize = nCells
        for c in 1:nCells
            totalCellListSize += size(cellPtIndices[c], 1)
        end

        write(f, "\nPOLYGONS $nCells $totalCellListSize\n")
        for c in 1:nCells
            ptCount = size(cellPtIndices[c], 1)
            str = "$ptCount"
            for pt in cellPtIndices[c]
                str = string(str, " ", pt-1 )
            end
            str = string(str, "\n")
            write(f, str)
        end

        # write(f, "CELL_TYPES $nCells\n")
        # for c in 1:nCells
        #     nPts = size(cellPtIndices[c], 1)
        #     cT = cellType[nPts]
        #     write(f, "$cT\n")
        # end

        #### Output SCALARS, VECTORS ####
        write(f, "\nCELL_DATA $nCells\n")
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
        write(f, "VECTORS U float\n")
        for c in 1:nCells
            Ux, Uy, Uz = cellPrimitives[c, 3:5]
            write(f, "$Ux $Uy $Uz\n")
        end
    end
end

function updateSolutionOutput(cellPrimitives, restartFile, meshPath, vtkCounter, createRestartFile)
    if createRestartFile
        println("Writing Restart File: $restartFile")
        writeRestartFile(cellPrimitives, restartFile)
    end
    solnName = "solution.$vtkCounter.vtk"
    println("Writing $solnName")
    outputVTK(meshPath, cellPrimitives, solnName)

end
