using Plots
using DelimitedFiles
using WriteVTK
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

function outputVTK(meshPath, cellPrimitives, fileName="solution")
    points, cellPtIndices = OpenFOAMMesh_findCellPts(meshPath)
    points = transpose(points)
    cells = Array{MeshCell, 1}(undef, 0)

    cellType = [ 1, 3, 5, 10, 14, 13, "ERROR", 12 ] # This array maps from number of points in a cell to the .vtk numeric cell type. Example: 8 pts -> "12", which is .vtk code for "VTK_HEXAHEDRON"
    # Corresponding .vtk cell types: [ "VTK_VERTEX", "VTK_LINE", "VTK_TRIANGLE", "VTK_TETRA", "VTK_PYRAMID", "VTK_WEDGE", "ERROR", "VTK_HEXAHEDRON" ]

    for i in eachindex(cellPtIndices)
        nPoints = size(cellPtIndices[i], 1)
        cT = cellType[nPoints]
        cell = MeshCell(VTKCellType(cT), cellPtIndices[i])
        push!(cells, cell)
    end

    file = vtk_grid(fileName, points, cells)
    file["P"] = cellPrimitives[:,1]
    file["T"] = cellPrimitives[:,2]
    file["U"] = transpose(cellPrimitives[:,3:5])
    
    return vtk_save(file)
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
        solnName = "solution.$vtkCounter"
        println("Writing $solnName")
        outputVTK(meshPath, cellPrimitives, solnName)
    end
end
