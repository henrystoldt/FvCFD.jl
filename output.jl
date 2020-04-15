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
        # display(points)
        # display(cellPtIndices)
        # println("$breakdown")
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

function writeNewOpenFOAMMesh(mesh::AdaptMesh, facesData; meshNewPath="./adapt/")
    home = pwd()
    #Find number for mesh
    cd(meshNewPath)

    maxNum = 0

    files = readdir()
    for item in files
        if occursin("mesh_", item)
            slnNumber = parse(Int, item[6:end])
            maxNum = max(maxNum, slnNumber)
        end
    end
    meshCount = maxNum + 1

    meshName = "mesh_$meshCount"
    newMeshPath = "./adapt/mesh_$meshCount/"
    println("Writing new mesh to $newMeshPath")
    mkdir(meshName)
    cd(meshName)

    cells = mesh.cells
    pointLocs = mesh.fPointLocs
    fPoints = mesh.fPoints
    faces = mesh.faces

    nPoints = size(pointLocs, 1)
    nCells = size(cells,1)
    nFaces = size(fPoints,1)
    nInternalFaces = facesData.nIntFaces

    note = "    note        \"nPoints:$nPoints  nCells:$nCells  nFaces:$nFaces  nInternalFaces:$nInternalFaces\";"


    writePointsFile(pointLocs)
    writeFacesFile(fPoints)
    writeOwnerFile(faces, note=note)
    writeNeighbourFile(faces, nInternalFaces, note=note)
    writeBoundaryFile(facesData)

    cd(home)

    return newMeshPath
end

function writePointsFile(pointLocs; fileName="points")
    open(fileName, "w") do f
        writeOpenFOAMHeader(f)
        write(f, "\r\n")

        nPoints = size(pointLocs, 1)
        write(f, "$nPoints\r\n")
        write(f, "(\r\n")

        for i in 1:nPoints
            x = pointLocs[i,1]
            y = pointLocs[i,2]
            z = pointLocs[i,3]
            write(f, "($x $y $z)\r\n")
        end

        write(f, ")\r\n")
        write(f, "\r\n")
        write(f, "\r\n")
        write(f, "// ************************************************************************* //")
    end

    return
end

function writeFacesFile(fPoints; fileName="faces")
    open(fileName, "w") do f
        writeOpenFOAMHeader(f, class="faceList", object="faces")
        write(f, "\r\n")

        nFaces = size(fPoints, 1)
        write(f, "$nFaces\r\n")
        write(f, "(\r\n")

        for i in 1:nFaces
            sizeOfFace = size(fPoints[i],1)
            if sizeOfFace == 3
                f1 = fPoints[i][1] - 1
                f2 = fPoints[i][2] - 1
                f3 = fPoints[i][3] - 1
                write(f, "$sizeOfFace($f1 $f2 $f3)\r\n")
            elseif sizeOfFace == 4
                f1 = fPoints[i][1] - 1
                f2 = fPoints[i][2] - 1
                f3 = fPoints[i][3] - 1
                f4 = fPoints[i][4] - 1
                write(f, "$sizeOfFace($f1 $f2 $f3 $f4)\r\n")
            elseif sizeOfFace == 5
                f1 = fPoints[i][1] - 1
                f2 = fPoints[i][2] - 1
                f3 = fPoints[i][3] - 1
                f4 = fPoints[i][4] - 1
                f5 = fPoints[i][5] - 1
                write(f, "$sizeOfFace($f1 $f2 $f3 $f4 $f5)\r\n")
            else
                println("Unexpected size of face while writing mesh files!")
                println("$breakdown")
            end
        end

        write(f, ")\r\n")
        write(f, "\r\n")
        write(f, "\r\n")
        write(f, "// ************************************************************************* //")
    end

    return
end

function writeOwnerFile(faces; note="", fileName="owner")
    open(fileName, "w") do f

        class = "labelList"
        object = "owner"
        writeOpenFOAMHeader(f, class="labelList", note=note, object="owner")
        write(f, "\r\n")

        nFaces = size(faces, 1)
        write(f, "$nFaces\r\n")
        write(f, "(\r\n")

        for i in 1:nFaces
            owner = faces[i][1] - 1
            write(f, "$owner\r\n")
        end

        write(f, ")\r\n")
        write(f, "\r\n")
        write(f, "\r\n")
        write(f, "// ************************************************************************* //")
    end

    return
end

function writeNeighbourFile(faces, nInternalFaces; note="", fileName="neighbour")
    open(fileName, "w") do f
        class = "labelList"
        object = "neighbour"
        writeOpenFOAMHeader(f, class="labelList", note=note, object="neighbour")
        write(f, "\r\n")

        nFaces = size(faces, 1)
        write(f, "$nInternalFaces\r\n")
        write(f, "(\r\n")

        for i in 1:nInternalFaces
            neighbour = faces[i][2] - 1
            if neighbour == -1
                println("Cell values have been applied incorrectly! Found a boundary cell in the neighbours file!")
                # println("Internal faces: $nInternalFaces")
                # println("Face at break $i")
                # display("Counter at break $counter")
                println("$breakdown")
            end
            write(f, "$neighbour\r\n")
        end

        write(f, ")\r\n")
        write(f, "\r\n")
        write(f, "\r\n")
        write(f, "// ************************************************************************* //")
    end

    return
end

function writeBoundaryFile(facesData; fileName="boundary")
    open(fileName, "w") do f
        class = "polyBoundaryMesh"
        object = "boundary"
        writeOpenFOAMHeader(f, class="polyBoundaryMesh", object="boundary")

        nBdry = facesData.nBdry
        write(f, "$nBdry\r\n")
        write(f, "(\r\n")
        nFaces = 0
        zeroGradCounter = 0

        for b in 1:nBdry
            nStart = facesData.bdryIndices[b] - 1
            if b == nBdry
                nFaces = facesData.nFaces - nStart
            else
                nFaces = facesData.bdryIndices[b+1] - nStart - 1
            end

            if facesData.bdryTypes[b] == emptyBoundary!
                write(f, "    defaultFaces\r\n")
                write(f, "    {\r\n")
                write(f, "        type            empty;\r\n")
                write(f, "        inGroups        1(empty);\r\n")
                write(f, "        nFaces          $nFaces;\r\n")
                write(f, "        startFace       $nStart;\r\n")
                write(f, "    }\r\n")
            elseif facesData.bdryTypes[b] == supersonicInletBoundary!
                write(f, "    inlet\r\n")
                write(f, "    {\r\n")
                write(f, "        type            patch;\r\n")
                write(f, "        nFaces          $nFaces;\r\n")
                write(f, "        startFace       $nStart;\r\n")
                write(f, "    }\r\n")
            elseif facesData.bdryTypes[b] == zeroGradientBoundary! && zeroGradCounter == 0
                write(f, "    outlet\r\n")
                write(f, "    {\r\n")
                write(f, "        type            patch;\r\n")
                write(f, "        nFaces          $nFaces;\r\n")
                write(f, "        startFace       $nStart;\r\n")
                write(f, "    }\r\n")
                zeroGradCounter += 1
            elseif facesData.bdryTypes[b] == symmetryBoundary!
                write(f, "    bottom\r\n")
                write(f, "    {\r\n")
                write(f, "        type            symmetryPlane;\r\n")
                write(f, "        inGroups        1(symmetryPlane);\r\n")
                write(f, "        nFaces          $nFaces;\r\n")
                write(f, "        startFace       $nStart;\r\n")
                write(f, "    }\r\n")
            elseif facesData.bdryTypes[b] == wallBoundary!
                write(f, "    obstacle\r\n")
                write(f, "    {\r\n")
                write(f, "        type            patch;\r\n")
                write(f, "        nFaces          $nFaces;\r\n")
                write(f, "        startFace       $nStart;\r\n")
                write(f, "    }\r\n")

            elseif facesData.bdryTypes[b] == zeroGradientBoundary! && zeroGradCounter == 1
                write(f, "    top\r\n")
                write(f, "    {\r\n")
                write(f, "        type            patch;\r\n")
                write(f, "        nFaces          $nFaces;\r\n")
                write(f, "        startFace       $nStart;\r\n")
                write(f, "    }\r\n")
                zeroGradCounter += 1
            #elseif facesData.bdryTypes[b] == supersonicInletBoundary!
            else
                println("Boundary type not recognized when writing boundary file! Update list of boundary types!")
                println("$breakdown")
            end
        end

        write(f, ")\r\n")
        write(f, "\r\n")
        write(f, "\r\n")
        write(f, "// ************************************************************************* //")
    end

    return
end

function writeOpenFOAMHeader(f; class="vectorField", note="", object="points")
    write(f, "/*--------------------------------*- C++ -*----------------------------------*\\\n")
    write(f, "  =========                 |\n")
    write(f, "  \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox\n")
    write(f, "   \\\\    /   O peration     | Website:  https://openfoam.org\n")
    write(f, "    \\\\  /    A nd           | Version:  6\n")
    write(f, "     \\\\/     M anipulation  |\n")
    write(f, "\\*---------------------------------------------------------------------------*/\n")
    write(f, "FoamFile\n")
    write(f, "{\n")
    write(f, "    version     2.0;\n")
    write(f, "    format      ascii;\n")
    write(f, "    class       $class;\n")
    if note != ""
        write(f, "$note\n")
    end
    write(f, "    location    \"constant/polyMesh\";\r\n")
    write(f, "    object      $object;\r\n")
    write(f, "}\r\n")
    write(f, "// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\r\n")
    write(f, "\r\n")

end
