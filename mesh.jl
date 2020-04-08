# Methods from Moukalled et al. FVM - OpenFOAM, Matlab
include("vectorFunctions.jl")
include("dataStructures.jl")

######################### Mesh/Cell Geometry ###########################
function crossProd(v1::Array{Float64, 1}, v2::Array{Float64, 1})
    x = v1[2]*v2[3] - v1[3]*v2[2]
    y = -(v1[1]*v2[3] - v1[3]*v2[1])
    z = v1[1]*v2[2] - v1[2]*v2[1]
    return [x,y,z]
end

function triangleCentroid(points)
    center = [ 0.0, 0.0, 0.0 ]
    nPts = size(points, 1)
    for pt in 1:nPts
        center = center .+ points[pt]
    end
    center /= nPts

    return center
end

#Alternate name for same function
geometricCenter = triangleCentroid

function triangleArea(points::Array{Array{Float64, 1}})
    side1 = points[2] .- points[1]
    side2 = points[3] .- points[1]
    fAVec = crossProd(side1, side2) ./ 2
    return fAVec
end

"""
    Calculates face area vector and centroid from the points that make up the face
    Points must be ordered sequentially

    How it works:
        1. Splits face into subtriangles
        2. Area and centroid is computed for each subtriangle
        3. Areas vectors are summed and returned
        4. The centroid returned is obtained from an area-weighted of sum of the subtriangle centroids
"""
function faceAreaCentroid(points::Array{Array{Float64, 1}})
    gC = geometricCenter(points)
    nPts = size(points, 1)

    fAVec = [ 0.0, 0.0, 0.0 ]
    centroid = [ 0.0, 0.0, 0.0 ]

    for i in 1:nPts
        if i < nPts
            subTriPts = [ gC, points[i], points[i+1] ]
        else
            subTriPts = [ gC, points[i], points[1] ]
        end

        triCentroid = triangleCentroid(subTriPts)
        subFAVec = triangleArea(subTriPts)

        fAVec += subFAVec
        centroid += triCentroid .* mag(subFAVec)
    end

    centroid /= mag(fAVec)

    return fAVec, centroid
end

"""
    Calculates cell volume (scalar) and centroid (vector) from the points and face area vectors (fAVecs) that make up the cell
        fAVecs can be computed using the faceAreaCentroids function

    How it works:
        1. Splits cell into polygonal pyramids, each incorporating a single face and the geometric center of the cell
        2. Computes volume and centroid of each sub-pyramid
        3. Resulting volume is sum of sub-pyramid volumes, centroid is the volume-weighted sum of sub-pyramid centroids
"""
function cellVolCentroid(points, fAVecs, faceCentroids)
    gC = geometricCenter(points)
    nFaces = size(fAVecs,1)

    vol = 0.0
    centroid = [ 0.0, 0.0, 0.0 ]

    for f in 1:nFaces
        cellCenterVec = faceCentroids[f] .- gC
        subPyrVol = abs(sum(fAVecs[f] .* cellCenterVec) / 3)
        subPyrCentroid = 0.75.*faceCentroids[f] .+ 0.25.*gC

        vol += subPyrVol
        centroid += subPyrCentroid .* subPyrVol
    end

    centroid /= vol

    return vol, centroid
end

# For each face, calculates the vectors from it's owner and neighbour cell centers to it's own center
function cellCentroidToFaceVec(faceCentroids::Array{Array{Float64, 1}}, cellCentroids::Array{Array{Float64, 1}})
    #TODO: Understand ordering to make sure we're accepting faces in the order we need to pass them back in??
    nFaces = size(faceCentroids, 1)

    cellToFaceVec = zeros(nFaces, 3)

    for f in 1:nFaces
        cellToFaceVec[f,:] = faceCentroids[f] - cellCentroids
    end

    return cellToFaceVec
end

######################### Utility Functions ###########################
# Returns basic info about the mesh: nCells, nFaces, nBoundaries, nBdryFaces
function unstructuredMeshInfo(mesh::Mesh)
    nCells = size(mesh.cells, 1)
    nFaces = size(mesh.faces, 1)
    nBoundaries = size(mesh.boundaryFaces, 1)

    # Count boundary faces
    nBdryFaces = 0
    for bdry in 1:nBoundaries
        nBdryFaces += size(mesh.boundaryFaces[bdry], 1)
    end

    return nCells, nFaces, nBoundaries, nBdryFaces
end

# Checks whether a string represents a number (using a regular expression)
function isNumber(str)
    re = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$"
    return occursin(re, str)
end

######################### Parse OpenFOAM Meshes ###########################
"""
    In the points, faces, owner, and neighbour files, the beginning of useful information in the file looks like this:

        ...
        object      faces;
        }
        // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


        5135
        (
        4(79 1516 1380 37)
        4(1381 1422 1249 39)
        ...

    The 5135 indicates how many lines of information follow it.
    This function would return 5135, and the line number of the first piece of information in the file ((line number of 5135) + 2)
"""
function OFFile_FindNItems(fileLines)
    nLines = size(fileLines, 1)
    lineCounter = 1
    itemCount = 0
    startLine = 0
    for line in fileLines
        if isNumber(strip(line))
            itemCount = parse(Int64, line)
            startLine = lineCounter+2
            break
        end
        lineCounter += 1
    end
    return startLine, itemCount
end

"""
    Takes the file path to an OpenFOAM mesh 'points' file, returns a 2D array of point coordinates:
            X   Y   Z
        P1  x1  y1  z1
        P2  x2  y2  z2
        ...
"""
function readOFPointsFile(filePath)
    f = open(filePath)
    pointsLines = readlines(f)
    startLine, pCount = OFFile_FindNItems(pointsLines)

    points = zeros(pCount, 3)
    for line in startLine:(startLine+pCount-1)
        pLine = pointsLines[line]

        bracketsRemoved = pLine[2:end-1]
        coords = split(bracketsRemoved)
        for d in 1:3
            points[line+1-startLine, d] = parse(Float64, coords[d])
        end
    end

    return points
end

"""
    Takes the file path to an OpenFOAM mesh 'faces' file, returns a 1D array of 1D arrays.
    Each entry in the upper level array represents a single face
    Each face-array contains the indices of the points that make up that face
    Note that all indices from the file are incremented by 1, to switch from OpenFOAM's 0-based indexing to Julia's 1-based indexing

    Example return value:
        [ [face1_pt1Index, face1_pt2Index, face1_pt3Index], [face2_pt1Index, face2_pt2Index, face2_pt3Index, face2_pt4Index], ... ]
"""
function readOFFacesFile(filePath)
    f = open(filePath)
    facesLines = readlines(f)
    startLine, fCount = OFFile_FindNItems(facesLines)

    faces = []
    for line in startLine:(startLine+fCount-1)
        fLine = facesLines[line]
        facePts = []

        bracketIndex = findfirst("(", fLine)[1]
        bracketsRemoved = fLine[bracketIndex+1:end-1]
        ptNumbers = split(bracketsRemoved)
        nPts = size(ptNumbers, 1)
        for p in 1:nPts
            # Add one to convert to 1-based indexing
            push!(facePts, parse(Int64, ptNumbers[p])+1)
        end
        push!(faces, facePts)
    end

    return faces
end

"""
    Takes the file path to an OpenFOAM 'owner' file, returns a 1D array of owner cell indices
    Example return value:
        [ face1_ownerCellIndex, face2_ownerCellIndex, ...  ]
"""
function readOFOwnerFile(filePath)
    f = open(filePath)
    ownerLines = readlines(f)
    startLine, oCount = OFFile_FindNItems(ownerLines)

    faceOwnerCells = zeros(Int64, oCount)
    for line in startLine:(startLine+oCount-1)
        # Add one to convert to 1-based indexing
        faceOwnerCells[line-startLine+1] = parse(Int64, ownerLines[line])+1
    end

    return faceOwnerCells
end

"""
    Takes the file path to an OpenFOAM 'neighbour' file, returns a 1D array of neighbour cell indices
    Note that all indices from the file are incremented by 1, to switch from OpenFOAM's 0-based indexing to Julia's 1-based indexing

    Example return value:
        [ face1_neighbourCellIndex, face2_neighbourCellIndex, ...  ]
"""
function readOFNeighbourFile(filePath)
    f = open(filePath)
    neighbourLines = readlines(f)
    startLine, nCount = OFFile_FindNItems(neighbourLines)

    faceNeighbourCells = zeros(Int64, nCount)
    for line in startLine:(startLine+nCount-1)
        # Add one to convert to 1-based indexing
        faceNeighbourCells[line-startLine+1] = parse(Int64, neighbourLines[line])+1
    end

    return faceNeighbourCells
end

"""
    Returns index of first line containing str
    Otherwise returns -1
"""
function findInLines(str, lines, startLine)
    nLines = size(lines, 1)
    for i in startLine:nLines
        if occursin(str, lines[i])
            return i
        end
    end
    return -1
end

"""
    Takes the file path to an OpenFOAM 'boundary' file, returns 3 arrays:
    boundaryNames:          Array of boundary names (strings)
    boundaryNumFaces:       Array containing the number of faces in each boundary
    boundaryStartFaces:     Array containing the index of the first face in the boundary (boundaries always occupy a continuous sequence of face indices, and are usually numbered last)
"""
function readOFBoundaryFile(filePath)
    f = open(filePath)
    bLines = readlines(f)
    startLine, bCount = OFFile_FindNItems(bLines)
    nLines = size(bLines, 1)

    boundaryNames = Array{String, 1}(undef, bCount)
    boundaryStartFaces = Array{Int64, 1}(undef, bCount)
    boundaryNumFaces = Array{Int64, 1}(undef, bCount)
    for i in 1:bCount
        bNameLine = findInLines("{", bLines, startLine)-1
        boundaryNames[i] = strip(bLines[bNameLine])

        bNFacesLine = findInLines("nFaces", bLines, startLine)
        boundaryNumFaces[i] = parse(Int64, split(bLines[bNFacesLine])[2][1:end-1])

        bStartFaceLine = findInLines("startFace", bLines, startLine)
        boundaryStartFaces[i] = parse(Int64, split(bLines[bStartFaceLine])[2][1:end-1])+1

        startLine = findInLines("}", bLines, startLine)+1
    end

    return boundaryNames, boundaryNumFaces, boundaryStartFaces
end

"""
    Input: Path to an OpenFOAM mesh FOLDER.
    Output: Calls readOFPoints/Faces/Owner/NeighbourFile and returns all of their results (reads all OpenFOAM's mesh file data into arrays)
"""
function readOpenFOAMMesh(polyMeshPath)
    pointsFilePath = "$polyMeshPath/points"
    points = readOFPointsFile(pointsFilePath)
    facesFilePath = "$polyMeshPath/faces"
    faces = readOFFacesFile(facesFilePath)
    ownerFilePath = "$polyMeshPath/owner"
    owner = readOFOwnerFile(ownerFilePath)
    neighbourFilePath = "$polyMeshPath/neighbour"
    neighbour = readOFNeighbourFile(neighbourFilePath)
    boundaryFilePath = "$polyMeshPath/boundary"
    boundaryNames, boundaryNumFaces, boundaryStartFaces = readOFBoundaryFile(boundaryFilePath)

    return points, faces, owner, neighbour, boundaryNames, boundaryNumFaces, boundaryStartFaces
end

"""
    Function used to find all the points in each cell.
    Not required for CFD, which is face-based, but required for cell-based .vtk file output.

    Input: Path to an OpenFOAM mesh FOLDER.
    Returns:
        points:         Array of points obtained from readOFPoints()
        cellPtIndices:  Array of arrays, where each subarrary represents a cell, and each entry in a cell's array is the index of one of its points
"""
function OpenFOAMMesh_findCellPts(polyMeshPath)
    points, OFfaces, owner, neighbour, boundaryNames, boundaryNumFaces, boundaryStartFaces = readOpenFOAMMesh(polyMeshPath)
    nCells = maximum(owner)
    nFaces = size(OFfaces, 1)
    nBoundaries = size(boundaryNames, 1)

    cellPtIndices = Array{Array{Int64, 1}, 1}(undef, nCells)    # OUTPUT, Array of arrays, where each subarrary represents a cell, and each entry in a cell's array is the index of one of its points
    cells = Array{Array{Int64, 1}, 1}(undef, nCells)            # Indices of faces that make up a cell

    ### Populate the cells array ###
    #TODO: This code to make the cells array is nearly duplicated in OpenFOAMMesh() below. Perhaps extract it into a function.
    nOwners = size(owner, 1)
    for f in 1:nOwners
        ownerCell = owner[f]
        if isassigned(cells, ownerCell)
            push!(cells[ownerCell], f)
        else
            cells[ownerCell] = [f,]
        end
    end

    nNeighbours = size(neighbour, 1)
    for f in 1:nNeighbours
        neighbourCell = neighbour[f]
        if isassigned(cells, neighbourCell)
            push!(cells[neighbourCell], f)
        else
            cells[neighbourCell] = [f,]
        end
    end
    # fAVecs, fCenters, faces, cells now complete

    ### Using the cells array, populate the cellPtIndices array ###
    for c in 1:nCells
        pts = Array{Int64,1}(undef, 0) # Points that make up the present cell

        # Add all points from face f to pts
        function addFace(f)
            for pt in OFfaces[f]
                push!(pts, pt)
            end
        end

        # Are all points in face f currently absent from the pts array (are the arrays OFfaces[f] and pts disjont?)
        function disjoint(f)
            for pt in OFfaces[f]
                if any(x->x==pt, pts)
                    return false
                end
            end
            return true
        end

        # First face to add to the pts array. Changing the index of this face (which must be <= Number of faces in smallest cell) can change in which order points are added to the array.
            # This sometimes happens to put points in the correct order for nice-looking .vtk output
            # Current best options for each mesh: Wedge - 3
        addFace(cells[c][1])

        # Add any disjoint faces
        for f in cells[c]
            if disjoint(f)
                addFace(f)
            end
        end

        # Add any leftover points
        for f in cells[c]
            for pt in OFfaces[f]
                if !any(x->x==pt, pts)
                    push!(pts, pt)
                end
            end
        end

        # Store the list of points that make up the current cell
        cellPtIndices[c] = pts
    end

    return points, cellPtIndices
end

"""
    Reads an OpenFOAM mesh and returns a Mesh object suitable for calculations in JuliaCFD.

    Mesh defined in dataStructures.jl, documented in dataStructureDefintions.md
"""
function OpenFOAMMesh(polyMeshPath)
    points, OFfaces, owner, neighbour, boundaryNames, boundaryNumFaces, boundaryStartFaces = readOpenFOAMMesh(polyMeshPath)
    nCells = maximum(owner)
    nFaces = size(OFfaces, 1)
    nBoundaries = size(boundaryNames, 1)

    cells = Array{Array{Int64, 1}, 1}(undef, nCells)                # Indices of faces that make up a cell
    cVols = zeros(nCells)                                           # Cell Volumes
    cCenters = Array{Array{Float64, 1}, 1}(undef, nCells)           # Cell Centroids
    cellSizes = Matrix{Float64}(undef, nCells, 3)                   # Cell sizes (x,y,z) direction

    faces = Array{Array{Int64, 1}, 1}(undef, nFaces)                # Indices of cells adjacent to each face always owner first, then neighbour
    fAVecs = Array{Array{Float64, 1}, 1}(undef, nFaces)             # Face Area Vectors
    fCenters = Array{Array{Float64, 1}, 1}(undef, nFaces)           # Face Centroids
    boundaryFaces = Array{Array{Int64, 1}, 1}(undef, nBoundaries)   # Indices of faces that make up each boundary

    # Calculate face area vectors and face centroids from list of points in each face
    for f in 1:nFaces
        fPts = [ points[pt,:] for pt in OFfaces[f] ]
        fAVecs[f], fCenters[f] = faceAreaCentroid(fPts)
    end
    #fAVecs and fCenters now complete

    # Build 'cells' and 'faces' arrays from info from 'owner' file
    nOwners = size(owner, 1)
    for f in 1:nOwners
        ownerCell = owner[f]
        if isassigned(cells, ownerCell)
            push!(cells[ownerCell], f)
        else
            cells[ownerCell] = [f,]
        end

        faces[f] = [ownerCell, -1]
    end

    # Finish building 'cells' and 'faces' arrays, adding info from 'neighbour' file
    nNeighbours = size(neighbour, 1)
    for f in 1:nNeighbours
        neighbourCell = neighbour[f]
        if isassigned(cells, neighbourCell)
            push!(cells[neighbourCell], f)
        else
            cells[neighbourCell] = [f,]
        end

        faces[f][2] = neighbourCell
    end
    # fAVecs, fCenters, faces, cells now complete

    # For each cell, use the face area vectors and face centroids for all of its faces, to calculate it's volume and centroid location
    for c in 1:nCells
        pts = []
        for f in cells[c]
            for pt in OFfaces[f]
                if !(points[pt,:] in pts)
                    push!(pts, points[pt,:])
                end
            end
        end
        fCs = [ fCenters[f] for f in cells[c] ]
        cell_fAVecs = [ fAVecs[f] for f in cells[c] ]

        cVols[c], cCenters[c] = cellVolCentroid(pts, cell_fAVecs, fCs)
    end
    # fAVecs, fCenters, faces, cells, cVols, cCenters now complete

    # Create boundaryFaces array
    for b in 1:nBoundaries
        startF = boundaryStartFaces[b]
        endF = startF + boundaryNumFaces[b] - 1
        boundaryFaces[b] = Array(startF:endF)
    end

    # Compute cell sizes in the x, y, z directions
    # In the past, cell sizes were used for a shitty CFL calculation, may be able to get rid of them now
    maxCoords = [ -1000000.0, -1000000.0, -1000000.0 ]
    minCoords = [ 1000000.0, 1000000.0, 1000000.0 ]
    for c in 1:nCells
        fill!(maxCoords, -10000000)
        fill!(minCoords, 10000000)

        pts = Vector{Vector{Float64}}()

        # Add points
        for f in cells[c]
            for pt in OFfaces[f]
                if !any(x->x==pt, pts)
                    push!(pts, points[pt, :])
                end
            end
        end

        for pt in pts
            for d in 1:3
                maxCoords[d] = max(maxCoords[d], pt[d])
                minCoords[d] = min(minCoords[d], pt[d])
            end
        end

        for d in 1:3
            cellSizes[c,d] = maxCoords[d] - minCoords[d]
        end
    end

    return Mesh(cells, cVols, cCenters, cellSizes, faces, fAVecs, fCenters, boundaryFaces)
end
