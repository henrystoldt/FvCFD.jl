# Methods from Moukalled et al. FVM - OpenFOAM, Matlab
include("vectorFunctions.jl")
include("dataStructures.jl")

######################### Mesh/Cell Geometry ###########################
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
    fAVec = cross(side1, side2) ./ 2
    return fAVec
end

#=
    Calculates face area vector and centroid from the points that make up the face
    Points must be ordered sequentially

    How it works:
        1. Splits face into subtriangles
        2. Area and centroid is computed for each subtriangle
        3. Areas vectors are summed and returned
        4. The centroid returned is obtained from an area-weighted of sum of the subtriangle centroids
=#
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

#=
    Calculates cell volume (scalar) and centroid (vector) from the points and face area vectors (fAVecs) that make up the cell
        fAVecs can be computed using the faceAreaCentroids function

    How it works:
        1. Splits cell into polygonal pyramids, each incorporating a single face and the geometric center of the cell
        2. Computes volume and centroid of each sub-pyramid
        3. Resulting volume is sum of sub-pyramid volumes, centroid is the volume-weighted sum of sub-pyramid centroids
=#
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

# For each face, calculates the vectors from its owner and neighbour cell centers to its own center
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
#=
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
=#
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

#=
    Takes the file path to an OpenFOAM mesh 'points' file, returns a 2D array of point coordinates:
            X   Y   Z
        P1  x1  y1  z1
        P2  x2  y2  z2
        ...
=#
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

#=
    Takes the file path to an OpenFOAM mesh 'faces' file, returns a 1D array of 1D arrays.
    Each entry in the upper level array represents a single face
    Each face-array contains the indices of the points that make up that face
    Note that all indices from the file are incremented by 1, to switch from OpenFOAM's 0-based indexing to Julia's 1-based indexing

    Example return value:
        [ [face1_pt1Index, face1_pt2Index, face1_pt3Index], [face2_pt1Index, face2_pt2Index, face2_pt3Index, face2_pt4Index], ... ]
=#
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
        nPts = length(ptNumbers)
        for p in 1:nPts
            # Add one to convert to 1-based indexing
            push!(facePts, parse(Int64, ptNumbers[p])+1)
        end
        push!(faces, facePts)
    end

    return faces
end

#=
    Takes the file path to an OpenFOAM 'owner' file, returns a 1D array of owner cell indices
    Example return value:
        [ face1_ownerCellIndex, face2_ownerCellIndex, ...  ]
=#
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

#=
    Takes the file path to an OpenFOAM 'neighbour' file, returns a 1D array of neighbour cell indices
    Note that all indices from the file are incremented by 1, to switch from OpenFOAM's 0-based indexing to Julia's 1-based indexing

    Example return value:
        [ face1_neighbourCellIndex, face2_neighbourCellIndex, ...  ]
=#
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

#=
    Returns index of first line containing str
    Otherwise returns -1
=#
function findInLines(str, lines, startLine)
    nLines = size(lines, 1)
    for i in startLine:nLines
        if occursin(str, lines[i])
            return i
        end
    end
    return -1
end

#=
    Takes the file path to an OpenFOAM 'boundary' file, returns 3 arrays:
    boundaryNames:          Array of boundary names (strings)
    boundaryNumFaces:       Array containing the number of faces in each boundary
    boundaryStartFaces:     Array containing the index of the first face in the boundary (boundaries always occupy a continuous sequence of face indices, and are usually numbered last)
=#
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

#=
    Input: Path to an OpenFOAM mesh FOLDER.
    Output: Calls readOFPoints/Faces/Owner/NeighbourFile and returns all of their results (reads all OpenFOAM's mesh file data into arrays)
=#
function readOpenFOAMMesh(polyMeshPath)
    println("Reading mesh: $polyMeshPath")
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

mutable struct Cell
    faceIndices::Vector{Int64}
    pointIndices::Vector{Int64}
end

#=
    Function used to find all the points in each cell.
    Not required for CFD, which is face-based, but required for cell-based .vtk file output.

    Input: Path to an OpenFOAM mesh FOLDER.
    Returns:
        points:         Array of points obtained from readOFPoints()
        cellPtIndices:  Array of arrays, where each subarrary represents a cell, and each entry in a cell's array is the index of one of its points
=#
function OpenFOAMMesh_findCellPts(polyMeshPath)
    pointLocations, pointIndicesByFace, faceOwnerCellIndices, faceNeighborCellIndices, boundaryNames, boundaryNumFaces, boundaryStartFaces = readOpenFOAMMesh(polyMeshPath)
    nCells = maximum(faceOwnerCellIndices)
    nFaces = size(pointIndicesByFace, 1)
    nBoundaries = size(boundaryNames, 1)

    # Output
    cells = Array{Cell, 1}(undef, nCells)
    for i in eachindex(cells)
        cells[i] = Cell(Vector{Int64}(undef, 0), Vector{Int64}(undef, 0))
    end

    ### Populate the cells array ###
    function addCellFaceIndices(adjacentCellList)
        nFaces = length(adjacentCellList)
        for f in 1:nFaces
            ownerCellIndex = adjacentCellList[f]
            push!(cells[ownerCellIndex].faceIndices, f)
        end
    end

    # Start by a just adding face indices, want to order the points appropriately before adding them
    addCellFaceIndices(faceOwnerCellIndices)
    addCellFaceIndices(faceNeighborCellIndices)

    function addAllNewPoints!(cell::Cell, faceIndex)
        for pointIndex in pointIndicesByFace[face2]
            if !any(x->x==pointIndex, cell.pointIndices)
                push!(cell.pointIndices, pointIndex)
            end
        end
    end

    function disjoint(f1Points, f2Points)
        allPoints = vcat(f1Points, f2Points)
        pointSet = Set(allPoints)
        return length(allPoints) == length(pointSet)
    end

    function intersection(points1, points2)
        return intersect(Set(points1), Set(points2))
    end

    function index(element, array)
        for i in eachindex(array)
            if element == array[i]
                return i
            end
        end
        return -1
    end            

    function addEdges!(cellPoints, facePoints, endFacePoints, endFacePointsSet, unusedFaces, pointOffset=4)
        oppositeFaceIndex = -1
        for i in 1:length(unusedFaces)
            fi = unusedFaces[i]
            fiPoints = pointIndicesByFace[fi]
            commonPoints = intersection(fiPoints, facePoints)

            if length(commonPoints) == 2
                # These two points form an edge connecting faces one and two
                points = collect(commonPoints)

                if points[1] in endFacePointsSet
                    p1 = points[1]
                    p2 = points[2]
                else
                    p1 = points[2]
                    p2 = points[1]
                end

                position = index(p1, endFacePoints) # Find position over the point on the face one side
                if position == -1
                    throw(ErrorException("Point $p1 not found in $endFacePoints"))
                end
                cellPoints[ position + pointOffset ] = p2 # The point on the face two side goes in the matching spot
            else
                oppositeFaceIndex = i # This is the face opposite to the one containing facePoints
            end
        end
        return oppositeFaceIndex
    end

    function populatePointIndices_Tet!(cell::Cell)
        # For a tetrahedron, the order of points is unimportant
        # And all points from an arbitrary face
        addAllNewPoints!(cell, cell.faceIndices[1])

        # Choose another face, it will contain the fourth point we need to complete the tetrahedron
        addAllNewPoints!(cell, cell.faceIndices[2])
    end

    function populatePointIndices_Pyramid!(cell::Cell)
        # Pyramids need to have the points that make up their square base ordered 1-4, with the tip of the pyramid coming last
        # Find the quad face
        for faceIndex in cell.faceIndices
            if length(pointIndicesByFace[faceIndex]) == 4
                addAllNewPoints!(cell, faceIndex)

                # There is only one square face, and all other faces contain the pyramid tip
                # Find any other face and add its points to make the pyramid tip the final point
                if faceIndex != cell.faceIndices[1]
                    addAllNewPoints!(cell, cell.faceIndices[1])
                else
                    addAllNewPoints!(cell, cell.faceIndices[2])
                end

                # All done
                break
            end
        end
    end

    function populatePointIndices_Wedge!(cell::Cell)
        # Wedges need the triangular faces numbered 1-3 and 4-6 respectively, where 1 is aligned with 4, 2 with 5, and 3 with 6
        # Using similar strategy it to that used for hex cells: Identify two end faces, then use the faces that connect them to ensure their points are aligned
        
        # First identify the two triangular faces, keep track of one and get rid of the other
        unusedFaces = deepcopy(cell.faceIndices)
        t1 = -1 # Triangular face one, we will base the ordering off this one

        for j in 1:2
            for i in eachindex(unusedFaces)
                faceIndex = unusedFaces[i]
                
                if length(pointIndicesByFace[faceIndex]) == 3
                    # Found a triangular face
                    if t1 == -1
                        t1 = faceIndex
                    end
                    deleteat!(unusedFaces, i)
                    break
                end
            end
        end

        if length(unusedFaces) != 3
            throw(ErrorException("Expected three quadrilateral faces remaining after removal of two triangular end faces, got $(length(unusedFaces))"))
        end

        t1Points = pointIndicesByFace[t1]
        t1PointsSet = Set(t1Points)
        cell.pointIndices = vcat(t1Points, [0, 0, 0]) # Still have to determine proper ordering of points in the other triangular end face, 0's will be replaced below

        # Now all the remaining faces (3) will be quadrilaterals connecting the two triangular end faces
        q1 = pop!(unusedFaces) # Select one at random
        q1Points = pointIndicesByFace[q1] 
        addEdges!(cell.pointIndices, q1Points, t1Points, t1PointsSet, unusedFaces, 3) # Use it to populate 2/3 of the unknown points

        q2 = pop!(unusedFaces) # Select another at random
        q2Points = pointIndicesByFace[q2]
        addEdges!(cell.pointIndices, q2Points, t1Points, t1PointsSet, unusedFaces, 3) # Use it to populate the last unknown point
    end

    function populatePointIndices_Hex!(cell::Cell)
        # Need to order point such that the points that make up a single quadrilateral face are ordered 1-4
            # The face opposite to it must have its points ordered 5-8, with point 1 connected to/aligned with point 5, 2 with 6, etc...
            # We will accomplish this by first picking a starting face at random, this face's points will be numbered 1-4
            # Then we will identify its opposite face and remove it
            # Then all the remaining faces form the connections between the first face we picked and its opposite Face
            # Pairs of these faces will contain 0 or 2 points in common
                # When zero, they are opposite each other
                # When two, they are adjacent and the two points they have in common form an edge connecting the two end faces
                    # This indicates that these two points should be aligned, and is used to properly order the remaining points

        # Start by choosing an arbitrary face
        unusedFaces = deepcopy(cell.faceIndices)
        f1 = pop!(unusedFaces)
        f1Points = pointIndicesByFace[f1]
        f1PointsSet = Set(f1Points)

        # Find the opposite face (no points in common) and remove it
        for i in 1:length(unusedFaces)
            fi = unusedFaces[i]
            fiPoints = pointIndicesByFace[fi]

            if disjoint(f1Points, fiPoints)
                deleteat!(unusedFaces, i)
                break
            end
        end

        ### Now find the correct orientation
            # Point f1_i needs to be aligned (spatially) with Point f2_i
            # We can check for this by using one of the other faces as a guide, since it forms part of the connection between f1 and f2
        cell.pointIndices = vcat(f1Points, [0, 0, 0, 0])

        # Pick an arbitrary other face
        f3 = pop!(unusedFaces)
        f3Points = pointIndicesByFace[f3]
        lastFace = addEdges!(cell.pointIndices, f3Points, f1Points, f1PointsSet, unusedFaces)
        
        # Add points from the edges of the face opposite f3
        lastFacePoints = pointIndicesByFace[unusedFaces[lastFace]]
        deleteat!(unusedFaces, lastFace)        
        noResult = addEdges!(cell.pointIndices, lastFacePoints, f1Points, f1PointsSet, unusedFaces)

        if noResult != -1
            throw(ErrorException("Failure to appropriately order hexahedral cell points for writing to .vtk"))
        end
    end

    # Now use the face indices to gather and order each cell's points appropriately
        # Correct point ordering is important for .vtk output: (See figure 2) https://vtk.org/wp-content/uploads/2015/04/file-formats.pdf
    for cell in cells
        nFaces = length(cell.faceIndices)
        if nFaces == 4
            populatePointIndices_Tet!(cell)
        elseif nFaces == 5
            quadFaceCount = 0

            for f in cell.faceIndices
                if length(pointIndicesByFace[f]) == 4
                    quadFaceCount += 1
                end
            end

            if quadFaceCount == 1
                populatePointIndices_Pyramid!(cell)
            elseif quadFaceCount == 3
                populatePointIndices_Wedge!(cell)
            else
                throw(ErrorException("Unrecognized cell type, cell: $cell. Expecting only hex, tet, wedge, or pyramid cells for vtk output."))
            end
        elseif nFaces == 6
            populatePointIndices_Hex!(cell)
        end
    end

    return pointLocations, cells
end

#=
    Reads an OpenFOAM mesh and returns a Mesh object suitable for calculations in JuliaCFD.

    Mesh defined in dataStructures.jl, documented in dataStructureDefintions.md
=#
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
