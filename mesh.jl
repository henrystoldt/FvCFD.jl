# Methods from Moukalled et al. FVM - OpenFOAM, Matlab
include("vectorFunctions.jl")

#3D only
function crossProd(v1::Array{Float64, 1}, v2::Array{Float64, 1})
    x = v1[2]*v2[3] - v1[3]*v2[2]
    y = -(v1[1]*v2[3] - v1[3]*v2[1])
    z = v1[1]*v2[2] - v1[2]*v2[1]
    return [x,y,z]
end

# Points must be ordered sequentially
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

# Points must be ordered sequentially
# Returns face area vector and centroid
# Splits face into subtriangles
    # Area and centroid is computed for each subtriangles
    # Areas vectors are summed and returned
    # The centroid returned is obtainde from an area-weighted of sum of the subtriangle areas
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

# Returns cell volume (scalar) and centroid (vector)
# fAVecs can be computed using the faceAreaCentroids function
# Splits cell into polygonal pyramids, each incorporating a single face and the geometric center of the cell
#   Computes volume and centroid of each sub-pyramid
#   Resulting volume is sum, centroid is the volume-weighted sum
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

function cellCentroidToFaceVec(points::Array{Array{Float64, 1}}, faceCentroids::Array{Array{Float64, 1}}, cellCentroids::Array{Array{Float64, 1}})
    # Calculates vectors from cell centroid to face centers
    #TODO: Understand ordering to make sure we're accepting faces in the order we need to pass them back in??

    nFaces = size(faceCentroids, 1)

    cellToFaceVec = zeros(nFaces, 3)

    for f in 1:nFaces
        cellToFaceVec[f,:] = faceCentroids[f] - cellCentroids
    end

    return cellToFaceVec
end

function unstructuredMeshInfo(mesh)
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = mesh
    nCells = size(cells, 1)
    nFaces = size(faces, 1)
    nBoundaries = size(boundaryFaces, 1)

    # Count boundary faces
    nBdryFaces = 0
    for bdry in 1:nBoundaries
        nBdryFaces += size(boundaryFaces[bdry], 1)
    end

    return nCells, nFaces, nBoundaries, nBdryFaces
end

function isNumber(str)
    re = r"^[+-]?([0-9]+([.][0-9]*)?|[.][0-9]+)$"
    return occursin(re, str)
end

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

# Returns index of first line containing str
# Otherwise returns -1
function findInLines(str, lines, startLine)
    nLines = size(lines, 1)
    for i in startLine:nLines
        if occursin(str, lines[i])
            return i
        end
    end
    return -1
end

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

# Supply absolute path to folder containing the OpenFOAM points, faces, cells, and boundaries files
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

function OpenFOAMMesh_findCellPts(polyMeshPath)
    points, OFfaces, owner, neighbour, boundaryNames, boundaryNumFaces, boundaryStartFaces = readOpenFOAMMesh(polyMeshPath)
    nCells = maximum(owner)
    nFaces = size(OFfaces, 1)
    nBoundaries = size(boundaryNames, 1)

    cellPtIndices = Array{Array{Int64, 1}, 1}(undef, nCells)
    cells = Array{Array{Int64, 1}, 1}(undef, nCells)

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

    for c in 1:nCells
        pts = Array{Int64,1}(undef, 0)

        function addFace(f)
            for pt in OFfaces[f]
                push!(pts, pt)
            end
        end
        function disjoint(f)
            for pt in OFfaces[f]
                if any(x->x==pt, pts)
                    return false
                end
            end
            return true
        end
        addFace(cells[c][5])

        # Add any faces with no points in common with existing faces
        for f in cells[c]
            if disjoint(f)
                addFace(f)
            end
        end

        # Add any remaining points
        for f in cells[c]
            for pt in OFfaces[f]
                if !any(x->x==pt, pts)
                    push!(pts, pt)
                end
            end
        end

        cellPtIndices[c] = pts
    end

    return points, cellPtIndices
end

function OpenFOAMMesh(polyMeshPath)
    points, OFfaces, owner, neighbour, boundaryNames, boundaryNumFaces, boundaryStartFaces = readOpenFOAMMesh(polyMeshPath)
    nCells = maximum(owner)
    nFaces = size(OFfaces, 1)
    nBoundaries = size(boundaryNames, 1)

    cells = Array{Array{Int64, 1}, 1}(undef, nCells)
    cVols = zeros(nCells)
    cCenters = Array{Array{Float64, 1}, 1}(undef, nCells)
    faces = Array{Array{Int64, 1}, 1}(undef, nFaces)
    fAVecs = Array{Array{Float64, 1}, 1}(undef, nFaces)
    fCenters = Array{Array{Float64, 1}, 1}(undef, nFaces)
    boundaryFaces = Array{Array{Int64, 1}, 1}(undef, nBoundaries)

    for f in 1:nFaces
        fPts = [ points[pt,:] for pt in OFfaces[f] ]
        fAVecs[f], fCenters[f] = faceAreaCentroid(fPts)
    end
    #fAVecs and fCenters now complete

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

    # Add boundaries
    for b in 1:nBoundaries
        startF = boundaryStartFaces[b]
        endF = startF + boundaryNumFaces[b] - 1
        boundaryFaces[b] = Array(startF:endF)
    end

    mesh = [ cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces ]

    return mesh
end
