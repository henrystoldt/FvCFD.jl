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
function triangleCentroid(points::Array{Array{Float64, 1}})
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
function cellVolCentroid(points::Array{Array{Float64, 1}}, fAVecs::Array{Array{Float64, 1}}, faceCentroids::Array{Array{Float64, 1}})
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
