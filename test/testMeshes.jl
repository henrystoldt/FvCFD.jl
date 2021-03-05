using fvCFD

# Creates a simple sample cell from Moukalled with one fully-defined cell and a centroid-only neighbour
# For laplacian test cases
function meshMoukalled813()
    f1Pts = [
        [ 2.5, 4.0, 1.0],
        [ 3.5, 1.5, 1.0],
        [ 3.5, 1.5, 0.0],
        [ 2.5, 4.0, 0.0]
    ]
    f2Pts = [
        [ 0.5, 3.0, 1.0],
        [ 2.5, 4.0, 1.0],
        [ 2.5, 4.0, 0.0],
        [ 0.5, 3.0, 0.0]
    ]
    f3Pts = [
        [ 0.0, 1.0, 1.0],
        [ 0.5, 3.0, 1.0],
        [ 0.5, 3.0, 0.0],
        [ 0.0, 1.0, 0.0]
    ]
    f4Pts = [
        [ 2.0, 0.0, 1.0],
        [ 0.0, 1.0, 1.0],
        [ 0.0, 1.0, 0.0],
        [ 2.0, 0.0, 0.0]
    ]
    f5Pts = [
        [ 2.0, 0.0, 1.0],
        [ 2.0, 0.0, 0.0],
        [ 3.5, 1.5, 0.0],
        [ 3.5, 1.5, 1.0]
    ]
    f6Pts = [
        [ 2.5, 4.0, 1.0],
        [ 0.5, 3.0, 1.0],
        [ 0.0, 1.0, 1.0],
        [ 2.0, 0.0, 1.0],
        [ 3.5, 1.5, 1.0],
    ]
    f7Pts = [
        [ 2.5, 4.0, 0.0],
        [ 3.5, 1.5, 0.0],
        [ 2.0, 0.0, 0.0],
        [ 0.0, 1.0, 0.0],
        [ 0.5, 3.0, 0.0]
    ]
    points = [
        [ 2.5, 4.0, 0.0],
        [ 3.5, 1.5, 0.0],
        [ 2.5, 4.0, 1.0],
        [ 3.5, 1.5, 1.0],
        [ 2.0, 0.0, 0.0],
        [ 0.0, 1.0, 0.0],
        [ 0.0, 1.0, 1.0],
        [ 2.0, 0.0, 1.0],
        [ 0.5, 3.0, 0.0],
        [ 0.5, 3.0, 1.0]
    ]
    cells = []

    facePoints = [ f1Pts, f2Pts, f3Pts, f4Pts, f5Pts, f6Pts, f7Pts ]
    fAVecs = Array{Array{Float64, 1}}(undef, 7)
    fCenters = Array{Array{Float64, 1}}(undef, 7)
    for i in 1:7
        fAVec, fCenter = fvCFD.faceAreaCentroid(facePoints[i])
        fAVecs[i] = fAVec
        fCenters[i] = fCenter
    end
    faces = [ [ 1,2 ] ]
    boundaryFaces = []

    cellVol, cellCentroid = fvCFD.cellVolCentroid(points, fAVecs, fCenters)
    cVols = [ cellVol, 1 ]
    cCenters = [ cellCentroid, [4.25, 3.5, 0.5 ]]

    return [ cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces ]
end
