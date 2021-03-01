using Test
include("../mesh.jl")
include("test.jl")
include("testMeshes.jl")

@testset "Face Geometry" begin
    # Test for a triangle
    pt1 = [ 0.0, 0.0, 0.0 ]
    pt2 = [ 1.0, 0.0, 0.0 ]
    pt3 = [ 0.0, 1.0, 0.0 ]
    points = [ pt1, pt2, pt3 ]

    centroid = [ 0.3333333333333, 0.333333333333, 0.0 ]
    centroid2 = triangleCentroid(points)
    @test almostEqual(centroid, centroid2, 9)

    area = [ 0.0, 0.0, 0.5 ]
    area2 = triangleArea(points)
    @test almostEqual(area, area2)

    area2, centroid2 = faceAreaCentroid(points)
    @test almostEqual(centroid, centroid2, 9)
    @test almostEqual(area, area2)

    # Test for a square
    pt4 = [ 1.0, 1.0, 0.0 ]
    points = [ pt1, pt2, pt4, pt3 ]
    area = [ 0.0, 0.0, 1.0 ]
    centroid = [ 0.5, 0.5, 0.0 ]
    area2, centroid2 = faceAreaCentroid(points)
    @test almostEqual(centroid, centroid2, 9)
    @test almostEqual(area, area2)

    # Test mesh from Moukalled, pg. 247
    cells, cVols, cCenters, faces, fAVecs, fCenters, boundaryFaces = meshMoukalled813()

    # Moukalled pg. 249
    @test almostEqual(cVols[1], 8.625)
    # Calculated online (the one in book is incorrect)
    @test almostEqual(cCenters[1], [1.754, 1.92, 0.5], 2)

    # Moukalled pg. 249
    fAVecs2 = [ [2.5, 1.0, 0.0], [-1.0, 2.0, 0.0], [-2.0, 0.5, 0.0], [-1.0, -2.0, 0.0], [1.5, -1.5, 0.0 ] ]
    for i in 1:5
        @test almostEqual(fAVecs[i], fAVecs2[i])
    end

    # @test almostEqual(mag(fAVecs[6]), 8.625)
    @test almostEqual(mag(fAVecs[7]), 8.625)
    @test almostEqual(mag(fAVecs[6]), 8.625)
end;

@testset "Cell Geometry" begin
    # Test for a cube
    pt1 = [ 0.0, 0.0, 0.0 ]
    pt2 = [ 1.0, 0.0, 0.0 ]
    pt3 = [ 1.0, 1.0, 0.0 ]
    pt4 = [ 0.0, 1.0, 0.0 ]
    pt5 = [ 0.0, 0.0, 1.0 ]
    pt6 = [ 1.0, 0.0, 1.0 ]
    pt7 = [ 1.0, 1.0, 1.0 ]
    pt8 = [ 0.0, 1.0, 1.0 ]

    points = [ pt1, pt2, pt3, pt4, pt5, pt6, pt7, pt8 ]
    faceCentroids = [ [0.5, 0.5, 0.0], [0.0, 0.5, 0.5], [ 0.5, 0.0, 0.5], [0.5, 0.5, 1], [1.0, 0.5, 0.5], [0.5, 1, 0.5] ]
    fAVecs = [ [0.0, 0.0, -1.0], [ -1.0, 0.0, 0.0 ], [ 0.0, -1.0, 0.0 ], [ 0.0, 0.0, 1.0 ], [ 1.0, 0.0, 0.0 ], [ 0.0, 1.0, 0.0 ] ]

    centroid = [ 0.5, 0.5, 0.5 ]
    volume = 1.0
    volume2, centroid2 = cellVolCentroid(points, fAVecs, faceCentroids)
    @test almostEqual(volume, volume2)
    @test almostEqual(centroid, centroid2)

end;

@testset "Reading OF Meshes" begin
    faces = readOFFacesFile("Test/OFshockTube_100/faces")
    @test length(faces) == 501
    @test faces[1] == [2, 103, 305, 204]

    points = readOFPointsFile("Test/OFshockTube_100/points")
    @test length(points)/3 == 404
    @test points[1,:] == [0, -1, -1]
end

@testset "Getting cell points" begin
    points, cells = OpenFOAMMesh_findCellPts("Test/OFshockTube_100")
    @test cells[1].pointIndices == [ 1, 102, 103, 2, 203, 304, 305, 204 ]
    @test cells[2].pointIndices == [ 2, 103, 305, 204, 3, 104, 306, 205 ]
end
