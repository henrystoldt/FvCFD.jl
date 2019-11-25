using Test
include("mesh.jl")
include("test.jl")

@testset "Cross Product" begin
    v1 = [ 1.0, 2.0, 3.0 ]
    v2 = [ 2.0, 2.0, 2.0 ]
    crossProduct = [ -2.0, 4.0, -2.0 ]
    crossProduct2 = crossProd(v1, v2)
    @test almostEqual(crossProduct, crossProduct2)
end;

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