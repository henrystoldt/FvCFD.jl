using Plots
using PyPlot
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
