#= See dataStructuresDefinitions.md for explanation
=#

mutable struct SolutionState
    cellState::Matrix{Float64}
    cellFluxes::Matrix{Float64}
    cellPrimitives::Matrix{Float64}
    fluxResiduals::Matrix{Float64}
    faceState::Matrix{Float64}
    faceFluxes::Matrix{Float64}
end

struct Mesh
    cells::Vector{Vector{Int64}}
    cVols::Vector{Float64}
    cCenters::Vector{Vector{Float64}}
    cellSizes::Matrix{Float64}

    faces::Vector{Vector{Int64}}
    fAVecs::Vector{Vector{Float64}}
    fCenters::Vector{Vector{Float64}}
    boundaryFaces::Vector{Vector{Int64}}
end
