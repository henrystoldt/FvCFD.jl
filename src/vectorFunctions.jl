######################### Vector Functions ########################
# Assumes sizes are equal
function dot(vec1, vec2)
    sum = 0.0
    @fastmath for i in eachindex(vec1)
        sum += vec1[i]*vec2[i]
    end
    return sum
end

# Assumes sizes are equal
function dot(vec1::Vector{Float64}, vec2::Vector{Float64})
    sum = 0.0
    @fastmath for i in eachindex(vec1)
        sum += vec1[i]*vec2[i]
    end
    return sum
end

function dot(vec::Vector{Float64}, matrix::Vector{Vector{Float64}})
    nRows = size(matrix, 1)
    result = Array{Float64, 1}(undef, nRows)
    @fastmath for i in 1:nRows
        result[i] = dot(vec, matrix[i])
    end
    return result
end

function dot(vec::Vector{Float64}, matrix::Matrix{Float64})
    nRows = size(matrix, 1)
    result = Vector{Float64}(undef, nRows)
    @fastmath for i in 1:nRows
        @views result[i] = dot(vec, matrix[i, :])
    end
    return result
end

function cross(v1::Array{Float64, 1}, v2::Array{Float64, 1})
    x = v1[2]*v2[3] - v1[3]*v2[2]
    y = -(v1[1]*v2[3] - v1[3]*v2[1])
    z = v1[1]*v2[2] - v1[2]*v2[1]
    return [x,y,z]
end

function mag(vec)
    #Returns 2-norm of vector
    sqrSum = 0.0
    @fastmath for i = 1:size(vec,1)
        sqrSum += vec[i]*vec[i]
    end
    return sqrt(sqrSum)
end

function normalize(vec)
    #Turns a vector into a unit vector
    @fastmath return vec ./ mag(vec)
end
