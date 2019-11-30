
######################### Vector Functions ########################
# Assumes sizes are equal
function dot(vec1, vec2)
    #Returns dot product of two vectors
    s1 = size(vec1, 1)
    s2 = size(vec2, 1)
    if s1 != s2
        return 0
    end
    sum = 0
    for i in 1:size(vec1, 1)
        sum += vec1[i]*vec2[i]
    end
    return sum
end

function mag(vec)
    #Returns 2-norm of vector
    sqrSum = 0
    for i = 1:size(vec,1)
        sqrSum += vec[i]*vec[i]
    end
    return sqrt(sqrSum)
end

function normalize(vec)
    #Turns a vector into a unit vector
    return vec ./ mag(vec)
end

function copyValues(fromIndex, toIndex, varArrays)
    for varArray in varArrays
        varArray[toIndex] = varArray[fromIndex]
    end
end

# From and to indices are assumed to be the first index
function copyValues(fromIndex, toIndex, varArrays::Array{Float64, 2})
    varArrays[toIndex, :] .= varArrays[fromIndex, :]
end

function setValues(value, indices, varArrays)
    for varArray in varArrays
        for i in indices
            if i <= size(varArray, 1)
                varArray[i] = value
            else
                println("Error, index $i out of bounds for array: $varArray")
            end
        end
    end
end
