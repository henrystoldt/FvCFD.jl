# For numbers of a magnitude less than ~1000
function almostEqual(iterable1, iterable2, allowedError=0.0000000000001)    
    if size(iterable1, 1) != size(iterable2, 1)
        s1 = size(iterable1, 1)
        s2  = size(iterable2, 1)
        println("Size mismatch: $s1 vs $s2")
        println(iterable1)
        println(iterable2)
        return false
    end

    if ndims(iterable1) >= 2
        return matrixAlmostEqual(iterable1, iterable2, allowedError)
    end

    for i in 1:size(iterable1, 1)
        if abs(iterable1[i] - iterable2[i]) > allowedError
            error = iterable1[i] - iterable2[i]
            println("Error: mismatch of $error at position $i")
            println(iterable1)
            println(iterable2)
            return false
        end
    end

    return true
end

function matrixAlmostEqual(m1, m2, allowedError=0.0000000000001)
    s1 = size(m1, 1)

    if size(m1) != size(m2)
        return false
    end

    for i in 1:s1
        if almostEqual(m1[i], m2[i], allowedError) != true
            return false
        end
    end

    return true
end