############## ddt functions defined in Assignment 3, Q1 ###############
function dxdt(x, y, z)
    return -10*x + 10*y
end

function dydt(x, y, z)
    return -x*z + 28*x - y
end

function dzdt(x, y, z)
    return x*y - (8/3)*z
end

# Modified function with b = 8/3 + 0.001
function dzdt2(x, y, z)
    return x*y - (8/3 + 0.001)*z
end