########## Functions specific to assignment Question 2 ###########
function dXdt(X, Y)
    return 1 - (2 + 1)*X + X*X*Y
end

function dYdt(X, Y)
    return 2*X - X*X*Y
end