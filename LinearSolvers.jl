__precompile__()

############# Interactivity ###############
#Function to get matrix from standard input
function getMatrix()
    println("Note: Enter the augmented matrix")
    print("Enter number rows: ")
    nRows = chomp(readline())
    nRows = parse(Int64, nRows)

    matrix = []
    for i in 1:nRows
        print("Enter Row $i: ")
        row = chomp(readline())
        row = split(row)
        
        if i == 1
            matrix = Array{Float64, 2}(undef, nRows, length(row))
        end

        for j in 1:length(row)
            matrix[i,j] = parse(Float64, row[j])
        end        
    end
    return matrix
end

function printMatrix(matrix, eqnRow=nothing)
    nRows = size(matrix, 1)
    
    if eqnRow == nothing
        # No effect on row ordering if one isn't passed in
        eqnRow = Array(1:nRows)
    end

    for i in 1:nRows
        #Print reo-ordered matrix
        row = eqnRow[i]
        println(matrix[row,:])
    end
    println("")
end

############ Gauss Elimination ###############
#Get max coefficient magnitude for each row
function getMaxMagnitudes(matrix)
    nRows = size(matrix, 1)
    maxElements = [ maximum(matrix[i, :]) for i in 1:nRows ]
    minElements = [ minimum(matrix[i, :]) for i in 1:nRows ]
    maxMagnitudes = [ max(abs(maxElements[i]), abs(minElements[i])) for i in 1:nRows ]
    return maxMagnitudes
end

#Choose pivot rows - scaled partial pivoting
function scaleRowsForColumn(matrix, col, maxMagnitudes, eqnRow, currRow)
    nRows = size(matrix, 1)
    maxScaledMag = 0
    maxRow = -1
    for origRow in currRow:nRows
        row = eqnRow[origRow]
        scaledMag = abs(matrix[row, col] / maxMagnitudes[row])
        if scaledMag > maxScaledMag
            maxScaledMag = scaledMag
            maxRow = origRow
        end
    end

    if maxScaledMag == 0
        error("There is no partial pivot possible to make diagonal elements non-zero. Error encountered in column $col")
    end

    return maxRow
end

#eqnRow isa list used to map row indices (which are pivoted) to their actual location in the matrix
function GaussElim_Eliminate!(matrix, eqnRow)
    nRows = size(matrix, 1)
    nCols = size(matrix, 2)

    # Exchange two rows by switching elements in the order vector eqnRow
    function pivot(row1, row2)
        temp = eqnRow[row1]
        eqnRow[row1] = eqnRow[row2]
        eqnRow[row2] = temp
    end

    maxMagnitudes = getMaxMagnitudes(matrix)

    # Elimination
    for origRow in 1:nRows
        pivotCol = origRow
        
        rowToPivot = scaleRowsForColumn(matrix, pivotCol, maxMagnitudes, eqnRow, origRow)
        
        if rowToPivot != origRow
            pivot(origRow, rowToPivot)
        end
        
        pivotRow = eqnRow[origRow]

        pivotElement = matrix[pivotRow, pivotCol]

        for rowToEliminate in (origRow+1):nRows
            pivotMultiple = matrix[eqnRow[rowToEliminate], pivotCol] / pivotElement
            for col in origRow:nCols
                matrix[eqnRow[rowToEliminate], col] -= matrix[pivotRow, col] * pivotMultiple
            end
        end
    end

    return matrix
end

function GaussElim_BackSubsitute!(matrix, eqnRow)
    nRows = size(matrix, 1)
    nCols = size(matrix, 2)

    x = zeros(nRows)
    for origRow in nRows:-1:1
        col = origRow
        row = eqnRow[origRow]
        
        #Starting at last row, find the value of each x
        #Can ignore all values except the forcing vector, which we update below to accumulate all known (numerical) values from x's that have already been solved for
        x[origRow] = matrix[row, nCols] / matrix[row, col]

        # Adjust the forcing vector of rows above current
        for adjustRow in (origRow-1):-1:1
            matrix[eqnRow[adjustRow], nCols] -= matrix[eqnRow[adjustRow], col] * x[origRow]
        end
    end
    return x
end

# Pass in augmented matrix
function Solve_GaussElim!(matrix)
    nRows = size(matrix, 1)
    nCols = size(matrix, 2)
    
    if nCols != nRows + 1
        ArgumentError("Expecting augmented matrix, where #Cols = #Rows + 1")
    end

    eqnRow = Array(1:nRows)
    
    matrix = GaussElim_Eliminate!(matrix, eqnRow)
    return GaussElim_BackSubsitute!(matrix, eqnRow)
end

############# Iterative Methods ############
function Solve_GaussSeidel!(matrix, minRes = 0.0000001, iterLimit=1000)
    return Solve_SOR!(matrix, 1, minRes, iterLimit)
end

function Solve_SOR!(matrix, omega=1.5, minRes=0.0000001, iterLimit=1000)
    nRows = size(matrix, 1)
    nCols = size(matrix, 2)

    # Initial guess is all zeros
    x = zeros(nRows)

    maxResidual = 1
    iterationCounter = 0
    while maxResidual > minRes && iterationCounter < iterLimit 
        iterationCounter += 1
        maxResidual = 0

        for row in 1:nRows
            residual = matrix[row, nCols]
            for col in 1:nCols-1
                residual -= matrix[row, col] * x[col]
            end
            if abs(residual) > maxResidual
                maxResidual = abs(residual)
            end

            x[row] += omega * residual / matrix[row, row]
        end
    end

    if iterationCounter == iterLimit
        println("ERROR: Convergence not achieved in $iterLimit iterations.")
    else
        for i in 1:length(x)
            if isnan(x[i])
                println("ERROR: Calculation diverged")
                break
            end
        end
    end

    return x
end

############# n-Diagonal Matrix Solver ############
# Pass in an augmented matrix which contains only the x-diagonal elements
# Ex:  [0 2 1 5; (tridiagonal)
      # 1 2 1 9;
      # 1 2 0 12 ]
# Solves using Thomas algorithm for an n-bandwidth matrix, auto-detects bandwidth
function Solve_Diag!(matrix)
    nRows = size(matrix, 1)
    nCols = size(matrix, 2)
    bandwidth = nCols - 1
    center = convert(Int32, nCols / 2)
    
    if bandwidth % 2 != 1
        ArgumentError("Bandwidth should always be an odd number, with the center elements being those which would be on the main diagonal of the full matrix")
    end

    # Eliminate
    elimRows = convert(Int32, (bandwidth - 1) / 2)
    for pivotRow in 1:nRows-1
        pivotElement = matrix[pivotRow, center]
        if pivotRow + elimRows > nRows
            elimRows = nRows - pivotRow
        end

        # Because the rows are offset by 1 for each row, elimRow is also the horizontal offset required
        for offset in 1:elimRows
            pivotMultiple = matrix[pivotRow + offset, center - offset] / pivotElement
            for col in 1:nCols - (1 + offset)
                matrix[pivotRow + offset, col] -= matrix[pivotRow, col + offset] * pivotMultiple
            end
            matrix[pivotRow + offset, nCols] -= matrix[pivotRow, nCols] * pivotMultiple
        end
    end

    # Back substitute
    elimRows = convert(Int32, (bandwidth - 1) / 2)
    x = zeros(nRows)
    for row in nRows:-1:1
        x[row] = matrix[row, nCols] / matrix[row, center]

        #Adjust forcing vector in upper rows to use this new information
        if row - elimRows < 1
            elimRows = row - 1
        end

        for elimRow in 1:elimRows
            matrix[row - elimRow, nCols] -= x[row] * matrix[row - elimRow, center + elimRow]
        end
    end

    return x    
end