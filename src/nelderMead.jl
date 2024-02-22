using Statistics

include("helperFunctions.jl")
include("sampleSpace.jl")

"""
    nelderMead(Xk, f::Function, pDict; alpha = 1.0, gamma = 2.0, beta = 0.5, delta = 0.5, verbose::Bool = false)

Implements the Nelder-Mead optimization algorithm, a heuristic search method used for finding the minimum of an objective function in a multidimensional space. The algorithm is well-suited for problems where the objective function does not have derivatives or the derivatives are not reliable.

# Arguments
- `Xk`: Initial simplex, represented as an `n x p` matrix, where `n` is the dimension of the space and `p` is the number of points in the simplex. Each column of `Xk` is a vertex of the simplex in the space.
- `f`: The objective function to be minimized. It must be a function that accepts a vector representing a point in the space and a dictionary `pDict` of parameters, and returns a scalar value representing the function's value at that point.
- `pDict`: A dictionary of parameters that will be passed to the objective function `f` at each evaluation. This allows for flexible specification of additional information required by `f`.
- `alpha`, `gamma`, `beta`, `delta`: Coefficients that control the reflection (`alpha`), extension (`gamma`), contraction (`beta`), and shrinking (`delta`) operations of the algorithm. These parameters allow the algorithm to be tuned for specific types of objective functions or optimization problems.
- `verbose`: A boolean flag that, when set to `true`, enables the printing of detailed information about the algorithm's progress at each step. This is useful for debugging or understanding the algorithm's behavior on specific problems.

# Key Components
- `fbest`: Keeps a record of the best function value found during the algorithm's execution. This is crucial for tracking the progress of the optimization and determining when the algorithm has converged to a solution.
- `Xk`: The current state of the simplex. The simplex is a geometric figure consisting of `p` points in `n`-dimensional space, which the algorithm adjusts iteratively to explore the space and hone in on the minimum of the objective function.
- `actions`: A dictionary that records the number of times each type of action has been taken during the algorithm's execution. The actions include `extend`, `insideContract`, `outsideContract`, `reflect`, `shrink`, `sort`, and `insertIntoSorted`. This provides insights into the behavior of the algorithm and can help in diagnosing issues or adjusting parameters for better performance.

# Returns
- `Xk`: The updated simplex after the execution of the algorithm, representing the final state of the search for the minimum.
- `fbest`: The best function value found by the algorithm, which is the minimum value of the objective function within the explored region of the space.
- `actions`: A detailed record of the actions taken during the algorithm's execution, useful for analysis and debugging.

# Notes
- The initial simplex `Xk` should be chosen carefully, as it can significantly influence the algorithm's efficiency and its ability to find the global minimum.
- The choice of parameters `alpha`, `gamma`, `beta`, and `delta` can affect the convergence rate and stability of the algorithm. These should be adjusted based on the nature of the objective function and the specific optimization problem.
- The `verbose` option can generate a lot of output for problems with many iterations; use it judiciously to avoid overwhelming the console with information.

# Examples
```julia
# Define an objective function
f(x, pDict) = sum((x .- pDict[:offset]).^2)

# Create an initial simplex
Xk = [1.0 2.0 3.0; 4.0 5.0 6.0]

# Define parameters for the objective function
pDict = Dict(:offset => [1.0, 2.0])

# Execute the Nelder-Mead algorithm
result = nelderMead(Xk, f, pDict, verbose=true)

# Access the results
println("Best point: ", result.Xk[:, 1])
println("Best function value: ", result.fb)
println("Actions taken: ", result.actions)
```
"""
function nelderMead(Xk, f::Function, pDict;
    alpha = 1.0,
    gamma = 2.0,
    beta = 0.5,
    delta = 0.5,
    verbose::Bool = false)

    fevals_NM = 0
    n, p = size(Xk)
    if p != n+1
        error("p = $p != n+1 = $(n+1) shouldn't happen in Nelder Mead")
    end
    actions = Dict(:extend => 0, :extensionFailure => 0, 
    :extensionSuccess => 0, :insideContract => 0, 
    :outsideContract => 0, :reflect => 0, :shrink => 0, 
    :sort => 0, :insertIntoSorted => 0)

    action = "unselected"
    # it is assumed that the simplex is sorted, best (most optimal) point first

    xb, xw, xsw = Xk[:, 1], Xk[:, p], Xk[:, p-1]

    xc = vec(mean(Xk[:, 1:p-1], dims=2))

    xr = reflect(xc, xw, alpha=alpha)

    action = "reflect"
    actions[:reflect] += 1

    F_xr, F_xb = f(xr, pDict, getGradientToo = false), f(xb, pDict, getGradientToo = false)
    fevals_NM += 2

    fbest = F_xb

    if F_xr < F_xb
        # better point than the best
        fbest = F_xr
        myprintln(verbose, "Extending good reflection.")
        xe = extend(xc, xw, gamma=gamma)
        actions[:extend] += 1
        F_xe = f(xe, pDict, getGradientToo = false)
        if F_xe < F_xr
            myprintln(verbose, "Extension a success!")
            actions[:extensionSuccess] += 1
            # extended point is even better! 
            # add it to the top of the simplex
            Xk = hcat(xe, Xk[:, 1:p-1])
            fbest = F_xe
            action = "extend"
        else 
            # extended point is a worse point than reflection
            # add reflected point to the top of the list
            myprintln(verbose, "Extension a failure. Reverting to reflection.")
            actions[:extensionFailure] += 1
            Xk = hcat(xr, Xk[:, 1:p-1])
        end
    else
        F_xsw = f(xsw, pDict, getGradientToo = false)
        fevals_NM += 1
        if F_xr < F_xsw 
            # point satisfies minimum required criterion for addition to simplex
            # just select it and move on
            myprintln(verbose, "Reflection better than second-worst point, so adding it to simplex")
            Xk = insertSortedSimplex(Xk, xr, f, pDict)
            fbest = f(Xk[:, 1], pDict, getGradientToo=false)
            actions[:insertIntoSorted] += 1
        else
            F_xw = f(xw, pDict, getGradientToo = false)
            fevals_NM += 1
            if F_xr < F_xw
                # this point is technically better than the worst point in the simplex, but not particularly useful. Can perform some contracts on it.
                myprintln(verbose, "Reflection only better than worst point. Trying outside contract.")
                xoc = outsideContract(xc, xr, beta=beta)
                actions[:outsideContract] += 1
                F_xoc = f(xoc, pDict, getGradientToo = false)
                fevals_NM += 1
                if F_xoc < F_xr
                    myprintln(verbose, "Outside Contract a success! Adding it to simplex.")
                    # choosing outside contracted point (it may or may not be a useful addition to the simplex)
                    action = "outsideContract"
                    Xk = insertSortedSimplex(Xk, xoc, f, pDict)
                    fbest = f(Xk[:, 1], pDict, getGradientToo=false)
                    actions[:insertIntoSorted] += 1
                    # display(Xk)
                else
                    # outside contract didn't help, but choosing reflected point anyway
                    myprintln(verbose, "Outside Contract a failure. Still adding the reflection to simplex.")
                    Xk = hcat(Xk[1:p-1], xr)
                end
            else
                # oh no the point is even worse than the current worst point
                # Let's try an inside contract
                myprintln(verbose, "Worse than the worst point. Trying an Inside Contract.")
                xic = insideContract(xc, xw, beta=beta)
                action = "insideContract"
                actions[:insideContract] += 1
                F_xic = f(xic, pDict, getGradientToo = false)
                fevals_NM += 1
                if F_xic < F_xw
                    myprintln(verbose, "Inside Contract better than worst point, so adding it.")
                    Xk = insertSortedSimplex(Xk, xic, f, pDict)
                    fbest = f(Xk[:, 1], pDict, getGradientToo=false)
                    actions[:insertIntoSorted] += 1
                else
                    # okay no improvment this time
                    # let's shrink the simplex
                    myprintln(verbose, "No improvement. Shrink the simplex.")
                    Xk = shrinkSortedSimplex(Xk, delta=delta)
                    action = "shrink"
                    actions[:shrink] += 1
                    Xk, Fk = sortSimplex(Xk, f, pDict)

                    fbest = Fk[1]
                    actions[:sort] += 1
                end
            end
        end             
    end

    return (Xk=Xk, fb=fbest, actions=actions)
end

function reflect(xc, xw;
    alpha = 1.0)
    xr = xc + alpha*(xc-xw)
    return xr
end

function extend(xc, xw;
    gamma = 2.0)
    xe = xc + gamma*(xc-xw)
    return xe
end

function outsideContract(xc, xr;
    beta = 0.5)
    xoc = xc + beta*(xr-xc)
    return xoc
end

function insideContract(xc, xw;
    beta = 0.5)
    xic = xc + beta*(xw-xc)
    return xic
end

function shrinkSortedSimplex(simplex;
    delta = 0.5)
    println("Shrinking the simplex")
    println("Previous size: $(size(simplex))")
    simplex += delta*(simplex[:, 1] .- simplex[:, :])
    println("New size $(size(simplex))")
    return simplex
end

function insertSortedSimplex(matrix, new_vector, f::Function, pDict)
    # values = [f(matrix[:, i], pDict, getGradientToo = false) for i in 1:size(matrix, 2)]
    # values = [f(col, pDict, getGradientToo=false) for col in eachcol(matrix)]
    values = [f(collect(col), pDict, getGradientToo=false) for col in eachcol(matrix)]

    new_value = f(new_vector, pDict, getGradientToo = false)

    # Find the insertion index
    insert_index = searchsortedfirst(values, new_value)

    # Insert the new vector and then trim the last column if necessary
    # println("Size of matrix before insertion = $(size(matrix))")
    matrix .= hcat(matrix[:, 1:insert_index-1], new_vector, matrix[:, insert_index:end-1])
    # println("Size of matrix after insertion = $(size(matrix))")
    return matrix
    
end

"""
    sortSimplex(simplex, f::Function, pDict) -> Tuple{Matrix, Vector}

Sorts a simplex based on the evaluation of a provided function `f` on each point (column) of the simplex, and returns both the sorted simplex and the corresponding sorted function values.

# Arguments
- `simplex::Matrix`: A matrix representing the simplex, where each column is a point in the simplex.
- `f::Function`: A function to be evaluated at each point of the simplex. The function should take two arguments: a point (column of the simplex) and `pDict`.
- `pDict`: Parameters to be passed to the function `f` along with each point of the simplex.

# Returns
- `Tuple{Matrix, Vector}`: A tuple containing the sorted simplex (Matrix) and the sorted function values (Vector).

# Details
The function evaluates `f` once for each point in the simplex to obtain the function values. It then sorts the simplex based on these values. The sorting is done by computing a permutation of indices that sorts the function values, and then applying this permutation to both the simplex and the function values array. This method ensures that the simplex points are aligned with their corresponding function evaluations in ascending order.

The function directly returns the sorted simplex and the sorted function values without requiring additional evaluations of `f` beyond the initial computation for each point.

# Example
```julia
sortedSimplex, F = sortSimplex(simplex, myFunction, pDict)
```
"""
function sortSimplex(simplex, f::Function, pDict)

    n, p = size(simplex)
    # Create an array to store the function values for each column of the simplex
    fValues = zeros(p)

    # Evaluate the function `f` for each column and store the results
    for i in 1:p
        fValues[i] = f(simplex[:, i], pDict, getGradientToo=false)
    end

    # Sort the simplex based on the evaluated function values
    # `sortperm` gives the permutation of indices that will sort the array
    sortedIndices = sortperm(fValues)

    # Apply the sorted indices to the simplex columns
    sortedSimplex = simplex[:, sortedIndices]

    # Sort the function values array using the sorted indices
    sortedFValues = vec(fValues[sortedIndices])

    # Return both the sorted simplex and the sorted function values
    return sortedSimplex, sortedFValues
end


"""
    simplexDiameter(simplex) -> Float64

Calculates the diameter of a simplex, defined as the maximum Euclidean distance between any two points (vertices) within the simplex.

# Arguments
- `simplex::Matrix`: A matrix representing the simplex, where each column corresponds to a point (vertex) in the simplex. The matrix should have dimensions `n x p`, where `n` is the dimensionality of the space the simplex resides in, and `p` is the number of points in the simplex.

# Returns
- `Float64`: The diameter of the simplex, which is the greatest distance between any two points in the simplex.

# Details
The function iterates over all unique pairs of points in the simplex, calculating the Euclidean distance between each pair. It keeps track of the maximum distance found during these comparisons, which is returned as the diameter of the simplex.

This measurement can be particularly useful in optimization algorithms and numerical methods that utilize simplices (e.g., the Nelder-Mead algorithm), as it provides a quantitative measure of the simplex's size. A decreasing diameter over iterations can indicate convergence towards a solution.

# Example
```julia
simplex = [1 2 3; 4 5 6]  # A 2x3 matrix representing a simplex in 2D space
diameter = simplexDiameter(simplex)
```
"""
function simplexDiameter(simplex)
    n, p = size(simplex)
    maxDiameter = 0.0
    for i in 1:p-1
        for j in i+1:p
            dist = norm(simplex[:, i] - simplex[:, j])
            maxDiameter = max(maxDiameter, dist)
        end
    end
    return maxDiameter
end

function createInitialSimplexFromOnePoint(x0; 
    deviationFactor=0.5)
    n = length(x0)
    # Generate n x (n+1) matrix using Halton sequence
    haltonMatrix = sampleSpaceHalton(n, n + 1)  # Assuming this function exists

    # Adjust the matrix to fit the deviation factor and center around x0
    adjustedSimplex = zeros(n, n + 1)
    for i in 1:n
        for j = 1:(n+1)
            # Scale and shift the Halton sequence points
            adjustedSimplex[i, j] = x0[i] * (1 - deviationFactor) + (x0[i] * 2 * deviationFactor) * haltonMatrix[i, j]
        end
    end

    # Ensure that x0 is included in the simplex
    adjustedSimplex[:, end] = x0  # Setting the last column to x0

    return adjustedSimplex
end

# function f(v)
#     x, y = v[1], v[2]
#     f = (x-1)*x + (y+1)*y
#     return f
# end

# simplex = Float64.([0 1 3;0 2 4])
# # Create a sample simplex matrix
# # simplex = Float64.([1 3 5; 2 4 6])  # 2x3 matrix, with vectors [1, 2], [3, 4], [5, 6] as columns
# println("originalSimplex:")
# display(simplex)
# v3 = [9, 12]
# # Define a new vector to insert
# new_vector = [2.5, 3.5]  # A new vector

# # Use the function to insert the new vector
# # insertSortedSimplex!(simplex, new_vector, f)
# v2 = [0.5, 2.3]
# # Display the updated simplex
# # insertSortedSimplex!(simplex, v2, f)
# # insertSortedSimplex!(simplex, v3, f)
# verbose = true
# # verbose = false
# simplex = insertSortedSimplex(simplex, [-0.75, -0.5], f)
# println("Latest Simplex:")
# display(simplex)
# simplex, action = nelderMead(simplex, f, verbose=verbose)
# println("Action performed: $(action)")
# println("Latest Simplex:")
# display(simplex)