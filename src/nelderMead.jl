using Statistics

include("helperFunctions.jl")
include("sampleSpace.jl")

"""
    nelderMead(Xk, f::Function, pDict; alpha = 1.0, gamma = 2.0, beta = 0.5, delta = 0.5, verbose::Bool = false)

Implements the Nelder-Mead optimization algorithm, a robust heuristic method designed for minimizing an objective function in a multidimensional space without requiring the function's derivatives. This implementation is particularly suited for scenarios where the objective function is non-linear, non-differentiable, or where derivatives are not practically obtainable.

# Arguments
- `Xk`: The initial simplex, represented as an `n x (n+1)` matrix, where `n` is the dimensionality of the space. The simplex is a geometric figure consisting of `n+1` points, each column representing a vertex in the multidimensional space.
- `f`: The objective function to be minimized. It must accept a vector (representing a point in the space), a dictionary `pDict` of parameters, and a flag `getGradientToo` (which should be set to `false`), returning a scalar value indicative of the function's value at the point.
- `pDict`: A dictionary of parameters that will be passed to the objective function `f` during its evaluation. This allows the user to specify additional details required by `f`.
- `alpha`, `gamma`, `beta`, `delta`: Tunable parameters controlling the reflection (`alpha`), extension (`gamma`), contraction (`beta`), and shrinkage (`delta`) operations of the algorithm. Adjusting these parameters can significantly impact the algorithm's performance and convergence behavior on different types of optimization problems.
- `verbose`: A boolean flag that, when set to `true`, activates verbose logging, providing detailed insights into each step of the algorithm's execution. This is invaluable for debugging or for a detailed examination of the algorithm's path through the solution space.

# Key Components
- `fbest`: Maintains the best (lowest) function value encountered during the execution of the algorithm. This is a critical component for tracking the optimization progress and for determining the convergence of the algorithm.
- `Xk`: Represents the current state of the simplex. The simplex is dynamically adjusted through a series of operations (reflection, extension, contraction, shrinkage) aimed at exploring the solution space and converging towards the minimum of the objective function.
- `actions`: A comprehensive dictionary that tallies the occurrences of each specific action (e.g., `extend`, `reflect`, `shrink`) undertaken during the algorithm's execution. This update introduces finer granularity in tracking, including successful and failed extensions and contractions (`extensionSuccess`, `insideContractFailure`, etc.), providing deeper insights into the algorithm's behavior and the efficiency of different operations.

# Returns
- `Xk`: The final simplex after the completion of the optimization process. This represents the last set of points considered by the algorithm in its search for the minimum.
- `fbest`: The best function value found by the algorithm, indicative of the minimum value of the objective function explored by the algorithm.
- `actions`: A detailed tally of all actions executed during the algorithm's run, offering valuable insights into the optimization process and the algorithm's interaction with the solution space.

# Notes
- The initial simplex `Xk` is pivotal to the algorithm's performance and its ability to converge to a global minimum. A well-chosen initial simplex can significantly enhance the efficiency and outcome of the optimization process.
- The parameters `alpha`, `gamma`, `beta`, and `delta` should be carefully tuned based on the characteristics of the objective function and the specific optimization problem at hand to ensure optimal performance and convergence of the algorithm.
- The verbose output can be quite extensive for problems involving a large number of iterations. It should be used judiciously to avoid clutter and to focus on key stages of the optimization process.

# Examples
```julia
# Define a quadratic objective function
f(x, pDict) = sum((x .- pDict[:offset]).^2)

# Set up an initial simplex
Xk = [1.0 2.0 3.0; 4.0 5.0 6.0]

# Specify parameters for the objective function
pDict = Dict(:offset => [1.0, 2.0])

# Execute the Nelder-Mead optimization
result = nelderMead(Xk, f, pDict, verbose=true)

# Display the optimization results
println("Optimal point: ", result.Xk[:, 1])
println("Optimal function value: ", result.fb)
println("Actions taken: ", result.actions)
```
"""
function nelderMead(Xk, f::Function, pDict;
    alpha = 1.0,
    gamma = 2.0,
    beta = 0.5,
    delta = 0.5,
    verbose::Bool = false)

    fevals_1NM = 0
    n, p = size(Xk)
    if p != n+1
        error("p = $p != n+1 = $(n+1) shouldn't happen in Nelder Mead")
    end
    actions = Dict(:extend => 0, :extensionFailure => 0, 
    :extensionSuccess => 0, :insideContract => 0, :insideContractFailure => 0, :insideContractSuccess => 0,
    :outsideContract => 0, :outsideContractFailure => 0,
    :outsideContractSuccess =>0, :reflect => 0, :shrink => 0, 
    :sort => 0, :insertIntoSorted => 0)

    action = "unselected"
    # it is assumed that the simplex is sorted, best (most optimal) point first

    xb, xw, xsw = Xk[:, 1], Xk[:, p], Xk[:, p-1]

    xc = vec(mean(Xk[:, 1:p-1], dims=2))

    xr = reflect(xc, xw, alpha=alpha)

    action = "reflect"
    actions[:reflect] += 1

    F_xr, F_xb = f(xr, pDict, getGradientToo = false), f(xb, pDict, getGradientToo = false)
    fevals_1NM += 2

    fbest = F_xb

    if F_xr < F_xb
        # better point than the best
        fbest = F_xr
        myprintln(verbose, "Extending good reflection.")
        xe = extend(xc, xw, gamma=gamma)
        actions[:extend] += 1
        F_xe = f(xe, pDict, getGradientToo = false)
        fevals_1NM += 1
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
        fevals_1NM += 1
        if F_xr < F_xsw 
            # point satisfies minimum required criterion for addition to simplex
            # just select it and move on
            myprintln(verbose, "Reflection better than second-worst point, so adding it to simplex")
            Xk = insertSortedSimplex(Xk, xr, f, pDict)
            fevals_1NM += (n+1)
            fbest = f(Xk[:, 1], pDict, getGradientToo=false)
            fevals_1NM += 1
            actions[:insertIntoSorted] += 1
        else
            F_xw = f(xw, pDict, getGradientToo = false)
            fevals_1NM += 1
            if F_xr < F_xw
                # this point is technically better than the worst point in the simplex, but not particularly useful. Can perform some contracts on it.
                myprintln(verbose, "Reflection only better than worst point. Trying outside contract.")
                xoc = outsideContract(xc, xr, beta=beta)
                action = "outsideContract"
                actions[:outsideContract] += 1
                F_xoc = f(xoc, pDict, getGradientToo = false)
                fevals_1NM += 1
                if F_xoc < F_xr
                    myprintln(verbose, "Outside Contract a success! Adding it to simplex.")
                    actions[:outsideContractSuccess] += 1
                    # choosing outside contracted point (it may or may not be a useful addition to the simplex)
                    Xk = insertSortedSimplex(Xk, xoc, f, pDict)
                    fevals_1NM += (n+1)
                    fbest = f(Xk[:, 1], pDict, getGradientToo=false)
                    fevals_1NM += 1
                    actions[:insertIntoSorted] += 1
                    # display(Xk)
                else
                    # outside contract didn't help, but choosing reflected point anyway
                    myprintln(verbose, "Outside Contract a failure. Still adding the reflection to simplex.")
                    actions[:outsideContractFailure] += 1
                    Xk = hcat(Xk[:, 1:p-1], xr)
                end
            else
                # oh no the point is even worse than the current worst point
                # Let's try an inside contract
                myprintln(verbose, "Worse than the worst point. Trying an Inside Contract.")
                xic = insideContract(xc, xw, beta=beta)
                action = "insideContract"
                actions[:insideContract] += 1
                F_xic = f(xic, pDict, getGradientToo = false)
                fevals_1NM += 1
                if F_xic < F_xw
                    myprintln(verbose, "Inside Contract better than worst point, so adding it.")
                    actions[:insideContractSuccess] += 1
                    Xk = insertSortedSimplex(Xk, xic, f, pDict)
                    fevals_1NM += (n+1)
                    fbest = f(Xk[:, 1], pDict, getGradientToo=false)
                    fevals_1NM += 1
                    actions[:insertIntoSorted] += 1
                else
                    # okay no improvment this time
                    # let's shrink the simplex
                    myprintln(verbose, "No improvement. Shrink the simplex.")
                    actions[:insideContractFailure] += 1
                    Xk = shrinkSortedSimplex(Xk, delta=delta)
                    action = "shrink"
                    actions[:shrink] += 1
                    Xk, Fk = sortSimplex(Xk, f, pDict)
                    fevals_1NM += (n+1)
                    fbest = Fk[1]
                    actions[:sort] += 1
                end
            end
        end             
    end

    return (Xk=Xk, fb=fbest, actions=actions, fevals_1NM=fevals_1NM)
end

"""
    reflect(xc, xw; alpha = 1.0)

Reflects the worst point `xw` of the simplex across the centroid `xc` to find a new point `xr`.

# Arguments
- `xc`: The centroid of the simplex excluding the worst point. It is calculated as the mean of all points except `xw`.
- `xw`: The worst point in the simplex, i.e., the point with the highest function value.
- `alpha`: Reflection coefficient. A value of `1.0` reflects the point directly across the centroid. Values greater than `1.0` move the reflected point further away from the centroid.

# Returns
- `xr`: The reflected point.

# Examples
```julia
xc = [1.0, 1.0]
xw = [2.0, 2.0]
xr = reflect(xc, xw, alpha = 1.0)
```
"""
function reflect(xc, xw;
    alpha = 1.0)
    xr = xc + alpha*(xc-xw)
    return xr
end


"""
    extend(xc, xw; gamma = 2.0)

Extends beyond the reflected point to attempt to find a better point `xe` in the direction away from the worst point `xw`.

# Arguments
- `xc`: The centroid of the simplex excluding the worst point.
- `xw`: The worst point in the simplex.
- `gamma`: Extension coefficient, typically greater than `1.0`. Specifies how far beyond the reflection the extension should go.

# Returns
- `xe`: The extended point.

# Examples
```julia
xc = [1.0, 1.0]
xw = [2.0, 2.0]
xe = extend(xc, xw, gamma = 2.0)
```
"""
function extend(xc, xw;
    gamma = 2.0)
    xe = xc + gamma*(xc-xw)
    return xe
end

"""
    outsideContract(xc, xr; beta = 0.5)

Performs an outside contraction between the centroid `xc` and the reflected point `xr`, creating a new point `xoc`.

# Arguments
- `xc`: The centroid of the simplex.
- `xr`: The reflected point.
- `beta`: Contraction coefficient, typically between `0.0` and `1.0`, to control the contraction strength.

# Returns
- `xoc`: The outside contracted point.

# Examples
```julia
xc = [1.0, 1.0]
xr = [0.5, 0.5]
xoc = outsideContract(xc, xr, beta = 0.5)
```
"""
function outsideContract(xc, xr;
    beta = 0.5)
    xoc = xc + beta*(xr-xc)
    return xoc
end


### Function: `insideContract`

"""
    insideContract(xc, xw; beta = 0.5)

Performs an inside contraction between the centroid `xc` and the worst point `xw`, producing a new point `xic`.

# Arguments
- `xc`: The centroid of the simplex.
- `xw`: The worst point in the simplex.
- `beta`: Contraction coefficient, influences the closeness of the contracted point to the centroid.

# Returns
- `xic`: The inside contracted point.

# Examples
```julia
xc = [1.0, 1.0]
xw = [2.0, 2.0]
xic = insideContract(xc, xw, beta = 0.5)
```
"""
function insideContract(xc, xw;
    beta = 0.5)
    xic = xc + beta*(xw-xc)
    return xic
end


"""
    shrinkSortedSimplex(simplex; delta = 0.5)

Shrinks the simplex towards its best point, scaling all points (except the best one) closer to the best point by a factor of `delta`.

# Arguments
- `simplex`: The current simplex, an `n x (n+1)` matrix where `n` is the dimensionality of the space.
- `delta`: Shrinkage coefficient, typically between `0.0` and `1.0`. A value of `0.5` halves the distance from each point to the best point.

# Returns
- The shrunk simplex.

# Examples
```julia
simplex = [1.0 2.0 3.0; 4.0 5.0 6.0]
shrinkSortedSimplex(simplex, delta = 0.5)
```
"""
function shrinkSortedSimplex(simplex;
    delta = 0.5)

    simplex += delta*(simplex[:, 1] .- simplex[:, :])

    return simplex
end

"""
    insertSortedSimplex(matrix, new_vector, f::Function, pDict)

Inserts a new vector into a sorted simplex matrix according to the objective function values, removing the worst point to maintain the simplex size.

# Arguments
- `matrix`: The current simplex, represented as an `n x (n+1)` matrix, where `n` is the dimensionality of the space. The simplex is assumed to be sorted in ascending order of the objective function values, with the first column being the best point and the last column the worst.
- `new_vector`: The new point to be inserted into the simplex. It is a vector of the same dimensionality as the other points in the simplex.
- `f`: The objective function used to evaluate the points in the simplex. It must accept a vector (representing a point in the space), a dictionary `pDict` of parameters, and a flag `getGradientToo` (which should be set to `false`), returning a scalar value indicative of the function's value at the point.
- `pDict`: A dictionary of parameters that will be passed to the objective function `f` during its evaluation.

The function first evaluates the objective function values for all columns in the existing simplex and for the new vector. It then finds the appropriate insertion index for the new vector based on its objective function value, ensuring that the simplex remains sorted. The new vector is inserted into the matrix at the found index, and the worst point (previously the last column) is removed.

# Returns
- The updated simplex matrix with the new vector inserted and the previous worst point removed.

# Examples
```julia
# Define the objective function
f(x, pDict) = sum((x .- pDict[:offset]).^2)

# Current simplex and new vector
matrix = [1.0 2.0 3.0; 4.0 5.0 6.0]
new_vector = [2.5, 5.5]

# Parameters for the objective function
pDict = Dict(:offset => [1.0, 2.0])

# Insert the new vector into the simplex
updated_simplex = insertSortedSimplex(matrix, new_vector, f, pDict)
```
"""
function insertSortedSimplex(matrix, new_vector, f::Function, pDict)

    values = [f(collect(col), pDict, getGradientToo=false) for col in eachcol(matrix)]

    new_value = f(new_vector, pDict, getGradientToo = false)

    # Find the insertion index
    insert_index = searchsortedfirst(values, new_value)

    # Insert the new vector and then trim the last column (worst point)
    matrix .= hcat(matrix[:, 1:insert_index-1], new_vector, matrix[:, insert_index:end-1])
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

"""
    createInitialSimplexFromOnePoint(x0; deviationFactor=0.5)

Generates an initial simplex around a given point `x0` using the Halton sequence for optimization algorithms. The generated simplex is adjusted to encompass a region defined by a `deviationFactor`, ensuring diverse starting points for the optimization.

# Arguments
- `x0`: The initial point from which to generate the simplex. It is a vector representing a position in the parameter space.
- `deviationFactor`: A scaling factor that determines the spread of the simplex points around `x0`. A higher deviation factor results in a wider spread of the simplex points.

The simplex is generated as an `n x (n+1)` matrix, where `n` is the dimensionality of the space derived from the length of `x0`. Each column of the matrix represents a vertex of the simplex. The vertices are generated using the Halton sequence to ensure low-discrepancy sampling of the parameter space, then scaled and shifted to create a simplex around `x0`.

# Returns
- An `n x (n+1)` matrix representing the initial simplex, with the last column explicitly set to `x0` to include it as a vertex of the simplex.

# Examples
```julia
# Initial point in a 2-dimensional space
x0 = [1.0, 2.0]

# Generate the initial simplex
initialSimplex = createInitialSimplexFromOnePoint(x0, deviationFactor=0.5)
```
"""
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