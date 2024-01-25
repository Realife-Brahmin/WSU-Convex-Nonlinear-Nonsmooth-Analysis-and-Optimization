using Statistics

include("helperFunctions.jl")
include("sampleSpace.jl")

function nelderMead(simplex, f::Function;
    alpha = 1.0,
    gamma = 2.0,
    beta = 0.5,
    delta = 0.5,
    verbose::Bool = false)
    n, p = size(simplex)
    action = "unselected"
    # it is assumed that the simplex is sorted, best (most optimal) point first
    xb, xw, xsw = simplex[:, 1], simplex[:, p], simplex[:, p-1] 
    xc = mean(simplex[:, 1:p-1], dims=2)
    xr = reflect(xc, xw, alpha=alpha)
    action = "reflect"
    F_xr, F_xb = f(xr), f(xb)
    if F_xr < F_xb
        # better point than the best
        xe = extend(xc, xw, gamma=gamma)
        F_xe = f(xe)
        if F_xe < F_xr
            # extended point is even better! 
            # add it to the top of the simplex
            simplex = hcat(xe, simplex[:, 2:p-1])
            action = "extend"
        else 
            # extended point is a worse point than reflection
            # add reflected point to the top of the list
            simplex = hcat(xr, simplex[:, 2:p-1])
        end
    else
        F_xsw = f(xsw)
        if F_xr < F_xsw 
            # point satisfies minimum required criterion for addition to simplex
            # just select it and move on
            !insertSortedSimplex!(simplex, xr, f)
        else
            F_xw = f(xw)
            if F_xr < F_xw
                # this point is technically better than the worst point in the simplex, but not particularly useful. Can perform some contracts on it.
                xoc = outsideContract(xc, xr, beta=beta)
                F_xoc = f(xoc)
                if F_xoc < F_xr
                    # choosing outside contracted point (it may or may not be a useful addition to the simplex)
                    action = "outsideContract"
                    insertSortedSimplex!(simplex, xoc, f)
                else
                    # outside contract didn't help, but choosing reflected point anyway
                    simplex = hcat(simplex[1:p-1], xr)
                end
            else
                # oh no the point is even worse than the current worst point
                # Let's try an inside contract
                xic = insideContract(xc, xw, beta=beta)
                action = "insideContract"
                F_xic = f(xic)
                if F_xic < F_xw
                    # technically okay
                    insertSortedSimplex!(simplex, xic, f)
                else
                    # okay no improvment this time
                    # let's shrink the simplex
                    shrinkSortedSimplex!(simplex, delta=delta)
                    action = "shrink"
                end
            end
        end             
    end

    return (simplex=simplex, action=action)
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

function shrinkSortedSimplex!(simplex;
    delta = 0.5)
    simplex += delta*(simplex[1, :] .- simplex[:, :])
    return simplex
end

"""
    insertSortedSimplex!(matrix::Matrix{T}, new_vector::Vector{T}, f::Function) where T

Insert a vector into a simplex (represented as a matrix with vectors as columns) 
such that the columns remain sorted based on the value of a given function `f`, 
applied to these vectors. After insertion, the last column (represening the worst point) is removed to maintain a constant size.

# Arguments
- `matrix::Matrix{T}`: The simplex matrix with vectors as columns.
- `new_vector::Vector{T}`: The vector to be inserted into the simplex.
- `f::Function`: A function that is applied to vectors to determine their order.

# Returns
- `Matrix{T}`: The updated simplex matrix with the new vector inserted and the 
last column removed if necessary to maintain size.

# Examples
```julia
# Define a function f() to apply to vectors (e.g., sum of elements)
f(v) = sum(v)

# Create a simplex matrix
simplex = [1 3 5; 2 4 6]  # 2x3 matrix, with columns as vectors

# Define a new vector to insert
new_vector = [7, 8]  # The new vector to be inserted

# Insert the new vector into the simplex
updated_simplex = insertSortedSimplex!(simplex, new_vector, f)

# The updated simplex is now [1 3 7; 2 4 8] if f is sum of elements
```
"""
function insertSortedSimplex!(matrix, new_vector, f::Function)
    values = [f(matrix[:, i]) for i in 1:size(matrix, 2)]
    new_value = f(new_vector)

    # Find the insertion index
    insert_index = searchsortedfirst(values, new_value)

    # Insert the new vector and then trim the last column if necessary
    new_matrix = hcat(matrix[:, 1:insert_index-1], new_vector, matrix[:, insert_index:end])
    if size(new_matrix, 2) > size(matrix, 2)
        return new_matrix[:, 1:end-1]  # Keep the matrix size constant by trimming the last column
    else
        return new_matrix
    end
end


# Define the function f()

f(v) = sum(v.^2)  # For example, f() could be the sum of the elements of the vector

# Create a sample simplex matrix
simplex = Float64.([1 3 5; 2 4 6])  # 2x3 matrix, with vectors [1, 2], [3, 4], [5, 6] as columns
println("originalSimplex:")
display(simplex)
v3 = [9, 12]
# Define a new vector to insert
new_vector = [2.5, 3.5]  # A new vector

# Use the function to insert the new vector
# insertSortedSimplex!(simplex, new_vector, f)
v2 = [0.5, 2.3]
# Display the updated simplex
# insertSortedSimplex!(simplex, v2, f)
# insertSortedSimplex!(simplex, v3, f)
simplex, action = nelderMead(simplex, f)
println("Action performed: $(action)")
println("Latest Simplex:")
display(simplex)



# xc = [1, 2, 3]
# xw = [4, 5, 6]
# xr = reflect(xc, xw, alpha=0.5)