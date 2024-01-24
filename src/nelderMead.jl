using Statistics

include("helperFunctions.jl")
include("sampleSpace.jl")

function NelderMead(simplex, f::Function;
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
            # simplex = insertSortedSimplex() 
        else
            F_xw = f(xw)
            if F_xr < F_xw
                # this point is technically better than the worst point in the simplex, but not particularly useful. Can perform some contracts on it.
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

function insertSortedSimplex!(matrix::Matrix{T}, new_vector::Vector{T}, f::Function) where T
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
f(v) = sum(v)  # For example, f() could be the sum of the elements of the vector

# Create a sample simplex matrix
simplex = Float64.([1 3 5; 2 4 6])  # 2x3 matrix, with vectors [1, 2], [3, 4], [5, 6] as columns

# Define a new vector to insert
new_vector = [2.5, 3.5]  # A new vector

# Use the function to insert the new vector
insertSortedSimplex!(simplex, new_vector, f)
v2 = [0.5, 2.3]
# Display the updated simplex
insertSortedSimplex!(simplex, v2, f)
println(updated_simplex)



xc = [1, 2, 3]
xw = [4, 5, 6]
# xr = reflect(xc, xw, alpha=0.5)