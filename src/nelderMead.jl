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
        myprintln(verbose, "Extending good reflection.")
        xe = extend(xc, xw, gamma=gamma)
        F_xe = f(xe)
        if F_xe < F_xr
            myprintln(verbose, "Extension a success!")
            # extended point is even better! 
            # add it to the top of the simplex
            simplex = hcat(xe, simplex[:, 2:p-1])
            action = "extend"
        else 
            # extended point is a worse point than reflection
            # add reflected point to the top of the list
            myprintln(verbose, "Extension a failure. Reverting to reflection.")
            simplex = hcat(xr, simplex[:, 1:p-1])
        end
    else
        F_xsw = f(xsw)
        if F_xr < F_xsw 
            # point satisfies minimum required criterion for addition to simplex
            # just select it and move on
            myprintln(verbose, "Reflection better than second-worst point, so adding it to simplex")
            simplex = insertSortedSimplex(simplex, xr, f)
        else
            F_xw = f(xw)
            if F_xr < F_xw
                # this point is technically better than the worst point in the simplex, but not particularly useful. Can perform some contracts on it.
                myprintln(verbose, "Reflection only better than worst point. Trying outside contract.")
                @show xoc = outsideContract(xc, xr, beta=beta)
                @show length(xoc)
                F_xoc = f(xoc)
                if F_xoc < F_xr
                    myprintln(verbose, "Outside Contract a success! Adding it to simplex.")
                    # choosing outside contracted point (it may or may not be a useful addition to the simplex)
                    action = "outsideContract"
                    simplex = insertSortedSimplex(simplex, xoc, f)
                    display(simplex)
                else
                    # outside contract didn't help, but choosing reflected point anyway
                    myprintln(verbose, "Outside Contract a failure. Still adding the reflection to simplex.")
                    simplex = hcat(simplex[1:p-1], xr)
                end
            else
                # oh no the point is even worse than the current worst point
                # Let's try an inside contract
                myprintln(verbose, "Worse than the worst point. Trying an Inside Contract.")
                xic = insideContract(xc, xw, beta=beta)
                action = "insideContract"
                F_xic = f(xic)
                if F_xic < F_xw
                    myprintln(verbose, "Inside Contract better than worst point, so adding it.")
                    simplex = insertSortedSimplex(simplex, xic, f)
                else
                    # okay no improvment this time
                    # let's shrink the simplex
                    myprintln(verobse, "No improvement. Shrink the simplex.")
                    simplex = shrinkSortedSimplex(simplex, delta=delta)
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

function shrinkSortedSimplex(simplex;
    delta = 0.5)
    simplex += delta*(simplex[1, :] .- simplex[:, :])
    return simplex
end

function insertSortedSimplex(matrix, new_vector, f::Function)
    values = [f(matrix[:, i]) for i in 1:size(matrix, 2)]
    new_value = f(new_vector)

    # Find the insertion index
    insert_index = searchsortedfirst(values, new_value)

    # Insert the new vector and then trim the last column if necessary
    matrix .= hcat(matrix[:, 1:insert_index-1], new_vector, matrix[:, insert_index:end-1])

    return matrix
    
end

function f(v)
    x, y = v[1], v[2]
    f = (x-1)*x + (y+1)*y
    return f
end

simplex = Float64.([0 1 3;0 2 4])
# Create a sample simplex matrix
# simplex = Float64.([1 3 5; 2 4 6])  # 2x3 matrix, with vectors [1, 2], [3, 4], [5, 6] as columns
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
verbose = true
# verbose = false
simplex = insertSortedSimplex(simplex, [-0.75, -0.5], f)
println("Latest Simplex:")
display(simplex)
# simplex, action = nelderMead(simplex, f, verbose=verbose)
# println("Action performed: $(action)")
# println("Latest Simplex:")
# display(simplex)