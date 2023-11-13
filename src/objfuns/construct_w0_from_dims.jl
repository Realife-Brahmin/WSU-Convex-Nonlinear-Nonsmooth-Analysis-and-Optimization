function construct_w0_from_dims(dims::Vector{Int})
    numWeights = getNumWeights(dims)
    w0 = ones(numWeights)
    return w0
end

function getNumWeights(dims::Vector{Int})
    # Initialize the sum
    numWeights = 0
    
    # Loop through the dimensions vector
    for i in 1:length(dims)-1
        # Multiply each pair of adjacent dimensions and add to the sum
        numWeights += dims[i] * dims[i+1]
    end

    return numWeights
end

# # Example usage:
# dim = [10, 10, 10, 1]
# result = construct_w0_from_dims(dim);
# println(result)  # Output should be 210
