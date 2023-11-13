function construct_w0_from_dims(dims::Vector{Int})
    # Initialize the sum
    numWeights = 0
    
    # Loop through the dimensions vector
    for i in 1:length(dims)-1
        # Multiply each pair of adjacent dimensions and add to the sum
        numWeights += dim[i] * dim[i+1]
    end

    w0 = ones(numWeights)
    return w0
end

# # Example usage:
# dim = [10, 10, 10, 1]
# result = construct_w0_from_dims(dim);
# println(result)  # Output should be 210
