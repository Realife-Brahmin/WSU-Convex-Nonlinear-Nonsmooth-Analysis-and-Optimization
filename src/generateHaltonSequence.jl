using Primes

"""
    generateHaltonSequence(key::Int, len::Int; discard::Int=0) -> Vector{Float64}

Generate a Halton sequence, a low-discrepancy sequence used in the field of numerical methods for integration and optimization. This function allows discarding a specified number of initial values in the sequence.

# Arguments
- `key::Int`: The base of the Halton sequence. A prime number is generally chosen for the base to ensure that the sequence covers the interval more uniformly.
- `len::Int`: The number of points in the Halton sequence to generate.
- `discard::Int=0` (optional): The number of initial values in the sequence to discard. Default is 0, meaning no values are discarded.

# Returns
- `Vector{Float64}`: A vector of length `len` containing the Halton sequence, after discarding the specified number of initial values.

# Details
The function generates a sequence of `len` points in the unit interval [0, 1), using the base `key`. The Halton sequence is known for its low-discrepancy properties, especially in multi-dimensional sampling contexts.

# Notes
The implementation converts the `key` to `Float64` to ensure arithmetic operations are performed in floating-point. It initializes an array `halton` of zeros with length `len`. The function then generates values for `len` indices, starting from `1 + discard` and ending at `len + discard`. For each number `k` in this range, the function calculates the Halton sequence value and stores it in the `halton` array at the corresponding index adjusted by the `discard` offset.

# Examples
```julia
# Generate a Halton sequence with base 2, 10 points, discarding the first 5 values
halton_seq = generateHaltonSequence(2, 10, discard=5)
```
"""
function generateHaltonSequence(key, len; discard::Int=0)
    key = Float64(key)
    halton = zeros(len)
    indices = 1+discard:len+discard 
    for k = indices 
        local num, map, bitNum = Float64(k), 0.0, 0.0
        while num > 0
            num, R = num ÷ key, num % key
            map += R*key^(-bitNum-1.0)
            bitNum += 1
        end
        halton[k-discard] = map
    end
    return halton
end

# halt = generateHaltonSequence(2, 62, discard=22)

"""
    sampleSpaceHalton(n::Int, p::Int; discard::Int=0) -> Matrix{Float64}

Generate a sample space using Halton sequences with varying bases. This function is useful for quasi-random sampling in high-dimensional spaces, particularly in numerical integration and optimization problems.

# Arguments
- `n::Int`: The dimension of the sample space. Each dimension will use a different prime number as the base for its Halton sequence.
- `p::Int`: The number of points to generate in each dimension of the sample space.
- `discard::Int=0` (optional): The number of initial values in each Halton sequence to discard. Default is 0, meaning no values are discarded.

# Returns
- `Matrix{Float64}`: A `n x p` matrix where each row represents a dimension of the sample space, populated with values from the corresponding Halton sequence.

# Details
The function creates a matrix `sampledSpace` of zeros with dimensions `n x p`. It then generates `n` distinct prime numbers using the first `n` primes greater than `10 * n`. Each of these primes serves as the base for a Halton sequence for each dimension of the sample space. The `generateHaltonSequence` function is called for each dimension with the respective prime number as the base and `p` as the number of points, along with the specified `discard` value.

# Examples
```julia
# Generate a 2-dimensional sample space with 100 points in each dimension, discarding the first 5 values of each Halton sequence
sample_space = sampleSpaceHalton(2, 100, discard=5)
```
"""
function sampleSpaceHalton(n, p;
    discard::Int = 0)

    sampledSpace = zeros(n, p)
    qs = primes(10*n)[1:n]
    for keyNum = 1:n
        local q = qs[keyNum]
        sampledSpace[keyNum, :] = generateHaltonSequence(q, p, discard=discard)
    end

    return sampledSpace
end
