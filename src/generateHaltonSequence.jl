using Primes

"""
    generateHaltonSequence(q::Int, p::Int; discard::Int=0) -> Vector{Float64}

Generate a Halton sequence, a low-discrepancy sequence used in the field of numerical methods for integration and optimization. This function allows discarding a specified number of initial values in the sequence.

# Arguments
- `q::Int`: The base of the Halton sequence. A prime number is generally chosen for the base to ensure that the sequence covers the interval more uniformly.
- `p::Int`: The number of points in the Halton sequence to generate.
- `discard::Int=0` (optional): The number of initial values in the sequence to discard. Default is 0, meaning no values are discarded.

# Returns
- `Vector{Float64}`: A vector of length `p` containing the Halton sequence, after discarding the specified number of initial values.

# Details
The function generates a sequence of `p` points in the unit interval [0, 1), using the base `q`. The Halton sequence is known for its low-discrepancy properties, especially in multi-dimensional sampling contexts.

# Notes
The implementation converts the `q` to `Float64` to ensure arithmetic operations are performed in floating-point. It initializes an array `halton` of zeros with length `p`. The function then generates values for `p` indices, starting from `1 + discard` and ending at `p + discard`. For each number `k` in this range, the function calculates the Halton sequence value and stores it in the `halton` array at the corresponding index adjusted by the `discard` offset.

# Examples
```julia
# Generate a Halton sequence with base 2, 10 points, discarding the first 5 values
halton_seq = generateHaltonSequence(2, 10, discard=5)
```
"""
function generateHaltonSequence(q, p; discard::Int=0)
    q = Float64(q)
    halton = zeros(p)
    indices = 1+discard:p+discard 
    for k = indices 
        local num, map, bitNum = Float64(k), 0.0, 0.0
        while num > 0
            num, R = num ÷ q, num % q
            map += R*q^(-bitNum-1.0)
            bitNum += 1
        end
        halton[k-discard] = map
    end
    return halton
end

halt = generateHaltonSequence(2, 62, discard=22)

