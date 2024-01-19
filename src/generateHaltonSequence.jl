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

halt = generateHaltonSequence(2, 62, discard=22)

function sampleSpaceWithHalton(n, p;
    discard::Int = 0)

    sampledSpace = zeros(n, p)
    qs = primes(10*n)[1:n]
    for keyNum = 1:n
        local q = qs[keyNum]
        sampledSpace[keyNum, :] = generateHaltonSequence(q, p, discard=discard)
    end

    return sampledSpace
end

n = 2
p = 62
discard = 20
ss = sampleSpaceWithHalton(n, p, discard=discard)
Plots.theme(:dao)
p1 = scatter(ss[1, :], ss[2, :],
        aspect_ratio=:equal,
        xlims=[0.0, 1.0],
        ylims=[0.0, 1.0],
        label=:none)
display(p1)

function sampleSpace(n, p; 
    method::String="Halton",
    discard::Int=0)

    space = zeros(n, p)
    if method == "Halton"
        sampledSpace = sampleSpaceWithHalton(n, p, discard=discard)
        qs = primes(10*n)[1:n]

    elseif method == "Latin Hypercube"
        @error "not implemented"
    elseif method == "Random"
        @error "not implemented"
    else
        @error "floc"
    end

end
