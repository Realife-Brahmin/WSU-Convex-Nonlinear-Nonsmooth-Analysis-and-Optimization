include("objective.jl")
include("../AugmentedLagrangian.jl")

using CSV
using DataFrames
using Parameters

begin
    rawDataFolder = "rawData/"
    filename = rawDataFolder * "PropData.csv"
    df0 = CSV.File(filename, header=false) |> DataFrame
    rename!(df0, [:D, :alpha, :n, :T, :Q])
    # finding the maximium of the unnormalized variables, against which the variables as well as the lower and upper bounds will be normalized
    maxVals0 = maximum.([df0.D, df0.alpha, df0.n])

    df = deepcopy(df0)
    df.D /= maxVals0[1] # Dmax0 
    df.alpha /= maxVals0[2] # alphamax0
    df.n /= maxVals0[3] # nmax0
    data = Matrix(df)

    M = size(data, 1)

    X = data[:, 1:3]
    T = data[:, 4]
    Q = data[:, 5]

    lb0 = [1, 0, 1]
    ub0 = [4, 10, 7]
    lb = lb0 ./ maxVals0
    ub = ub0 ./ maxVals0

    pointNum = rand(1:M)
    # ratio = rand()
    # x0 = lb*ratio + (1-ratio)*ub
    x0 = X[pointNum, :]
    xk = x0
    n = length(x0)
    T0 = T[pointNum]
end

"""
    constructEx(N::Int) -> Vector{Vector{Int}}

Construct an array of power tuples for all terms in a polynomial in three variables (x, y, z) up to degree `N`.

# Arguments
- `N::Int`: The maximum total degree of the polynomial.

# Returns
- `Vector{Vector{Int}}`: An array of vectors, where each vector `[a, b, c]` represents a term `x^a * y^b * z^c` in the polynomial.

# Notes
This function generates all possible power combinations for a polynomial of three variables, `x`, `y`, and `z`. It iterates through all combinations of powers that sum up to values from `0` to `N`. Each combination `[a, b, c]` represents the powers of `x`, `y`, and `z` respectively, making up a term in the polynomial. The function is specifically tailored for polynomials in exactly three variables.

# Example
```julia
N = 2
powers = constructEx(N)
# Output will include combinations like [0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [1, 1, 0], etc., up to the degree N.
```
"""
function constructEx(N)
    pw = Vector{Int}[]
    for j = 0:N
        for k = 0:j
            for l = 0:k
                # This generates the tuple (j-k, k-l, l) which corresponds to the powers of x, y, z
                exElem = [(j - k), (k - l), l]
                push!(pw, exElem)
            end
        end
    end
    return pw
end

"""
    findFitCoeffs(N::Int, X::Matrix{Float64}, T::Vector{Float64}) -> (Vector{Float64}, Vector{Vector{Int}})

Compute the coefficients for a polynomial fit of degree `N` for a given set of data points and corresponding output values.

# Arguments
- `N::Int`: The degree of the polynomial.
- `X::Matrix{Float64}`: A matrix where each row represents a data point with three variables (e.g., Degree, alpha, and n_RPM for a propeller).
- `T::Vector{Float64}`: A vector containing the output values (e.g., Thrust or Torque) corresponding to each row in `X`.

# Returns
- `(Vector{Float64}, Vector{Vector{Int}})`: A tuple where the first element is a vector of coefficients for the polynomial terms, and the second element is a vector of vectors, each representing the powers of the polynomial terms.

# Notes
This function generates a matrix of coefficients (`AT`) and a vector (`BT`) for solving the linear system `AT * a = BT` to find the polynomial coefficients. The polynomial terms are generated for three variables up to the specified degree `N`. The function is tailored for problems involving exactly three input variables.

# Example
```julia
N = 2
X = [1.0 2.0 3.0; 4.0 5.0 6.0; 7.0 8.0 9.0]
T = [1, 2, 3]
a, pw = findFitCoeffs(N, X, T)
# `a` will contain the coefficients of the polynomial fit, and `pw` will show the corresponding powers of each term.
```
"""
function findFitCoeffs(N, X, T)
    pw = constructEx(N)
    R = length(pw)
    AT = zeros(R, R)
    BT = zeros(R)
    for r = 1:R
        Tk = T[r]
        x, y, z = X[r, :]
        temp = x.^pw[r][1] .* y.^pw[r][2] .* z.^pw[r][3]
        BT[r] = sum(Tk .* temp)
        for c = 1:R
            p = pw[r] + pw[c]
            AT[r, c] = sum( x.^p[1] .* y.^p[2] .* z.^p[3] )
        end
    end
    a = AT\BT
    return a, pw
end

"""
    computePolynomialEstimate(X, pw, a)

Estimates the values of a polynomial function based on provided coefficients over a set of input data points.

# Arguments
- `X::Union{AbstractVector, AbstractMatrix}`: Input data points. Each row represents a single set of variables (x, y, z).
- `pw::Vector{Vector{Int}}`: Powers of each variable in the polynomial terms.
- `a::Vector`: Coefficients for each polynomial term.

# Returns
- `T_est::Vector{Float64}`: Estimated values of the polynomial at each data point in `X`.

# Details
This function computes the polynomial value for each data point specified in `X`. The polynomial is defined by the terms in `pw` and their corresponding coefficients in `a`.

# Notes
- `X` can either be a matrix of multiple data points or a single vector for one data point. If it's a vector, it is reshaped internally for consistent processing.
- If `X` is provided as a single vector, it is treated as a single data point.

# Examples
```julia
X = [1 2 3; 4 5 6]  # Two points in R^3
pw = [[1, 0, 0], [0, 1, 1]]  # Corresponds to x and y*z
a = [2, 3]  # Coefficients for x and y*z

# Compute polynomial estimates
T_est = computePolynomialEstimate(X, pw, a)
```
"""
function computePolynomialEstimate(X, pw, a)
    R = length(a)
    M = size(X, 1)
    if ndims(X) == 1  # X is a single vector treated as one data point
        M = 1  # There is just one data point
        X = reshape(X, 1, :)  # Reshape X to be a 1xN matrix for consistency
    end
    T_est = zeros(M)
    for i = 1:M
        x, y, z = X[i, :]
        for r = 1:R
            temp = x .^ pw[r][1] .* y .^ pw[r][2] .* z .^ pw[r][3]
            T_est[i] += a[r] .* temp
        end
    end
    return T_est
end

"""
    computePolynomialEstimateTermwise(Xi, pw, a)

Calculates the individual term contributions of a polynomial for a single data point.

# Arguments
- `Xi::AbstractVector`: A single data point in R^3, typically represented as [x, y, z].
- `pw::Vector{Vector{Int}}`: Powers for each variable in the polynomial terms.
- `a::Vector`: Coefficients corresponding to each term in the polynomial.

# Returns
- `T_est_1toR::Vector{Float64}`: The computed values of each polynomial term for the given data point `Xi`.

# Details
This function is used to compute the contribution of each term of the polynomial individually for a single data point, which is useful for debugging or understanding the influence of each term.

# Notes
- The input `Xi` must always be a vector of length 3, representing a single point in R^3. This function is specifically designed for term-by-term evaluation where each term's contribution is isolated.
- This function can be used to analyze or visualize how each polynomial term affects the overall polynomial evaluation, which can be useful in polynomial regression diagnostics.

# Examples
```julia
Xi = [2, 3, 4]  # A single data point
pw = [[2, 0, 0], [0, 1, 1]]  # Powers for x^2 and y*z
a = [1, -1]  # Coefficients for x^2 and y*z

# Calculate individual term contributions
T_est_terms = computePolynomialEstimateTermwise(Xi, pw, a)
```
"""
function computePolynomialEstimateTermwise(Xi, pw, a)

    R = length(a)

    T_est_1toR = zeros(R)
    x, y, z = Xi
    for r = 1:R
        temp = x .^ pw[r][1] .* y .^ pw[r][2] .* z .^ pw[r][3]
        T_est_1toR[r] += a[r] .* temp
    end

    return T_est_1toR
end

"""
    findOptimalPolynomialDegrees(M::Int, X::Matrix{Float64}, T::Vector{Float64}, Q::Vector{Float64}) -> Tuple

Finds the optimal polynomial degrees for fitting polynomial models to two sets of data (Thrust T and Torque Q) and returns the best degrees and their corresponding coefficients and power vectors.

# Arguments
- `M::Int`: Number of data points.
- `X::Matrix{Float64}`: A matrix of input data where each row represents a data point and columns represent variables (e.g., Degree, alpha, n_RPM for a propeller).
- `T::Vector{Float64}`: A vector of target values for Thrust.
- `Q::Vector{Float64}`: A vector of target values for Torque.

# Returns
- `NT_best::Int`: Optimal polynomial degree for Thrust.
- `NQ_best::Int`: Optimal polynomial degree for Torque.
- `aT_best::Vector{Float64}`: Coefficients of the best fitting polynomial for Thrust.
- `aQ_best::Vector{Float64}`: Coefficients of the best fitting polynomial for Torque.
- `pwT_best::Vector{Vector{Int}}`: Powers of each term in the best fitting polynomial for Thrust.
- `pwQ_best::Vector{Vector{Int}}`: Powers of each term in the best fitting polynomial for Torque.

# Notes
The function iterates over polynomial degrees from 1 to 6 for both Thrust and Torque, computes the total discrepancy between the actual data and the polynomial estimates, and identifies the degree which results in the minimum discrepancy. The discrepancy is calculated as the sum of squared differences between the estimated and actual values. This process is conducted independently for both Thrust and Torque.

# Example
```julia
M = 100
X = rand(M, 3)  # Example input matrix with 3 variables
T = rand(M)     # Random Thrust values
Q = rand(M)     # Random Torque values
NT_best, NQ_best, aT_best, aQ_best, pwT_best, pwQ_best = findOptimalPolynomialDegrees(M, X, T, Q)
# Now `NT_best` and `NQ_best` hold the optimal polynomial degrees for Thrust and Torque, respectively,
# and `aT_best`, `aQ_best`, `pwT_best`, `pwQ_best` contain the coefficients and power vectors for the best fits.
```
"""
function findOptimalPolynomialDegrees(M, X, T, Q)
    minDiscT = Inf
    NT_best = 1
    aT_best = []
    pwT_best = []
    minDiscQ = Inf
    NQ_best = 1
    aQ_best = []
    pwQ_best = []

    for (NT, NQ) in zip(1:6, 1:6)
        discTTotal, discQTotal = 0.0, 0.0
        aT, pwT = findFitCoeffs(NT, X, T)
        aQ, pwQ = findFitCoeffs(NQ, X, Q)

        discT = sum((computePolynomialEstimate(X, pwT, aT) - T).^2)
        discQ = sum((computePolynomialEstimate(X, pwQ, aQ) - Q).^2)
        discTTotal += discT
        discQTotal += discQ

        avgDiscT = discTTotal / M
        avgDiscQ = discQTotal / M

        # Update best N values if the current average discrepancy is lower
        if avgDiscT < minDiscT
            minDiscT = avgDiscT
            NT_best = NT
            aT_best = aT
            pwT_best = pwT
        end
        if avgDiscQ < minDiscQ
            minDiscQ = avgDiscQ
            NQ_best = NQ
            aQ_best = aQ
            pwQ_best = pwQ
        end
    end

    myprintln(true, "Best polynomial fits for Thrust T and Torque Q have been found to be $(NT_best) and $(NQ_best) respectively.")

    return NT_best, NQ_best, aT_best, aQ_best, pwT_best, pwQ_best
end

NT, NQ, aT, aQ, pwT, pwQ = findOptimalPolynomialDegrees(M, X, T, Q)

function propellorObj(x, p;
    getGradientToo::Bool=true)

    if length(x) != 3
        @error "propellorObj expects a length 3 vector"
    end

    @unpack p = n, lb, ub, R, aQ, pwQ
    Q = 0

    if n != 3
        @error "Discrepancy between x and pDict[:n]"
    end

    # R = length(aQ)
    for r = 1:R
        Q += aQ[r] * prod(x .^ pwQ[r, :])
    end

    f = x[1] * x[3] * Q

    if !getGradientToo
        return f
    elseif getGradientToo
        # compute g yourself
        return f, g # do
    else
        @error("floc")
    end

    @error("floc")
end

mE = 1

function propellorEcons(x, p;
    getGradientToo::Bool=true)

    @unpack n, mE, aT, NT, T0 = p

    if length(x) != 3 || n != 3
        @error "propellorEcons expects a length 3 vector"
    end

    mE = 1
    cE = zeros(mE)
    cE[1] = computePolynomialEstimate(NT, x, aT) - T0
    
    if !getGradientToo
        return cE
    elseif getGradientToo
        hE = zeros(mE, n)
        # do : compute hE
        return cE, hE
    else
        @error("floc")
    end

    @error("floc")

end

mI = 6
y0 = rand(mI)
yk = y0

function propellorIcons(x, p;
    getGradientToo::Bool=true)

    @unpack n, mI, lb, ub, = p

    if length(x) != 3 || n != 3
        @error "propellorIcons expects a length 3 vector"
    end

    mI = 6
    cI = zeros(mI)
    cI[1] = ub[1] - x[1]
    cI[2] = ub[2] - x[2]
    cI[3] = ub[3] - x[3]
    cI[4] = x[1] - lb[1]
    cI[5] = x[2] - lb[2]
    cI[6] = x[3] - lb[3]

    if !getGradientToo
        return cI
    elseif getGradientToo
        n = length(x)
        hI = zeros(mI, n)
        hI[1:3, :] = -I(3)
        hI[4:6, :] = I(3)
        return cI, hI
    else
        @error("floc")
    end

    @error("floc")

end

m = mE+mI
w0 = vcat(x0, y0)
wk = w0

objective = propellorObj
objectiveOriginal = propellorObj
objectiveString = "propellorObj"
problemType = "Constrained"
econ = propellorEcons
icon = propellorIcons
T0 = T[20]

pDictALP = Dict(:n=>n, :m=>m, :mE=>mE, :econ=>econ, :mI=>mI, :icon=>icon, :NT=>NT, :pwT => pwT, :aT=>aT, :NQ=>NQ, :pwQ=>:pwQ, :aQ=>aQ, :lb=>lb, :ub=>ub, :T0=>T0)

pr = generate_pr(objective, w0, params=pDictALP, problemType=problemType; objectiveString=objectiveString);
