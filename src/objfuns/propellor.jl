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
    ratio = rand()
    x0 = lb*ratio + (1-ratio)*ub
    xk = x0
    n = length(x0)
end

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

function computePolynomialEstimate(X, pw, a)
    R = length(a)
    M = size(X, 1)
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

function findOptimalPolynomialDegrees(M, X, T, Q)
    minDiscT = Inf
    NT_best = 1
    aT_best = []
    minDiscQ = Inf
    NQ_best = 1
    aQ_best = []

    for (NT, NQ) in zip(1:6, 1:6)
        discTTotal, discQTotal = 0.0, 0.0
        aT = findFitCoeffs(NT, X, T)
        aQ = findFitCoeffs(NQ, X, Q)

        discT = sum((computePolynomialEstimate(NT, X, aT) - T).^2)
        discQ = sum((computePolynomialEstimate(NQ, X, aQ) - Q).^2)
        discTTotal += discT
        discQTotal += discQ

        avgDiscT = discTTotal / M
        avgDiscQ = discQTotal / M

        # Update best N values if the current average discrepancy is lower
        if avgDiscT < minDiscT
            minDiscT = avgDiscT
            NT_best = NT
            aT_best = aT
        end
        if avgDiscQ < minDiscQ
            minDiscQ = avgDiscQ
            NQ_best = NQ
            aQ_best = aQ
        end
    end

    myprintln(true, "Best polynomial fits for Thrust T and Torque Q have been found to be $(NT_best) and $(NQ_best) respectively.")

    return NT_best, NQ_best, aT_best, aQ_best
end

NT, NQ, aT, aQ = findOptimalPolynomialDegrees(M, X, T, Q)
pwT, pwQ = constructEx.([NT, NQ])

function propellorObj(x, p;
    getGradientToo::Bool=true)

    if length(x) != 3
        @error "propellorObj expects a length 3 vector"
    end

    @unpack p = n, lb, ub, R, aQ, pw
    Q = 0

    if n != 3
        @error "Discrepancy between x and pDict[:n]"
    end

    # R = length(aQ)
    for r = 1:R
        Q += aQ[r] * prod(x .^ pw[r, :])
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
