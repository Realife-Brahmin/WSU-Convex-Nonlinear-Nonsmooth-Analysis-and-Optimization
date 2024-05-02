include("objective.jl")
include("../AugmentedLagrangian.jl")

using CSV
using DataFrames
using Parameters

rawDataFolder = "rawData/"
filename = rawDataFolder * "PropData.csv"
df = CSV.File(filename, header=false) |> DataFrame
rename!(df, [:D, :alpha, :n, :T, :Q])
df.D /= maximum(df.D)
df.alpha /= maximum(df.alpha)
df.n /= maximum(df.n)

data = Matrix(df)
x0 = rand(2)
xk = x0
n = length(x0)

function propellorObj(x, p;
    getGradientToo::Bool=true)

    if length(x) != 2
        @error "alpTestFunction01 expects a length 2 vector"
    end

    Q = 0
    # b = p[:b]
    m = length(p[:b])
    for r = 1:m
        Q += p[:b][r] * prod(x .^ p[:pw][r, :])
    end
    f = x[1]*x[3]*Q

    if !getGradientToo
        return f
    elseif getGradientToo
        # compute g yourself
        return f, g
    else
        @error("floc")
    end

    @error("floc")
end

M = size(df, 1)

function constructEx(N)
    ex = Vector{Int}[]
    for j = 0:N
        for k = 0:j
            for l = 0:k
                # This generates the tuple (j-k, k-l, l) which corresponds to the powers of x, y, z
                exElem = [(j - k), (k - l), l]
                push!(ex, exElem)
            end
        end
    end
    return ex
end

X = data[:, 1:3]
T = data[:, 4]
Q = data[:, 5]

function findFitCoeffs(N, M, X, T)
    ex = constructEx(N)
    R = length(ex)
    AT = zeros(R, R)
    BT = zeros(R)
    for r = 1:R
        Tk = T[r]
        x, y, z = X[r, :]
        temp = x.^ex[r][1] .* y.^ex[r][2] .* z.^ex[r][3]
        BT[r] = sum(Tk .* temp)
        for c = 1:R
            p = ex[r] + ex[c]
            AT[r, c] = sum( x.^p[1] .* y.^p[2] .* z.^p[3] )
        end
    end
    a = AT\BT
    return a
end

function computePoynomialEstimate(N, X, a)
    R = length(a)
    ex = constructEx(N)
    estimate = 0.0
    x, y, z = X
    for r = 1:R
        temp = x .^ ex[r][1] .* y .^ ex[r][2] .* z .^ ex[r][3]
        estimate += a[r] * sum(temp)
    end
    return estimate
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
        aT = findFitCoeffs(NT, M, X, T)
        aQ = findFitCoeffs(NQ, M, X, Q)

        for i = 1:M
            Xi, Ti, Qi = X[i, :], T[i], Q[i]
            discT = (computePoynomialEstimate(NT, Xi, aT) - Ti)^2
            discQ = (computePoynomialEstimate(NQ, Xi, aQ) - Qi)^2
            discTTotal += discT
            discQTotal += discQ
        end

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

    return NT_best, NQ_best, aT_best, aQ_best
end

NT, NQ, aT, aQ = findOptimalPolynomialDegrees(M, X, T, Q)


# function cE01(x, p;
#     getGradientToo::Bool=true)

#     if length(x) != 2
#         @error "cE01 expects a length 2 vector"
#     end

#     mE = 1
#     cE = zeros(mE)
#     cE[1] = x[1]^2 + x[2]^2 - 1

#     if !getGradientToo
#         return cE
#     elseif getGradientToo
#         n = length(x)
#         hE = zeros(mE, n)
#         hE[1, 1] = 2*x[1]
#         hE[1, 2] = 2*x[2]
#         return cE, hE
#     else
#         @error("floc")
#     end

#     @error("floc")

# end

mE = 0

function propellorIcons(x, p;
    getGradientToo::Bool=true)

    if length(x) != 2
        @error "propellorIcons expects a length 2 vector"
    end

    mI = 2
    cI = zeros(mI)
    cI[1] = x[1]^3 + x[2]
    cI[2] = x[1]^2 + 2*x[2]^2 - 1.1
    if !getGradientToo
        return cI
    elseif getGradientToo
        n = length(x)
        hI = zeros(mI, n)
        hI[1, 1] = 3*x[1]^2
        hI[1, 2] = 1
        hI[2, 1] = 2*x[1]
        hI[2, 2] = 4*x[2]
        return cI, hI
    else
        @error("floc")
    end

    @error("floc")

end

mI = 2
y0 = rand(mI)
yk = y0
m = mE+mI
w0 = vcat(x0, y0)
wk = w0

objective = propellorObj
objectiveOriginal = propellorObj
objectiveString = "propellorObj"
problemType = "Constrained"
econ = nothing
icon = propellorIcons
;
# pDictALP = Dict(:n=>n, :m=>m, :mE=>mE, :econ=>econ, :mI=>mI, :icon=>icon)
# pr = generate_pr(objective, w0, params=pDictALP, problemType=problemType; objectiveString=objectiveString)

# objectiveUnc = ALOBJ
# subroutineCall = true
# muk = 1.0
# m = mE+mI
# lambdak = ones(m)
# addendum = Dict(:subroutineCall => subroutineCall, :lambda => lambdak, :mu => muk, :objective => objective, :objectiveUnc => objectiveUnc)
# pDictUnc = merge(deepcopy(pDictALP), addendum)
# # xk = x0
# F, G = ALOBJ(wk, pDictUnc, getGradientToo=true)

# using JuMP, Ipopt

# model = Model(Ipopt.Optimizer)
# # x0 = [-1.2, 1.0]
# # x0 = [1, 2]
# slackifyInequalities = false
# slackifyInequalities = true
# # n = length(x0)
# @variable(model, x[i=1:n], start = x0[i])
# @variable(model, y[i=1:mI], start = y0[i])
# @objective(model, Min, (x[1] - 1)^2 + (x[2] - 2)^2)
# @constraint(model, x[1]^2 + x[2]^2 - 1 == 0) # regular equality constraint

# if slackifyInequalities
#     @constraint(model, x[1]^3 + x[2] - y[1]^2 == 0) 
#     @constraint(model, x[1]^2 + 2 * x[2]^2 - 1.1 - y[2]^2 == 0)
# else
#     @constraint(model, x[1]^3 + x[2] >= 0)
#     @constraint(model, x[1]^2 + 2 * x[2]^2 - 1.1 >= 0)
# end

# optimize!(model)
# xopt = [value(x[i]) for i ∈ 1:n]
# yopt = [value(y[i]) for i ∈ 1:mI]
# f_optimal = objective_value(model)
# println("Optimal Variables x: ", xopt)
# println("Optimal Slack Variables y: ", yopt)
# println("Optimal objective value: ", f_optimal)

