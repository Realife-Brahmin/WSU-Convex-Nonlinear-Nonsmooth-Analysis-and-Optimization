include("objective.jl")

"""
Generalized n-dim rosenbrock function with steepness parameter p[1]
"""
function rosenbrock(x::Vector{Float64}, 
    pDict; 
    getGradientToo::Bool=true, 
    verbose::Bool=false)

    # p = p[:params]
    @unpack p = pDict
    scale = p[1]
    n = length(x)
    f = 0.0
    g = zeros(Float64, n)

    for k = 1:n-1
        T = x[k+1] - x[k]^2
        S = 1.0 - x[k]
        f += scale * T^2 + S^2
        if getGradientToo
            g[k] = -4 * scale * x[k] * T - 2 * S
            if k > 1
                g[k] += 2 * scale * (x[k] - x[k-1]^2)
            end
        end
    end
    
    if getGradientToo
        g[n] = 2 * scale * (x[n] - x[n-1]^2)
        return f, g
    else
        return f
    end
end


objective = rosenbrock;
n = 10
p = Float64.([n])
pDict = Dict(:p => p)
x0 = collect(1.0/p[1]:1.0/p[1]:1.0)
pr = generate_pr(objective, x0, params=pDict)
