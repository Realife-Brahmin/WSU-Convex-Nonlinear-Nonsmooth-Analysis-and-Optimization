include("objective.jl")

function TestFunction3(x::Vector{Float64}, p;
    getGradientToo::Bool=true)
    
    p = p.params
    beta = p[1]
    n = length(x)
    t = zeros(Float64, n)
    s = zeros(Float64, n, n)
    
    for i in 1:n
        for j in 1:n
            s[i, j] = (j + 1 + beta) * (x[j]^(i) - 1/((j + 1)^(i)))
        end
        t[i] = sum(s[i, :])
    end

    f = sum(t.^2)

    if getGradientToo
        g = zeros(Float64, n)
        for j in 1:n
            for i in 1:n
                g[j] += (t[i] * (i) * (j + 1 + beta)) * x[j]^(i - 1)
            end
        end
        return f, g
    else
        return f
    end
end

objective = TestFunction3
x0 = sort!(rand(10).^2, rev=true)
params = Float64.([10])

pr = generate_pr(objective, x0, params=params)