FuncParam = NamedTuple{(:params, :data), Tuple{Vector{Float64}, Matrix{Float64}}}

function TestFunction1(x::Vector{Float64}, 
    p::FuncParam; 
    getGradientToo::Bool=true)

    p = p.params

    a = 4 - 2.1 * x[1]^2 + (1/3) * x[1]^4
    c = 4 * x[2]^2 - 4
    f = a * x[1]^2 + x[1] * x[2] + c * x[2]^2 + p[1]

    if getGradientToo
        g = zeros(Float64, 2)
        g[1] = 2*x[1]*a + x[1]^2 * ((4/3) * x[1]^3 - 4.2*x[1]) + x[2]
        g[2] = x[1] + 2*x[2]*c + 8*x[2]^3
        return f, g
    else
        return f
    end
end


function TestFunction2(x::Vector{Float64}, p::FuncParam; getGradientToo::Bool=true)
    
    p = p.params
    a, b, c = p
    den = 2 * length(x)
    f = sum(a* x.^4 + b* x.^2 + c*x) / den + 40

    if getGradientToo
        g = (4a* x.^3 + 2b* x + c*ones(length(x))) / den
        return f, g
    else
        return f
    end
end


function TestFunction3(x::Vector{Float64}, p::FuncParam;
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

"""
Generalized n-dim rosenbrock function with steepness parameter p[1]
"""
function rosenbrock(x::Vector{Float64}, 
    p::FuncParam; 
    getGradientToo::Bool=true, 
    verbose::Bool=false)

    p = p.params
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

