function TestFunction1(x::Vector{Float64}, p::Float64; getGradientToo::Bool=true)

    a = 4 - 2.1 * x[1]^2 + (1/3) * x[1]^4
    c = 4 * x[2]^2 - 4
    f = a * x[1]^2 + x[1] * x[2] + c * x[2]^2 + p

    if getGradientToo
        g = zeros(Float64, 2)
        g[1] = 2*x[1]*a + x[1]^2 * ((4/3) * x[1]^3 - 4.2*x[1]) + x[2]
        g[2] = x[1] + 2*x[2]*c + 8*x[2]^3
        return f, g
    else
        return f
    end
end


function TestFunction2(x::Vector{Float64}, p::Vector{Float64}; getGradientToo::Bool=true)
    
    a, b, c = p
    den = 2 * length(x)
    f = sum(a .* x.^4 + b .* x.^2 + c .* x) / den + 40

    if getGradientToo
        g = (4*a .* x.^3 + 2*b .* x + c) / den
        return f, g
    else
        return f
    end
end


function TestFunction3(x::Vector{Float64}, p::Float64; getGradientToo::Bool=true)
    
    beta = p
    n = length(x)
    t = zeros(Float64, n)
    s = zeros(Float64, n, n)
    
    for i in 1:n
        for j in 1:n
            s[i, j] = (j + 1 + beta) * (x[j]^(i) - (j + 1)^(-i))
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
