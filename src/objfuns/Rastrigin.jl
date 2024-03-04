include("objective.jl")

function Rastrigin(x::Vector{Float64}, p; getGradientToo::Bool=false)
    A = 10
    n = length(x)
    sum_term = sum(x[i]^2 - A*cos(2*pi*x[i]) for i in 1:n)
    f = A*n + sum_term + p[1]

    if getGradientToo
        g = [2*x[i] + 2*pi*A*sin(2*pi*x[i]) for i in 1:n]
        return f, g
    else
        return f
    end
end
