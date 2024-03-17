include("objective.jl")

function Rastrigin(x::Vector{Float64}, p;           
    getGradientToo::Bool=true,
    verbose::Bool=false)

    A = 10
    n = length(x)
    # sum_term = sum(x[i]^2 - A*cos(2*pi*x[i]) for i in 1:n)
    sum_term = sum(x.^2 - A*cos.(2*π*x))

    f = A*n + sum_term

    if getGradientToo
        # g = [2*x[i] + 2*pi*A*sin(2*pi*x[i]) for i in 1:n]
        g = 2x + 2π*A*sin.(2π*x)
        return f, g
    else
        return f
    end

end


x0 = rand(7)
# x = x0
# p = Dict()
# f, g = Rastrigin(x, p)

params = Dict();

objective = Rastrigin;

pr = generate_pr(objective, x0, params=params);