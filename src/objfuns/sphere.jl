include("objective.jl")

function sphere(x::Vector{Float64}, p;
    getGradientToo::Bool=true,
    verbose::Bool=false)

    f = dot(x, x)
    if getGradientToo
        g = 2x
        return f, g
    else
        return f
    end
end

objective = sphere;
# x0 = collect(0.1:0.1:1)
x0 = Float64.([1, 2])

pr = generate_pr(objective, x0)