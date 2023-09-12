# module objective

using DataFrames
using Symbolics

# export objFun

function objFun(df::DataFrame;
    getGradientToo::Bool=true, 
    getFunctionsToo::Bool=true,
    getVariables::Bool=true, 
    verbose::Bool=false)

    N = size(df, 1)
    @variables t, A₀, A, τ, ω, α, ϕ
    x = [A₀, A, τ, ω, α, ϕ] 
    # @show typeof(x)
    # @show eltype(x)
    fsym = A₀ + A*exp(-t/τ)sin((ω+α*t)t + ϕ)
    f = sum([(df.V[i] - substitute(fsym, t => df.t[i]))^2 for i in 1:N])

    if getGradientToo
        # Calculate the gradient for each squared difference
        gradients = [Symbolics.gradient(f[i], x) 
        for i ∈ eachindex(f)]   
        # Sum up the gradients to get the overall gradient
        ∇f = sum(gradients)
        if getFunctionsToo
            fnum = build_function(f, x, expression=Val{false})
            ∇fnum = build_function(∇f, x, expression=Val{false})
            return f, ∇f, fnum, ∇fnum, x
        else
            return f, ∇f, x
        end
    else
        if getFunctionsToo
            fnum = build_function(f, x, expression=Val{false})
            return f, fnum, x
        else
            return f, x
        end
    end
end

# end