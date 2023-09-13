# module objective

using DataFrames
using Symbolics

# export objFun
"""
    objFun(df::DataFrame; getGradientToo=true, getFunctionsToo=true,
        getVariables=true, verbose=false)

Compute the objective function for a decaying sinusoidal model.

# Arguments

- `df::DataFrame`: A DataFrame with columns `t` and `V` where `t` represents time in seconds and `V` represents voltage in millivolts.

# Keyword Arguments

- `getGradientToo::Bool=true`: If true, the function will also compute and return the gradient of the objective function.

- `getFunctionsToo::Bool=true`: If true, the function will also return the numerical (compiled) versions of the symbolic objective function (and its gradient if `getGradientToo` is true).

- `getVariables::Bool=true`: Not used in the provided function code.

- `verbose::Bool=false`: If true, will print out detailed information during the computation. Not implemented in the provided function code.

# Returns

- `f`: The symbolic representation of the objective function.

- `∇f` (optional): The symbolic gradient of the objective function. Only returned if `getGradientToo=true`.

- `fnum` (optional): The compiled numerical function of the objective function. Only returned if `getFunctionsToo=true`.

- `∇fnum` (optional): The compiled numerical function of the gradient. Only returned if `getGradientToo=true` and `getFunctionsToo=true`.

- `x`: The symbolic vector of parameters used in the model.

# Example

```julia
df = DataFrame(t=[0.1, 0.2, 0.3], V=[0.5, 0.4, 0.3])
f, ∇f, x = objFun(df)
"""
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
    f = mean([(df.V[i] - substitute(fsym, t => df.t[i]))^2 for i in 1:N])

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

function findDirection(pr::NamedTuple, x_now::Vector{Float64}, ∇f;
    verbose::Bool=false)::Vector{Float64}
    method = pr.method
    N = length(x_now)
    if method == "GradientDescent"
        Bₖ = I(n_now)
        ∇f_now = ∇f(x_now)
        pₖ = -Bₖ*∇f_now
    else 
        @error "Currently not formulated for this method"
    end

    return pₖ
end

function linesearch(pr::NamedTuple, x0::Vector{Float64};
    verbose::Bool=false)::Float64

end
# end