FuncParam = NamedTuple{(:params, :data), Tuple{Vector{Float64}, Matrix{Float64}}}

"""
    dampedSHM(x::Vector{Float64}, p::FuncParam; verbose::Bool=false, log::Bool=true, getGradientToo::Bool=true)

Computes the mean squared deviation of the predicted damped simple harmonic motion from the given data. The motion is modeled using a combination of exponential decay and sinusoidal functions.

# Arguments
- `x::Vector{Float64}`: A vector of parameters for the motion model. They are in the order [A₀, A, τ, ω, α, ϕ].

- `p::FuncParam`: A named tuple encapsulating the parameters and data. It has two fields:
- `params::Vector{Float64}`: A vector of additional parameters (not currently used in the function).
- `data::Matrix{Float64}`: A matrix where columns correspond to different features, and rows are individual observations. The last column is the target value y, while the other columns represent the feature values (in this case, time values).

# Keyword Arguments
- `verbose::Bool`: (default `false`) If set to true, the function prints additional information during its execution.

- `log::Bool`: (default `true`) Placeholder for a logging functionality (not currently implemented in the function).

- `getGradientToo::Bool`: (default `true`) Determines whether to compute and return the gradient along with the function value.

# Returns
- If `getGradientToo` is `true`, returns a tuple `(f, g)` where `f` is the computed function value, and `g` is the gradient.

- If `getGradientToo` is `false`, returns only the function value `f`.

# Example
```julia
data_matrix = [...]
params_vector = [...]
p = FuncParam(params=params_vector, data=data_matrix)
x_values = [A₀_value, A_value, τ_value, ω_value, α_value, ϕ_value]
f, g = dampedSHM(x_values, p)
```
# Notes
The function computes the predicted damped harmonic motion using the model:
    ŷ(t) = A₀ + A * exp(-t/τ) * sin((ω+α*t) * t + ϕ)
    and then finds the mean squared deviation from the provided data.
    
    
    FuncParam = NamedTuple{(:params, :data), Tuple{Vector{Float64}, Matrix{Float64}}}
    
    function dampedSHM(x::Vector{Float64}, p::FuncParam; verbose::Bool=false, log::Bool=true, getGradientToo::Bool=true)
"""
function dampedSHM(x::Vector{Float64}, 
    p::FuncParam;
    verbose::Bool=false,
    log::Bool=true,
    getGradientToo::Bool=true)

    data = p.data
    M = size(data, 1)
    nf = size(data, 2) - 1
    y = data[:, nf+1]
    xf = data[:, 1:nf]
    t = xf

    n = length(x) # note that df's x (time) is different from function parameter x
    A₀, A, τ, ω, α, ϕ = x
    f = 0.0;
    g = zeros(Float64, n)
    for k = 1:M
        tₖ = t[k]
        yₖ = y[k]
        expₖ = exp(-tₖ/τ)
        Sₖ = expₖ*sin((ω+α*tₖ)tₖ + ϕ)
        Cₖ = expₖ*cos((ω+α*tₖ)tₖ + ϕ)

        ŷₖ = A₀ + A*Sₖ
        Δyₖ = ŷₖ - yₖ
        f += (1/M)*Δyₖ^2
        if getGradientToo
            g += (2/M)* Δyₖ * [1, Sₖ, A*tₖ*(τ^-2)*Sₖ, A*tₖ*Cₖ, A*tₖ^2*Cₖ, A*Cₖ]
        end
    end
    
    if getGradientToo
        return f, g
    else
        return f
    end
end