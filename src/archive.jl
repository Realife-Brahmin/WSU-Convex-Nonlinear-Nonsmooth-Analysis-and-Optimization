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

- `fobj`: The symbolic representation of the objective function.

- `∇fobj` (optional): The symbolic gradient of the objective function. Only returned if `getGradientToo=true`.

- `fnum` (optional): The compiled numerical function of the objective function. Only returned if `getFunctionsToo=true`.

- `∇fnum` (optional): The compiled numerical function of the gradient. Only returned if `getGradientToo=true` and `getFunctionsToo=true`.

- `x`: The symbolic vector of parameters used in the model.

# Example

```julia
df = DataFrame(t=[0.1, 0.2, 0.3], V=[0.5, 0.4, 0.3])
fobj, ∇fobj, x = objFun(df)
"""
function objFun(df::DataFrame;
    getGradientToo::Bool=true, 
    getFunctionsToo::Bool=true,
    getVariables::Bool=true, 
    verbose::Bool=false)

    N = size(df, 1)
    @variables t, A₀, A, τ, ω, α, ϕ
    x = [A₀, A, τ, ω, α, ϕ] 

    fsig = A₀ + A*exp(-t/τ)sin((ω+α*t)t + ϕ)
    fobj = mean([(df.V[i] - substitute(fsig, t => df.t[i]))^2 for i in 1:N])

    if getGradientToo
        # Calculate the gradient for each squared difference
        gradients = [Symbolics.gradient(fobj[i], x) 
        for i ∈ eachindex(fobj)]   
        # Sum up the gradients to get the overall gradient
        ∇fobj = sum(gradients)
        if getFunctionsToo
            fnum = build_function(fobj, x, expression=Val{false})
            ∇fnum = build_function(∇fobj, x, expression=Val{false})
            return fobj, ∇fobj, fnum, ∇fnum, x
        else
            return fobj, ∇fobj, x
        end
    else
        if getFunctionsToo
            fnum = build_function(fobj, x, expression=Val{false})
            return fobj, fnum, x
        else
            return fobj, x
        end
    end
end

"""
    evaluateFunction(pr, x::Vector{Float64}, t::Float64; kwargs...)

Dynamically evaluates a function specified by the `pr.objective` string, given the vector `x`, parameter `t`, and any additional keyword arguments.

# Arguments
- `pr`: An object (typically a NamedTuple) containing configurations. Specifically:
    * `pr.objective`: A string that specifies the name of the function to be evaluated.
- `x::Vector{Float64}`: The input vector for the function.
- `t::Float64`: A parameter required by the function.

# Keyword Arguments
- `kwargs...`: Additional keyword arguments that might be required by the function being evaluated.

# Returns
- The result of the function evaluation.

# Example
```julia
pr = (objective="myFunctionName", ...)
x_values = [1.0, 2.0, 3.0]
t_value = 0.5
result = evaluateFunction(pr, x_values, t_value, arg1=value1, arg2=value2)
"""
function evaluateFunction(pr, x::Vector{Float64}, t::Float64; kwargs...)
    # Convert the string name to a function symbol
    func_sym = Symbol(pr.objective)
    
    # Dynamically call the function with the provided arguments and keyword arguments
    value = eval(:($func_sym($x, $t; $kwargs...)))
    
    return value
end