include("src/helperFunctions.jl")



"""
        create_algorithm_settings(; kwargs...)

Create a dictionary to store algorithm settings for optimization methods.

# Arguments
- `method::String="QuasiNewton"`: Specifies the optimization method.
- `maxiter::Int=Int(1e4)`: Maximum number of iterations.
- `gtol::Float64=1e-12`: Gradient tolerance for convergence.
- `dftol::Float64=1e-15`: Function value difference tolerance.
- `dxtol::Float64=1e-10`: Argument difference tolerance.
- `lambda::Int=1`: Lambda parameter.
- `lambdaMax::Int=100`: Maximum lambda value.
- `linesearch::String="StrongWolfe"`: Type of line search to use.
- `c1::Float64=1e-4`: Parameter for line search.
- `c2::Float64=0.9`: Another parameter for line search.
- `progress::Int=100`: Progress report frequency.

Additional keyword arguments can be passed to override the default settings.

# Returns
- `Dict{Symbol, Any}`: A dictionary with the algorithm settings.

# Examples
```julia
alg = create_algorithm_settings()
alg_custom = create_algorithm_settings(method="GradientDescent", maxiter=500)
```
"""
function create_algorithm_settings(;
        # method="QuasiNewton",
        # method = "ConjugateGradientDescent",
        # method = "GradientDescent",
        method = "NelderMead",
        # method = "GeneticAlgorithm",
        maxiter=Int(1e4),
        gtol=1e-12,
        dftol=1e-15,
        dxtol=1e-10,
        lambda=1,
        lambdaMax=100,
        linesearch="StrongWolfe",
        c1=1e-4,
        c2=0.9,
        progress=100,
        kwargs...
)

        # Create a dictionary with default values
        alg_settings = Dict(
                :method => method,
                :maxiter => maxiter,
                :gtol => gtol,
                :dftol => dftol,
                :dxtol => dxtol,
                :lambda => lambda,
                :lambdaMax => lambdaMax,
                :linesearch => linesearch,
                :c1 => c1,
                :c2 => c2,
                :progress => progress
        )

        # Update dictionary with any keyword arguments provided
        for (key, value) in kwargs
                alg_settings[key] = value
        end

        # Logic for method-specific modifications, similar to your original struct
        if alg_settings[:method] == "GradientDescent"
                # Specific adjustments for GradientDescent
                alg_settings[:progress] = 100
                
        elseif alg_settings[:method] == "QuasiNewton"
                # Adjustments for QuasiNewton
                alg_settings[:linesearch] = "StrongWolfe"
                alg_settings[:progress] = 1

        elseif alg_settings[:method] == "ConjugateGradientDescent"
                # Adjustments for ConjugateGradientDescent
                alg_settings[:linesearch] = "StrongWolfe"
                alg_settings[:c2] = 0.5
                alg_settings[:progress] = 1

        elseif alg_settings[:method] == "TrustRegion"
                # Adjustments for TrustRegion
                alg_settings[:linesearch] = "QuasiNewton-SR1"
                alg_settings[:c2] = -1.0
                alg_settings[:progress] = 1

        elseif alg_settings[:method] == "NelderMead"
                # Adjustments for NelderMead
                alg_settings[:linesearch] = "NA"
                alg_settings[:gtol] = "NA"
                alg_settings[:c1] = "NA"
                alg_settings[:c2] = "NA"
                alg_settings[:lambda] = "NA"
                alg_settings[:lambdaMax] = "NA"
                alg_settings[:DeltaTol] = 1e-12
                alg_settings[:Delta] = 100.0
                alg_settings[:progress] = 10
                alg_settings[:alpha] = 1.0
                alg_settings[:beta] = 0.5
                alg_settings[:gamma] = 2.0
                alg_settings[:delta] = 0.5

        elseif alg_settings[:method] == "GeneticAlgorithm"
                alg_settings[:maxiter] = 10000
                # alg_settings[:maxiter] = 10
                alg_settings[:progress] = alg_settings[:maxiter]/10
                # alg_settings[:progress] = 1
                alg_settings[:linesearch] = "NA"
                alg_settings[:gtol] = "NA"
                alg_settings[:c1] = "NA"
                alg_settings[:c2] = "NA"
                alg_settings[:lambda] = "NA"
                alg_settings[:lambdaMax] = "NA"
                alg_settings[:fvalRepeatTol] = alg_settings[:maxiter]/10
                alg_settings[:popSize] = 10
                alg_settings[:delta] = 0.3 # probability of mutation for each dimension of a point
                alg_settings[:Dist] = randn # probability distribution to choose value from for mutation
                # randn() has a mean of zero and a stddev of 1.
                alg_settings[:deviation] = 0.1 # magnitude of mutation allowed
                alg_settings[:parentsSurvive] = true
                # alg_settings[:parentsSurvive] = false
        else
                @warn "Bad condition."
        end

        return alg_settings
end

# Example usage
alg = create_algorithm_settings()



