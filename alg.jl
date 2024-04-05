include("src/helperFunctions.jl")


function create_algorithm_settings(;
        # problemType = "Unconstrained",
        problemType = "ECQP",
        # method = "ConjugateGradientDescent",
        # method = "GeneticAlgorithm",
        # method = "GradientDescent",
        # method = "NelderMead",
        method="QuasiNewton",
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
                :problemType => problemType,
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

        if alg_settings[:problemType] == "LP"
                @error "Haven't accounted what to do for LPs"
                
        elseif alg_settings[:problemType] == "ECQP"
                alg_settings[:method] = "ProjectedGradientCG"

        elseif alg_settings[:problemType] == "QP" 
                alg_settings[:method] = "ASQP"

        elseif alg_settings[:problemType] == "Unconstrained"
                # don't change method from whatever I specified

        else
                @error "Unknown method"
        end

        # Logic for method-specific modifications, similar to your original struct
        if alg_settings[:method] == "ConjugateGradientDescent"
                # Adjustments for ConjugateGradientDescent
                alg_settings[:linesearch] = "StrongWolfe"
                alg_settings[:c2] = 0.5
                alg_settings[:progress] = 1

        elseif alg_settings[:method] == "GradientDescent"
                # Specific adjustments for GradientDescent
                alg_settings[:progress] = 100

        elseif alg_settings[:method] == "GeneticAlgorithm"
                alg_settings[:maxiter] = 10000
                # alg_settings[:maxiter] = 10
                alg_settings[:progress] = alg_settings[:maxiter] / 10
                # alg_settings[:progress] = 1
                alg_settings[:linesearch] = "NA"
                alg_settings[:gtol] = "NA"
                alg_settings[:c1] = "NA"
                alg_settings[:c2] = "NA"
                alg_settings[:lambda] = "NA"
                alg_settings[:lambdaMax] = "NA"
                alg_settings[:fvalRepeatTol] = alg_settings[:maxiter] / 10
                alg_settings[:popSize] = 10
                alg_settings[:delta] = 0.3 # probability of mutation for each dimension of a point
                alg_settings[:Dist] = randn # probability distribution to choose value from for mutation
                # randn() has a mean of zero and a stddev of 1.
                alg_settings[:deviation] = 0.1 # magnitude of mutation allowed
                alg_settings[:parentsSurvive] = true
                # alg_settings[:parentsSurvive] = false

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

        elseif alg_settings[:method] == "ProjectedGradientCG"
                alg_settings[:linesearch] = "NA"
                alg_settings[:gtol] = "NA"
                alg_settings[:c1] = "NA"
                alg_settings[:c2] = "NA"
                alg_settings[:lambda] = "NA"
                alg_settings[:lambdaMax] = "NA"
                alg_settings[:tol] = 1e-6
                alg_settings[:progress] => progress
                
        elseif alg_settings[:method] == "QuasiNewton"
                # Adjustments for QuasiNewton
                alg_settings[:lambdaMax] = 1
                alg_settings[:linesearch] = "StrongWolfe"
                alg_settings[:progress] = 1


        elseif alg_settings[:method] == "TrustRegion"
                # Adjustments for TrustRegion
                alg_settings[:linesearch] = "QuasiNewton-SR1"
                alg_settings[:c2] = -1.0
                alg_settings[:progress] = 1

        else
                @warn "Bad condition."
        end

        return alg_settings
end

# Example usage
alg = create_algorithm_settings()



