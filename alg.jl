include("src/helperFunctions.jl")

"""
        AlgorithmSettings(; kwargs...)

A mutable struct for holding the algorithm settings for an optimization solver.

# Fields
- `method::String`: The optimization method to be used. Supported methods include "GradientDescent", "QuasiNewton", and "ConjugateGradientDescent".
- `maxiter::Int`: Maximum number of iterations for the solver.
- `gtol::Float64`: Gradient tolerance.
- `dftol::Float64`: Change in function value tolerance.
- `dxtol::Float64`: Change in solution vector tolerance.
- `lambda::Int`: Initial scaling parameter.
- `lambdaMax::Int`: Maximum scaling parameter.
- `linesearch::String`: Type of line search to be used, default is "StrongWolfe".
- `c1::Float64`: Armijo condition constant.
- `c2::Float64`: Wolfe condition constant.
- `progress::Int`: Interval for logging progress.

# Constructor

The constructor allows for customization of the algorithm settings through keyword arguments, with default values provided for each field.

# Notes
- Based on the chosen method, the constructor will automatically set some fields to suitable values, e.g., for "GradientDescent", it will set progress to 100.
- A warning will be raised if an unsupported method is provided.

# Example
```julia-repl
algSettings = AlgorithmSettings(method = "GradientDescent", maxiter = 10000)
"""
mutable struct AlgorithmSettings
        method::String
        maxiter::Int
        gtol::Float64
        dftol::Float64
        dxtol::Float64
        lambda::Int
        lambdaMax::Int
        linesearch::String
        c1::Float64
        c2::Float64
        progress::Int

        # Constructor with default values
        function AlgorithmSettings(; 
                # method="ConjugateGradientDescent",
                # method = "GradientDescent",
                method = "QuasiNewton",
                maxiter=Int(1e5),
                # maxiter=Int(1e4),
                # maxiter=Int(1e3),
                # maxiter = Int(5),
                # maxiter = Int(3),
                # maxiter = Int(1),
                gtol=1e-10,
                dftol=1e-12,
                dxtol=1e-10,
                lambda=1,
                lambdaMax=100,
                linesearch="StrongWolfe",
                c1=1e-4,
                c2=0.9,
                progress=100)

                alg = new(method, maxiter, gtol, dftol, dxtol, lambda, lambdaMax, linesearch, c1, c2, progress)

                if alg.method == "GradientDescent"
                        myprintln(true, "Method chosen: GradientDescent", log=false)
                        alg.progress = 100
                elseif alg.method == "QuasiNewton"
                        myprintln(true, "Method chosen: QuasiNewton", log=false)
                        alg.linesearch = "StrongWolfe"
                        # alg.progress = 100
                        alg.progress = 1
                elseif alg.method == "ConjugateGradientDescent"
                        myprintln(true, "Method chosen: ConjugateGradientDescent", log=false)
                        alg.linesearch = "StrongWolfe"
                        alg.c2 = 0.5
                        alg.progress = 1
                else
                        @warn "Bad condition."
                end

                return alg
        end
end

alg = AlgorithmSettings();
