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
                # maxiter=Int(1e5),
                # maxiter = Int(5),
                maxiter = Int(3),
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
                        myprintln(true, "Method chosen: GradientDescent")
                        alg.progress = 100
                elseif alg.method == "QuasiNewton"
                        myprintln(true, "Method chosen: QuasiNewton")
                        alg.linesearch = "StrongWolfe"
                        alg.progress = 1
                elseif alg.method == "ConjugateGradientDescent"
                        myprintln(true, "Method chosen: ConjugateGradientDescent")
                        alg.linesearch = "StrongWolfe"
                        alg.c2 = 0.5
                        alg.progress = 5
                else
                        @warn "Bad condition."
                end

                return alg
        end
end

alg = AlgorithmSettings();
#     # Example usage:
#     alg = AlgorithmSettings()
#     println(alg)

#     alg2 = AlgorithmSettings(method="GradientDescent", progress=10)
#     println(alg2)
