include("types.jl")
include("nelderMead.jl")

function optimizeNM(pr; 
    verbose::Bool=false, 
    verbose_ls::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    log_txt = log_path*"log_"*string(pr.objective)*"_"*pr.alg.method*"_"*string(pr.alg.maxiter)*".txt"

    if isfile(log_txt)
        rm(log_txt)
    end # remove logfile if present for the run

    
    solverState = SolverStateType()

    fevals = 0
    alpha, beta, gamma, delta = pr.alg.alpha, pr.alg.beta, pr.alg.gamma, pr.alg.delta

    progress = pr.alg.progress
    maxiter = pr.alg.maxiter

    x0 = pr.x0
    xk = x0
    myprintln(verbose, "Starting with initial point x = $(xk).", log_path=log_txt)
    f = pr.objective
    p = pr.p
    fk = f(x0, p, getGradientToo=false)
    myprintln(verbose, "which has fval = $(fk)", log_path=log_txt)

    X00 = createInitialSimplexFromOnePoint(x0, deviationFactor=deviationFactor) # this simplex is currently unsorted

    X0 = sortSimplex(X00, f, p)
    solState = SolStateNMType(Xk=X0)

    # @pack! solState = fk

    keepIterationsGoing = true



    # function SolStateNMType(; k=0, Xkm1=zeros(0, 0), Xk=zeros(0, 0),
    #     Fkm1=zeros(0), Fk=zeros(0), Delta=0.0)
    #     return Dict(:k => k, :Xkm1 => Xkm1, :Xk => Xk, :Fkm1 => Fkm1, :Fk => Fk, :Delta => Delta)
    # end

    while keepIterationsGoing

        @unpack k = solverState
        
        printOrNot = verbose && ((k - 1) % progress == 0)
        printOrNot_ls = printOrNot & verbose_ls

        Xkp1, action = nelderMead(Xk, f, p)
        k += 1

        if k >= maxiter
            keepIterationsGoing = false
        elseif diameterSimplex < sizeTol
            keepIterationsGoing = false
        end

    end
end