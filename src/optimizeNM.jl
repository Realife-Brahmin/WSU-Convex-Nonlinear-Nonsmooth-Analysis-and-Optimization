include("types.jl")

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
    solState = SolStateNMType(xk=pr.x0)

    fevals = 0
    alpha, beta, gamma, delta = pr.alg.alpha, pr.alg.beta, pr.alg.gamma, pr.alg.delta

    progress = pr.alg.progress
    maxiter = pr.alg.maxiter

    x0 = pr.x0
    xk = x0
    myprintln(verbose, "Starting with initial point x = $(xk).", log_path=log_txt)
    obj = pr.objective
    p = pr.p
    fk = obj(x0, p, getGradientToo=false)
    myprintln(verbose, "which has fval = $(fk)", log_path=log_txt)
    @pack! solState = fk

end