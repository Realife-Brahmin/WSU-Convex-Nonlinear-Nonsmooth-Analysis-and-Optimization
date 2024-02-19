include("types.jl")
include("nelderMead.jl")

function optimizeNM(pr; 
    verbose::Bool=false, 
    verbose_ls::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    log_txt = log_path*"log_"*string(pr.objective)*"_"*pr.alg[:method]*"_"*string(pr.alg[:maxiter])*".txt"

    if isfile(log_txt)
        rm(log_txt)
    end # remove logfile if present for the run

    solverState = SolverStateNMType()

    @unpack alpha, beta, gamma, delta = pr.alg

    progress = pr.alg[:progress]
    maxiter = pr.alg[:maxiter]
    DeltaTol = pr.alg[:DeltaTol]

    x0 = pr.x0
    n = length(x0)
    xk = x0
    # doing this even though NM requires multiple f evals
    myprintln(verbose, "Starting with initial point x = $(xk).", log_path=log_txt)
    f = pr.objective
    pDict = pr.p
    fk = f(x0, pDict, getGradientToo=false)
    myprintln(verbose, "which has fval = $(fk)", log_path=log_txt)
    @unpack fevals = solverState
    fevals += 1
    @pack! solverState = fevals

    X00 = createInitialSimplexFromOnePoint(x0, deviationFactor=0.1) # this simplex is currently unsorted

    X0, F0 = sortSimplex(X00, f, pDict)
    Delta0 = simplexDiameter(X0)
    @unpack actions = solverState
    actions[:sort] += 1
    fevals += (n+1)
    Fk = F0
    @pack! solverState = fevals, actions

    solState = SolStateNMType(Xk=X0, Fk=Fk, Delta=Delta0)

    keepIterationsGoing = true

    while keepIterationsGoing

        @unpack k = solverState
        @unpack Xk, Fk = solState

        Xkm1, Fkm1 = Xk, Fk
        @pack! solState = Xkm1, Fkm1
        
        printOrNot = verbose && ((k - 1) % progress == 0)
        printOrNot_ls = printOrNot & verbose_ls

        Xkp1, action = nelderMead(Xk, f, pDict)
        k += 1

        Deltak = simplexDiameter(Xk)
        
        Xk = Xkp1

        @pack! solState = Deltak, Xk

        if k >= maxiter
            keepIterationsGoing = false
        elseif Deltak < DeltaTol
            keepIterationsGoing = false
        end

        @pack! solverState = k

    end
end