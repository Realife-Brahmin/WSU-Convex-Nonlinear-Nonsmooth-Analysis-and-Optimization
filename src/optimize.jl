# module optimize

# export optimize

function optimize(pr; 
    verbose::Bool=false, 
    log::Bool=true,
    log_path::String="./logging/",
    itrStart::Int64=1)

    log_txt = log_path*"log_"*string(pr.objective)*"_"*pr.alg.method*"_"*pr.alg.linesearch*"_"*string(pr.alg.maxiter)*".txt"
    if isfile(log_txt)
        rm(log_txt)
    end # remove logfile if present for the run
    # Initial settings
    fevals = 0
    gevals = 0
    dftol = pr.alg.dftol
    progress = pr.alg.progress
    maxiter = pr.alg.maxiter
    x0 = pr.x0
    x = x0

    myprintln(verbose, "Starting with initial point x = $(x).", log_path=log_txt)
    obj = pr.objective
    p = pr.p
    M = max(size(p.data, 1), 1)
    fnext = 1e10
    fₖ = obj(x0, p, getGradientToo=false)
    fevals += 1
    n = length(x)
    itr = 1
    fvals, αvals = [zeros(Float64, maxiter) for _ in 1:2]
    backtrackVals = zeros(Int64, maxiter, 1)
    xvals = zeros(Float64, n, maxiter)
    
    myprintln(true, "Begin with the solver:", log_path=log_txt)
    
    while abs(fnext - fₖ) ≥ dftol && itr ≤ maxiter
        printOrNot = verbose && (itr % progress == 0)
        # printOrNot = false
        myprintln(printOrNot, "Iteration $(itr):", log_path=log_txt)
        fₖ, ∇fₖ = obj(x, p)
        fevals += 1
        gevals += 1
        pₖ = findDirection(pr, ∇fₖ)
        α, x, fnext, backtrackNum, fevals_ls, gevals_ls = linesearch(pr, x, pₖ, itrStart=itrStart, verbose=printOrNot)
        
        myprintln(printOrNot, "Iteration $(itr): x = $(x) is a better point with new fval = $(fnext).", log_path=log_txt)

        fevals += fevals_ls
        gevals += gevals_ls
        fvals[itr] = fnext
        αvals[itr] = α
        backtrackVals[itr] = backtrackNum
        xvals[:, itr] = x
        itr += 1
    end
    
    if itr > maxiter
        converged = false
        statusMessage = "Failed to converge despite $(maxiter) iterations! 😢"
        myprintln(true, statusMessage, log_path=log_txt)
        @warn statusMessage
    else
        converged = true
        statusMessage = "Convergence achieved in $(itr) iterations 😄"
        myprintln(true, statusMessage, log_path=log_txt)
        # truncating arrays as they weren't filled to capacity
        fvals, αvals, backtrackVals = [arr[1:itr-1] for arr in (fvals, αvals, backtrackVals, xvals)]
        xvals = xvals[:, 1:itr-1]
    end
    
    res = (converged=converged, statusMessage=statusMessage, fvals=fvals, αvals=αvals, backtrackVals=backtrackVals, xvals=xvals, M=M, fevals=fevals, gevals=gevals)

    return res
end


function findDirection(pr::NamedTuple, ∇fnow::Vector{Float64};
    verbose::Bool=false)::Vector{Float64}
    method = pr.alg.method
    n = length(∇fnow)
    if method == "GradientDescent"
        # Bₖ = I(n)
        # pₖ = -Bₖ*∇fnow
        pₖ = -∇fnow
    elseif method == "ConjugateGradientDescent"
        @error "Currently not formulated for this method"
    elseif method == "QuasiNewton"
        @error "Currently not formulated for this method"
    else
        @error "Currently not formulated for this method"
    end

    return pₖ
end

function linesearch(pr::NamedTuple, xnow::Vector{Float64}, 
    pₖ::Vector{Float64};
    itrMax::Int64=50,
    itrStart::Int64=1,
    verbose::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")
    
    fevals_ls = 0
    gevals_ls = 0
    obj = pr.objective
    p = pr.p
    isStrongWolfe = (pr.alg.linesearch == "StrongWolfe")
    c₁ = pr.alg.c1
    c₂ = pr.alg.c2
    ρ = 0.5
    β = ρ^(itrStart-1)
    xnext = copy(xnow)
    if pr.alg.method == "GradientDescent"
        fₖ = obj(xnow, p, getGradientToo=false)
        ∇fₖ = -pₖ
    else
        fₖ, ∇fₖ = obj(xnow, p)
        gevals_ls += 1
    end
    fevals_ls += 1
    
    fnext = fₖ
    log_txt = log_path*"log_"*string(pr.objective)*"_"*pr.alg.method*"_"*pr.alg.linesearch*"_"*string(pr.alg.maxiter)*".txt"
    itr_search_for_α = itrStart-1

    while itr_search_for_α ≤ itrMax
        xnext .= xnow .+ β .* pₖ
        myprintln(verbose, "Let's try shifting x to $(xnext)", log_path=log_txt)
        
        comparison_val = fₖ + c₁ * β * dot(∇fₖ, pₖ)
        fnext = obj(xnext, p, getGradientToo=false)
        fevals_ls += 1
        
        if fnext ≤ comparison_val
            myprintln(verbose, "Armijo condition satisfied for β = $(β)", log_path=log_txt)
            if isStrongWolfe
                fnext, ∇fnext = obj(xnext, p)
                gevals_ls += 1
                fevals_ls += 1
                if abs(dot(∇fnext, pₖ)) > c₂*abs(dot(∇fₖ, pₖ))
                    myprintln(false, "Curvature condition NOT satisfied for β = $(β)", log_path=log_txt)
                    β *= ρ
                    itr_search_for_α += 1
                else
                    break
                end
            else
                break
            end
        else
            myprintln(false, "Armijo condition NOT satisfied for β = $(β)", log=log)
            β *= ρ
            itr_search_for_α += 1
        end
    end

    if itr_search_for_α > itrMax
        @error "Line Search failed at point x = $(xnext) despite $(itr_search_for_α) iterations."
    end

    α = β
    return (α=α, x=xnext, f=fnext, backtracks=itr_search_for_α, fevals=fevals_ls, gevals=gevals_ls) 
end

function linesearchSW(pr::NamedTuple, xnow::Vector{Float64}, 
    pₖ::Vector{Float64};
    itrMax::Int64=50,
    itrStart::Int64=1,
    verbose::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")
    
    fevals_ls = 0
    gevals_ls = 0
    obj = pr.objective
    p = pr.p
    isStrongWolfe = (pr.alg.linesearch == "StrongWolfe")
    itr_search_for_α = 1
    ϕ(α) = obj(xnow + α * pₖ, p, getGradientToo=false)
    dϕ(α) = dot(obj(xnow + α * pₖ, p)[2], pₖ)

    # Initial values
    α0 = 1.0
    ϕ0 = ϕ(0.0)
    dϕ0 = dϕ(0.0)

    # Perform the StrongWolfe line search
    α, ϕα = strongWolfe(ϕ, dϕ, α0, ϕ0, dϕ0)

    # Update x using the found α
    xnext = xnow + α * pₖ
    # fnext = obj(xnext, p, getGradientToo=false)
    # α = β
    return (α=α, x=xnext, f=ϕα, backtracks=itr_search_for_α, fevals=fevals_ls, gevals=gevals_ls) 
end

function strongWolfe(ϕ, dϕ, α0, ϕ0, dϕ0)
    @warn "Unwritten function, Returns non-useful values"
    α = α0
    fnext = ϕ(α)
    return (α=α, ϕα=fnext)
end
# end