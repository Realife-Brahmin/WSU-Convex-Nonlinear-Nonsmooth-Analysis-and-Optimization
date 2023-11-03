function ArmijoBackracking(pr::NamedTuple, xk::Vector{Float64}, 
    pₖ::Vector{Float64};
    ρ = 0.5,
    itrMax::Int64=50,
    itrStart::Int64=1,
    verbose::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")
    
    fevals_ls = 0
    gevals_ls = 0
    obj = pr.objective
    p = pr.p
    c₁ = pr.alg.c1
    β = ρ^(itrStart-1)
    
    xnext = copy(xk)
    if pr.alg.method == "GradientDescent"
        fₖ = obj(xk, p, getGradientToo=false)
        ∇fₖ = -pₖ
    else
        fₖ, ∇fₖ = obj(xk, p)
        gevals_ls += 1
    end
    fevals_ls += 1
    
    fnext = fₖ
    log_txt = log_path*"log_"*string(pr.objective)*"_"*pr.alg.method*"_"*pr.alg.linesearch*"_"*string(pr.alg.maxiter)*".txt"
    itr_search_for_α = itrStart-1

    while itr_search_for_α ≤ itrMax
        xnext = xk + β * pₖ
        myprintln(verbose, "Let's try shifting x to $(xnext)", log_path=log_txt)
        
        comparison_val = fₖ + c₁ * β * pₖ'*∇fₖ
        fnext = obj(xnext, p, getGradientToo=false)
        fevals_ls += 1
        
        if fnext ≤ comparison_val
            myprintln(verbose, "Armijo condition satisfied for β = $(β)", log_path=log_txt)
            break;
        else
            myprintln(false, "Armijo condition NOT satisfied for β = $(β)", log=log, log_path=log_txt)
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