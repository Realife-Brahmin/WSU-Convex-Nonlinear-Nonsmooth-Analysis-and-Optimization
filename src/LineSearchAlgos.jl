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

function StrongWolfeBisection(pr::NamedTuple, xk, pk; 
    α_min=0, α_max=1, 
    itrMax=50, tol=1e-8,
    itrStart::Int64=1,
    verbose::Bool=false,
    log::Bool=true,
    log_path::String="./logging/")

    log_txt = log_path*"log_"*string(pr.objective)*"_"*pr.alg.method*"_"*pr.alg.linesearch*"_"*string(pr.alg.maxiter)*".txt"
    c₁ = pr.alg.c1
    c₂ = pr.alg.c2
    obj = pr.objective
    p = pr.p
    ϕ0, ϕprime0 = ((f, g) -> (f, pk'*g))(obj(xk, p)...)
    fevals_ls = 1
    gevals_ls = 1
    iteration = 0
    xmid = xk
    ϕ_mid = ϕ0
    α_mid = 0
    strongWolfeSatisfied = false
    while iteration < itrMax && abs(α_max - α_min) > tol
        α_mid = 0.5 * (α_min + α_max)
        xmid = xk + α_mid * pk
        myprintln(verbose, "Let's try shifting x to $(xmid)", log_path=log_txt)
        # Compute function value and gradient at α_mid
        ϕ_mid = obj(xmid, p, getGradientToo=false)
        fevals_ls += 1
        # Check the Armijo condition
        if ϕ_mid > ϕ0 + c₁ * α_mid * ϕprime0
            α_max = α_mid
        else
            myprintln(verbose, "StrongWolfe1 condition satisfied for α = $(α_mid)", log_path=log_txt)

            ϕmid, ϕprime_mid = ((f, g) -> (f, pk'*g))(obj(xmid, p)...)
            fevals_ls += 1
            gevals_ls += 1
            # Check the curvature condition
            if abs(ϕprime_mid) <= c₂ * abs(ϕprime0)
                myprintln(verbose, "StrongWolfe2 condition satisfied for α = $(α_mid)", log_path=log_txt)
                strongWolfeSatisfied = true

                ans = (α=α_mid, x=xmid, f=ϕ_mid, gmag=ϕprime_mid, backtracks=iteration, fevals=fevals_ls, gevals=gevals_ls)

                return ans
            else

                myprintln(verbose, "StrongWolfe2 condition NOT satisfied for α = $(α_mid)", log_path=log_txt)
                α_min = α_mid
            end
        end
        
        iteration += 1
    end
    
    if abs(α_max - α_min) > tol
        @warn "Strong Wolfe line search did not converge in $itrMax iterations, Returning best α with tolerance gap of $(abs(α_max - α_min)) between αmax and αmin"
    end
    
    # ans = (α=α_mid, x=xmid, f=ϕ_mid, backtracks=iteration, fevals=fevals_ls, gevals=gevals_ls)
    ans = (α=α_mid, x=xmid, f=ϕ_mid, gmag=ϕprime_mid, backtracks=iteration, fevals=fevals_ls, gevals=gevals_ls)

    return ans
end


