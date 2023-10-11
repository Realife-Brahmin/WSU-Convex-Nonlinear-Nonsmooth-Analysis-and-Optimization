function strongWolfe(f, xnow, p, pk, α0, ϕ0, ϕprime0; pr=pr)

    αk = α0
    ϕnext = ϕ0
    interpolationSteps = 0
    if !sufficientDecrease(f, xnow, p, pk, αk, ϕ0, ϕprime0, c₁=pr.alg.c1)
        ϕ_α0 = f(xnow+α0*pk, p)
        α1 = - (ϕprime0 * α0^2) / ( 2* (ϕ_α0 - ϕ0 - ϕprime0*α0) )
        αk = α1
        interpolationSteps += 1
    end

    ϕ_α1 = f(xnow+α1*pk, p)
    αnow = α1
    αprev = α0
    ϕnow = ϕ_α1
    ϕprev = ϕ_α0
    while !sufficientDecrease(f, xnow, p, pk, αnow, ϕ0, ϕprime0, c₁=pr.alg.c1)
        a, b = 1/( αprev^2 * αnow^2 * (αnow - αprev) ) * [αprev^2 - αnow^2; -αprev^3 αnow^3] * [ϕnow - ϕ0 - ϕprime0*αnow, ϕprev - ϕ0 - ϕprime0*αprev]

        αprev = αnow
        αnow = (-b + sqrt(b^2 - 3a*ϕprime0))/(3a)

        ϕprev = ϕnow
        ϕnow = f(xnow+αnow*pk, p)
        interpolationSteps += 1
    end

    soln = (αk=αk, ϕnext=ϕnext, interpolationSteps=interpolationSteps)
    return soln
end

function sufficientDecrease(f, xnow, p, pk, αk, ϕ0, ϕprime0; c₁=1e-4)
    ϕ_αk = f(xnow+αk*pk, p)
    if ϕ_αk ≤ ϕ0 + c₁*αk*ϕprime0
        return true
    else
        return false
    end
end

function strongWolfeBisection(f, g, xnow, pk, α_low, α_high, ϕ0, ϕprime0; c₁=1e-4, c₂=0.9, maxiter=100)
    iteration = 0
    
    while iteration < maxiter
        α_mid = 0.5 * (α_low + α_high)
        
        # Compute function value and gradient at α_mid
        ϕ_mid = f(xnow + α_mid * pk)
        ϕprime_mid = dot(g(xnow + α_mid * pk), pk)
        
        # Check the Armijo condition
        if ϕ_mid > ϕ0 + c₁ * α_mid * ϕprime0
            α_high = α_mid
        else
            # Check the curvature condition
            if abs(ϕprime_mid) <= c₂ * abs(ϕprime0)
                return α_mid
            else
                α_low = α_mid
            end
        end
        
        iteration += 1
    end
    
    throw(ConvergenceException("Strong Wolfe line search did not converge in $maxiter iterations"))
end

# Example use:
f(x) = sum((x - ones(Float64, length(x))).^2)
g(x) = 2 * (x - ones(Float64, length(x)))
xnow = [2.0, 2.0]
pk = [-1.0, -1.0]

α_low = 0.0
α_high = 1.0
ϕ0 = f(xnow)
ϕprime0 = dot(g(xnow), pk)

α = strongWolfeBisection(f, g, xnow, pk, α_low, α_high, ϕ0, ϕprime0)
