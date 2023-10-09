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