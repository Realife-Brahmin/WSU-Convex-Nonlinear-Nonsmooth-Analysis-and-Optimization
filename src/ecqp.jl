using Parameters

function solveECQP(pr)
    
    if pr.problemType != "ECQP"
        @error "Cannot solve Problem of Type $(pr.problemType) as an Equality Constrained Quadratic Programming Problem"
    else
        @unpack G, A, c, b = pr
    end

    # q(x) = 1//2*xT*G*x + xT*c
    # s.t.
    # A*x = b
    n = size(G, 1)
    m = size(A, 1)
    KKT0 = [G -transpose(A); A zeros(m, m)] # the matrix containing relating x*, λ* with G, c, A, b

    RHS0 = [-c; b]

    wOpt = KKT0\RHS0

    xOpt, lambdaOpt = wOpt[1:n], wOpt[n+1:end]
    
    return xOpt, lambdaOpt
end