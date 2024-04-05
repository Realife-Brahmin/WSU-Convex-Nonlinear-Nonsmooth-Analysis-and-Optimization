function solveForNextPGCGIterate(xk, gk, dk, rk, num, G, A, AAT;
        verbose=false)

        den = transpose(dk)*G*dk

        alphak = num/den
        xkp1 = xk + alphak*dk
        rkp1 = rk + alphak*G*dk
        vkp1 = AAT\(A*rkp1)
        gkp1 = rkp1 - transpose(A)*vkp1
        betak = transpose(rkp1)*gk/num

        dkp1 = -gkp1 + betak*dk

        # adding fevals and actions as oargs for consistency across my opt algos
        fevals_1PGCG = 0
        actions_1PGCG = Dict()

        return xkp1, gkp1, dkp1, rkp1, fevals_1PGCG, actions_1PGCG
end