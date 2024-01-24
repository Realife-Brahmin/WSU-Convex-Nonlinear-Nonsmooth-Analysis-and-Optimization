include("helperFunctions.jl")

function NelderMead(simplex, f::Function;
    alpha = 1.0,
    gamma = 2.0,
    beta = 0.5,
    delta = 0.5,
    verbose::Bool = false)
    n, p = size(simplex)
    action = "unselected"
    # it is assumed that the simplex is sorted, best (most optimal) point first
    xb, xw, xsw = simplex[:, 1], simplex[:, p], simplex[:, p-1] 
    xc = mean(simplex[:, 1:p-1])
    xr = reflect(xc, xw, alpha=alpha)
    if f(xr) < f(xb)
        # better point than the best
        xe = extend(xc, xw)
        
    return (simplex=simplex, action=action)
end

function reflect(xc, xw;
    alpha = 1.0)
    xr = xc + alpha*(xc-xw)
    return xr
end

function extend(xc, xw;
    gamma = 2.0)
    xe = xc + gamma*(xc-xw)
    return xe
end

xc = [1, 2, 3]
xw = [4, 5, 6]
# xr = reflect(xc, xw, alpha=0.5)