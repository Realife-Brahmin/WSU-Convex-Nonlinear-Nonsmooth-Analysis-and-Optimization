def rosenbrock(x,p,nargout):
    """
    Generalized n-dim rosenbrock function with steepness parameter p[0]
    """
    scale=p[0]
    n=len(x)
    f=0
    if nargout>1:
        g=np.zeros((n,1))
    for k in range(n-1):
        T=x[k+1,0]-x[k,0]**2
        S=1-x[k,0]
        f+=scale*T**2+S**2
        if nargout>1:
            g[k,0]=-4*scale*x[k,0]*T-2*S
            if k>0:
                g[k,0]+=2*scale*(x[k,0]-x[k-1,0]**2)
    if nargout>1:
        g[n-1,0]=2*scale*(x[n-1,0]-x[n-2,0]**2)
        return f,g
    else:
        return f
    
