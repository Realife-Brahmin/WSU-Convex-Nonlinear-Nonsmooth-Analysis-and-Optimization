alg = (method = "GradientDescent",
        # maxiter = Int(1e5),
        maxiter = Int(1e4),
        # maxiter = Int(5),
        # maxiter = Int(50),
        ngtol = 1e-10,
        dftol = 1e-12,
        dxtol = 1e-10,
        lambda = 1,
        lambdaMax = 100,
        linesearch = "Armijo",
        # linesearch = "StrongWolfe",
        c1 = 1e-4, # Pg 33 (3.1 Step Length)
        c2 = 0.9,
        # progress = 50
        progress = 1
        );