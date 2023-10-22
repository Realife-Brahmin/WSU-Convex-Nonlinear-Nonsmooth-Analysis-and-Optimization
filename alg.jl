alg = (
        # method = "GradientDescent",
        method = "QuasiNewton",
        maxiter = Int(1e5),
        # maxiter = Int(1e4),
        # maxiter = Int(5),
        # maxiter = Int(50),
        gtol = 1e-10,
        dftol = 1e-12,
        # dftol = 1e-10,
        # dftol = 1e-8,
        dxtol = 1e-10,
        lambda = 1,
        lambdaMax = 100,
        # linesearch = "Armijo",
        linesearch = "StrongWolfe",
        c1 = 1e-4, # Pg 33 (3.1 Step Length)
        # c1 = 1e-2,
        c2 = 0.9,
        progress = 100
        );

if alg.method == "QasiNewton"
        alg.progress = 5
elseif alg.method == "GradientDescent"
        alg.progress = 100
end