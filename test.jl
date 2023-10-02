using Optim
rosenbrock(x) =  (1.0 - x[1])^2 + 100.0 * (x[2] - x[1]^2)^2
result = Optim.optimize(rosenbrock, zeros(2), GradientDescent(), linesearch=LineSearches.Armijo)