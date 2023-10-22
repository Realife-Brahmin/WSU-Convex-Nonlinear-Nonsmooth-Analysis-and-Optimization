# using BenchmarkTools
# using DifferentialEquations
# using ForwardDiff: GradientConfig, Chunk, gradient!
using Plots
using Revise 
using Symbolics: gradient, jacobian, solve_for, simplify

# function rosenbrock(x, y)
#     a = one(eltype(x))
#     # b = 100 * a
#     result = zero(eltype(x))
#     # for i in 1:length(x)-1
#     #     result += (a - x[i])^2 + b*(x[i+1] - x[i]^2)^2
#     # end
#     fxy = (a+x.^2).^2 + 10(y-x.^2)^2
#     3 fxy
# end

xvals = -1:0.1:1;
yvals = -1:0.1:1;
# out = similar(x);
# cdfg1 = GradientConfig(rosenbrock, x, Chunk{1}());
# cdfg4 = GradientConfig(rosenbrock, x, Chunk{4}());
# cdfg10 = GradientConfig(rosenbrock, x, Chunk{10}());
# @time gradient!(out, rosenbrock, x, y, cdfg1)
@variables x y;
f(x, y) = (1+x^2)^2 + 10(y-x^2)^2;
display(f(x, y))
z = [f(i, j) for i in xvals, j in yvals];
surface!(xvals, yvals, z);

∇f(x, y) = gradient(f(x, y), [x, y]);
simplify(∇f(x, y))
# x★, y★ = Symbolics.solve_for(∇f(x, y), [x, y])
# solve_for(∇f ~ 0, [x, y])
∇²f(x, y) = simplify(jacobian(∇f(x, y), [x, y])); 

display(∇²f(x, y))

# I found these points manually
# yet to figure out how to make julia do that lol
x0, y0 = 0, 0;
x1, y1 = -sqrt(1/19), 1/19;
x2, y2 = sqrt(1/19), 1/19;
xlist = [x0, x1, x2]
ylist = [y0, y1, y2]
f.(xlist, ylist)
∇f.(xlist, ylist)
∇²f.(xlist, ylist)
