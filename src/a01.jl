using BenchmarkTools
# using DifferentialEquations
using ForwardDiff
using Revise
using Symbolics

@variables x y;
f(x, y) = (1+x^2)^2 + 10(y-x^2)^2;

∇f(x, y) = Symbolics.gradient(f(x, y), [x, y])
println(∇f(x, y))

# x★, y★ = Symbolics.solve_for(∇f(x, y), [x, y])
solve(∇f)
∇²f(x, y) = Symbolics.jacobian(∇f(x, y), [x, y]); 

∇²f(x, y)

∇f(0, 0)
∇²f(0, 0)