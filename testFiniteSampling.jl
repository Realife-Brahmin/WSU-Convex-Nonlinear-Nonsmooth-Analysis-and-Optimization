include("src/sampleSpace.jl")
n = 2
p = 62
# p = 5
# method = "Random"
method = "Halton"
# method = "Latin Hypercube"
discard = 20
seed = 1234
ss = sampleSpace(n, p, method=method, seed=seed, discard=discard)
plotSampledSpace(ss, method, view=true, savePlot=false)