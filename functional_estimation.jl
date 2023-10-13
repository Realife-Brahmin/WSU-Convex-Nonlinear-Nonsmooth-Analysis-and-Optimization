include("setup.jl")

myprintln(true, "You are currently using $(Threads.nthreads()) threads.")
myprintln(true, "Your machine has a total of $(Sys.CPU_THREADS) available threads.")

# inClassExample()
