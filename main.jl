# main.jl
include("setup.jl")

# global const JULIA_NUM_THREADS = 4;
println("You are currently using $(Threads.nthreads()) threads.")
println("Your machine has a total of $(Sys.CPU_THREADS) available threads.")

functionName = "dampedSHM";

# functionName = "dampedSHM_Parallel"
# functionName = "TestFunction1"

pr = generate_pr(functionName);


verbose = false
# verbose = true;
logging = true;

if logging
        if !isdir("./logging")
                println("Creating logging directory since it doesn't exist.")
                mkdir("./logging")
        else 
                println("No need to create logging directory, it already exists.")
                initialize_logging(overwrite=true)
        end
end

# @profile begin
# @btime begin
        @time begin
                res = optimize(pr, verbose=verbose, itrStart=7)
        end
        
showresults(res)
# plotresults(pr, res)



# ProfileView.view();

# Open a file in write mode
# f = open("./logging/profile_results.txt", "w")

# Redirect the output of Profile.print to the file
# Profile.print(f, mincount=5000)

# Close the file
# close(f)