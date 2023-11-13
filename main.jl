include("setup.jl")

println("You are currently using $(Threads.nthreads()) threads.")
println("Your machine has a total of $(Sys.CPU_THREADS) available threads.")

verbose = false
# verbose = true;
verbose_ls = false;
# verbose_ls = true;
# verbose_ls = verbose & verbose_ls
logging = true;
profiling = false;
# benchmarking = false;
benchmarking = true;

warmStart = true

# functionName = "dampedSHM";
# functionName = "drag"; functionName == "drag" ? verbose = false : verbose = verbose  
# functionName = "rosenbrock";
# functionName = "sphere";
# functionName = "TestFunction1";
# functionName = "TestFunction2";
# functionName = "TestFunction3";
# functionName = "signalDenoise";
functionName = "nnloss";

pr = include("src/objfuns/"*String(functionName)*".jl")

# Call the function with the initial problem setup
# res = @btime begin
@time begin
    res = warm_start_optimize(pr, verbose=verbose, verbose_ls=verbose_ls)
end

showresults(res)

if functionName == "nnloss"
    w_trained = res.xvals[:, end]
    classify = true
    p = pr.p
    p.params[:classify] = classify
    # pr = replace_field(pr, :x0, w_trained)
    # pr = replace_field(pr, :p, p)
    y_pred_tr = vec(nnloss(w_trained, p))
end
y_tr = pr.p.params[:classData]

# Binarize predictions based on the threshold
threshold = 0.5
binary_predictions = Int.(y_pred_tr .> threshold)
binary_actual = Int.(y_tr .> threshold)  # Convert to Int if not already

# Calculate the elements of the confusion matrix
true_positives = sum((binary_predictions .== 1) .& (binary_actual .== 1))
false_positives = sum((binary_predictions .== 1) .& (binary_actual .== 0))
true_negatives = sum((binary_predictions .== 0) .& (binary_actual .== 0))
false_negatives = sum((binary_predictions .== 0) .& (binary_actual .== 1))

# Create the confusion matrix
confusion_matrix_tr = [true_negatives false_negatives;
                    false_positives true_positives]

# Plot the confusion matrix
pTrain = Plots.heatmap(["Actual Negative", "Actual Positive"], ["Predicted Negative", "Predicted Positive"], confusion_matrix, color=:blues, aspect_ratio=1,
title = "Fit of FFNN on the Training Data")

# df_train, df_test = MLJ.partition(df, 0.7, rng=123);
# yTrain, XTrain = MLJ.unpack(df_train, ==(:Diagnosis));
yTest, XTest = MLJ.unpack(df_test, ==(:Diagnosis));

# trainData = Matrix(XTrain)
# classData = Vector(yTrain)

trainData = Matrix(XTest);
classData = Vector(yTest);
# trainData = Matrix(XTrain)[1, :]
# classData = Vector(yTrain)[1]

# params = Dict(:trainData=>trainData, :classData=>classData, :dims=>dims, :classify=>classify)

p.params[:trainData] = trainData
p.params[:classData] = classData

y_pred_test = vec(nnloss(w_trained, p))

y_test = pr.p.params[:classData]

binary_predictions = Int.(y_pred_test .> threshold)
binary_actual = Int.(y_test .> threshold)  # Convert to Int if not already

# Calculate the elements of the confusion matrix
true_positives = sum((binary_predictions .== 1) .& (binary_actual .== 1))
false_positives = sum((binary_predictions .== 1) .& (binary_actual .== 0))
true_negatives = sum((binary_predictions .== 0) .& (binary_actual .== 0))
false_negatives = sum((binary_predictions .== 0) .& (binary_actual .== 1))

# Create the confusion matrix
confusion_matrix_test = [true_negatives false_negatives;
                    false_positives true_positives]

# Plot the confusion matrix
pTest = Plots.heatmap(["Actual Negative", "Actual Positive"], ["Predicted Negative", "Predicted Positive"], confusion_matrix, color=:blues, aspect_ratio=1,
title = "Fit of FFNN on the Testing Data")
# objective = nnloss;

# pr = generate_pr(objective, w0, params=params)
# plotresults(res, savePlot=true)