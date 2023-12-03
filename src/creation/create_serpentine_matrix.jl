function create_serpentine_speed_matrix()
    very_low_speed = 0.001
    high_speed = 10
    mx, my = 256, 256
    matrix = very_low_speed * ones(mx, my)
    amplitude = mx / 20  # Amplitude of the sine wave
    frequency = 2  # How many "waves" you want in the matrix
    thickness = 10
    # Create a serpentine path by varying the x-coordinate with a sine function
    for y in 1:my
        # Sine function for smooth curves
        x = round(Int, amplitude * sin(frequency * 2 * pi * y / my) + mx / 2)
        # Set the width of the "snake"
        x_start = max(1, x - thickness)
        x_end = min(mx, x + thickness)
        matrix[x_start:x_end, y] .= high_speed
    end
    return matrix
end

speedData = create_serpentine_speed_matrix()
using Plots
Plots.heatmap(speedData, yflip = true)

# Now, to save the matrix to a CSV file
using CSV, DataFrames

func_dir = @__DIR__
root_dir = dirname(dirname(func_dir))
rawData_dir = joinpath(root_dir, "rawData")
speedMatrixID = "SpeedData_Serpentine"
ext = ".csv"
# Save the matrix to a CSV file
filename_csv = joinpath(rawData_dir, speedMatrixID*ext)
CSV.write(filename_csv, DataFrame(speedData, :auto), header=false)
