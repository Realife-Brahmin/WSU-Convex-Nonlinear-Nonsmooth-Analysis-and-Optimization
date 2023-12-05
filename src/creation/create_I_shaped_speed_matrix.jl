function create_I_shaped_speed_matrix()
    very_low_speed = 0.001
    mx, my = 256, 256
    matrix = very_low_speed*ones(mx, my)
    xc, yc = mx/2, my/2
    width = mx/10
    length = 8*my/10

    matrix[round.(Int, xc-width/2:xc+width/2), round.(Int, yc-length/2:yc+length/2)] .= 10
    return matrix'
end


SpeedData_I = create_I_shaped_speed_matrix()
using Plots
Plots.heatmap(SpeedData_I, yflip = true)
speedMatrixID = "SpeedData_I"
ext = ".csv"
func_dir = @__DIR__
root_dir = dirname(dirname(func_dir))
rawData_dir = joinpath(root_dir, "rawData")
filename_csv = joinpath(rawData_dir, speedMatrixID*ext)
# Save the matrix to a CSV file
# CSV.write("SpeedData_Spiral.csv", DataFrame(spiral_speed_matrix, :auto))
CSV.write(filename_csv, DataFrame(SpeedData_I, :auto), header=false)