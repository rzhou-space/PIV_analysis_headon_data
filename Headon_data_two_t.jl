# Install ReadVTK with: import Pkg; Pkg.add("ReadVTK")

using ReadVTK

# Loading VTI files into Julia

function read_vti( filename, typ=Float64 )
    vtk_file       = VTKFile( filename )
    vtk_cell_data  = get_cell_data(vtk_file)
    vtk_data_array = vtk_cell_data["intensities"]
    data_size, _   = ReadVTK.get_wholeextent(vtk_data_array.vtk_file.xml_file, false)
    vtk_data       = get_data( vtk_data_array ); 
    data_reshaped  = reshape(vtk_data, (data_size .- 1 )...); 
    return typ.( data_reshaped ) 
end


# Extract the i-th layer from vti files.

function extract_layer(folder_path, i)
    # variable folder_path should be in type string.
    #n_files = length(readdir(folder_path))
    
    # Storing data from a single layer over time in an Array. 
    #layer_data = []
    layer_data = Matrix{Float64}[]
    
    for file in readdir(folder_path)
        # file is in type String. Open the file.
        cell_data = read_vti(folder_path * "\\" * file)
        # Choose the i-th layer from the cell data and append it 
        # to layer_data.
        push!(layer_data, cell_data[:,:,i])
        
    end
    return layer_data
end

    

cell_data = read_vti("F:\\Downloads\\headon_data\\tp0.vti")

using PyPlot

# Data from layer 2. 

l_data = extract_layer("F:\\Downloads\\headon_data", 2)

imshow(l_data[:,:,10])

# reduce((x,y) -> cat(x,y,dims=3), l_data)

size(l_data)

imshow(l_data[10])

# Completely destoryed the data structure. 

#=
using DelimitedFiles


writedlm("txt_file.txt", l_data)
=#

# JSON works BUT: the data type will NOT be saved and is therefore not
# efficient for further analysis!

#= 
using JSON

layer = 2
cell_images = l_data
data = Dict("layer" => layer, "cell_images" => cell_images)

open("foo.json", "w") do f
    JSON.print(f, data)
end

load_data = JSON.parsefile("C:\\Users\\Rutian Zhou\\PIV\\foo.json")
=#

load_data["cell_images"][10]

imshow(transpose(load_data["cell_images"][10]))  # The image is transformed! :/

using JLD

save("l_data.jdl", "data", l_data)

using HDF5

h5open("mydata_"*string(2)*".h5", "w") do file
    write(file, "A", l_data)  # alternatively, say "@write file A"
end

c = h5open("mydata.h5", "r") do file
    read(file, "A")
end

imshow(c[:,:,10])

# Doing PIV for two time points. 

# Could be useful to generate the dynamic animation after PIV maps.

using Plots
for i in 1:20
    IJulia.clear_output(true)
    x = range(0,2*pi,1000); y = sin.(i*x)
    Plots.display(Plots.plot(x,y, color="red"))
    sleep(0.1)
end
println("Done!")


