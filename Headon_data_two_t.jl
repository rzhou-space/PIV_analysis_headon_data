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


time0 = read_vti("E:\\PhD\\Headon_results\\headon_data\\tp0.vti")

time0[:, :, 1]

# Extract the single i-th layer from vti files and save the image data
# over all time points locally into a .h5 file.

using HDF5

function extract_layer(folder_path::String, i::Int64)
    
    # variable folder_path should be in type string.
    n_files = length(readdir(folder_path))
    
    # Storing data from a single layer over time in an Array. 
    layer_data = Matrix{Float64}[]
    
    #for file in readdir(folder_path)
    for t in 0:n_files-1
        # file is in type String. Open the file.
        cell_data = read_vti(folder_path*"\\tp"*string(t)*".vti")
        # Choose the i-th layer from the cell data and append it 
        # to layer_data.
        push!(layer_data, cell_data[:,:,i])
    end
    
    # Convert the layer_data to 3D array structure
    # (width, height, timepoint) dimensional.
    l_data = reduce((x,y) -> cat(x,y,dims=3), layer_data)

    # Save the data in .h5 file with package HDF5.jl
    h5open("headon_layer_"*string(i)*".h5", "w") do file
    write(file, "data", l_data)  
    end

end


string(1)

for i in 0:4-1
    print(i)
end

sort!(readdir("E:\\PhD\\Headon_results\\headon_data"))



using PyPlot

# Data from layer 2. 
extract_layer("E:\\PhD\\Headon_results\\headon_data", 2)

using HDF5

function read_h5(folder_path::String)
    h5open(folder_path, "r") do file
        read(file, "data")
    end
end

l_data = read_h5("headon_layer_5.h5")

size(l_data)[3]

imshow(l_data[:,:,10])

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

# Doing PIV for two time points. 

using Images

img1 = Float64.(Gray.(l_data[:, :, 5]))
img2 = Float64.(Gray.(l_data[:, :, 7]))

begin
subplot( 1, 2, 1 ); imshow( img1, cmap="gray" );
subplot( 1, 2, 2 ); imshow( img2, cmap="gray" );
end

using multi_quickPIV

VF, SN = multi_quickPIV.PIV( img1, img2 );

using PyPlot
U = VF[ 1, :, : ]; 
V = VF[ 2, :, : ]; 
# Relevant parameters that define the positions of the interrogation and search region
default_params = multi_quickPIV.setPIVParameters(); 
IA_size = multi_quickPIV._isize( default_params )[1:2]; # interrogation area size (in pixels)
IA_step = multi_quickPIV._step( default_params )[1:2]; # distance between consecutive interrogation areas (in pixels)
xgrid = [ (x-1)*IA_step[2] + div(IA_size[2],2) for y in 1:size(U,1), x in 1:size(U,2) ]; 
ygrid = [ (y-1)*IA_step[1] + div(IA_size[1],2) for y in 1:size(U,1), x in 1:size(U,2) ]; 

default_params 

begin 
imshow( img1 )
#PyPlot.quiver( xgrid, ygrid, V, -1 .* U, color="red" )
PyPlot.quiver( xgrid, ygrid, U, V, color="red" )
end

using Plots, PyPlot, multi_quickPIV

default_params = multi_quickPIV.setPIVParameters() 
IA_size = multi_quickPIV._isize( default_params )[1:2] # interrogation area size (in pixels)
IA_step = multi_quickPIV._step( default_params )[1:2] # distance between consecutive interrogation areas (in pixels)

for t in 1:2
    
    #IJulia.clear_output(true)
    PyPlot.clf()
    
    img1 = Float64.(Gray.(l_data[:, :, t]))
    img2 = Float64.(Gray.(l_data[:, :, t+1]))
    
    VF, SN = multi_quickPIV.PIV(img1, img2)
    
    U = VF[ 1, :, : ]
    V = VF[ 2, :, : ]
    # Relevant parameters that define the positions of the interrogation and search region
    xgrid = [ (x-1)*IA_step[2] + div(IA_size[2],2) for y in 1:size(U,1), x in 1:size(U,2) ]; 
    ygrid = [ (y-1)*IA_step[1] + div(IA_size[1],2) for y in 1:size(U,1), x in 1:size(U,2) ]; 
    
    #imshow(img1)
    
    #Plots.display(Plots.quiver(xgrid, ygrid, quiver=(V, -1 .* U) ))
    PyPlot.quiver( xgrid, ygrid, V, -1 .* U )
    PyPlot.pause(0.1)
    
end


using PyCall
@pyimport matplotlib.animation as anim
using PyPlot


A = randn(20,20,20)

fig = figure(figsize=(3,3))

function make_frame(i)
    #PyPlot.clf()
    #PyPlot.quiver(i, i, i, 1, color="red")
    if i==0
        print("Tutu :/")
    end
    imshow(Float64.(A[:, :, i+1]))
    imshow(A[:, :, i+2])
end
myanim = anim.FuncAnimation(fig, make_frame, frames=19,
                            repeat = false, interval=500)
myanim[:save]("test.gif", writer="pillow")


A = randn(20,20,20)

using PyCall
@pyimport matplotlib.animation as anim
using PyPlot
using Images
using multi_quickPIV

l_data = read_h5("headon_layer_5.h5")
#imshow(Gray.(l_data[:,:,1]))
default_params = multi_quickPIV.setPIVParameters() 
IA_size = multi_quickPIV._isize( default_params )[1:2]
IA_step = multi_quickPIV._step( default_params )[1:2]

fig = PyPlot.figure(figsize=(10, 10))

function make_frame(i) # Due to Python function i begins with 0!
    
    PyPlot.clf()
    
    img1 = Float64.(Gray.(l_data[:, :, i+1]))
    img2 = Float64.(Gray.(l_data[:, :, i+2]))
    
    VF, SN = multi_quickPIV.PIV(img1, img2)
    
    U = VF[ 1, :, : ]
    V = VF[ 2, :, : ]
    
    xgrid = [ (x-1)*IA_step[2] + div(IA_size[2],2) for y in 1:size(U,1), x in 1:size(U,2) ] 
    ygrid = [ (y-1)*IA_step[1] + div(IA_size[1],2) for y in 1:size(U,1), x in 1:size(U,2) ]
    
    PyPlot.quiver( xgrid, ygrid, V, -1 .* U )
end

myanim = anim.FuncAnimation(fig, make_frame, frames=size(l_data,3)-1, 
                            interval=400)
myanim[:save]("test3.gif", writer="pillow")
#myanim.save("test2.gif", writer="pillow", args=["-loop", '1']))

using PyCall
@pyimport matplotlib.animation as anim
using PyPlot
import IJulia

A = randn(20,20,20,2)

fig, axes = PyPlot.subplots(nrows=1, ncols=2, figsize=(7, 2.5))
ax1, ax2 = axes

function make_frame(i)
    ax1.clear()
    ax2.clear()
    ax1.imshow(A[:,:,i+1, 1])
    ax2.imshow(A[:,:,i+1, 2])
end

myanim = anim.FuncAnimation(fig, make_frame, frames=size(A,3), interval=20, blit=false)
myanim[:save]("test2.gif", writer="pillow")


# Could be useful to generate the dynamic animation after PIV maps.

using Plots
#@gif
for i in 1:20
    IJulia.clear_output(true)
    x = range(0,2*pi,1000); y = sin.(i*x)
    Plots.display(Plots.plot(x,y, color="red"))
    sleep(0.1)
end
println("Done!")

using GLMakie

using FileIO, GeometryTypes, Colors

xs = range(0, 5)

using GLMakie, multi_quickPIV

time = Observable(0)

img1 = @lift(Float64.(Gray.(l_data[:, :, $time+1])))
img2 = @lift(Float64.(Gray.(l_data[:, :, $time+2])))

VF, SN = multi_quickPIV.PIV(img1, img2)

U = VF[ 1, :, : ]
V = VF[ 2, :, : ]
# Relevant parameters that define the positions of the interrogation and search region
default_params = multi_quickPIV.setPIVParameters()
IA_size = multi_quickPIV._isize( default_params )[1:2] # interrogation area size (in pixels)
IA_step = multi_quickPIV._step( default_params )[1:2] # distance between consecutive interrogation areas (in pixels)
xgrid = [ (x-1)*IA_step[2] + div(IA_size[2],2) for y in 1:size(U,1), x in 1:size(U,2) ]
ygrid = [ (y-1)*IA_step[1] + div(IA_size[1],2) for y in 1:size(U,1), x in 1:size(U,2) ]

PyPlot.quiver( xgrid, ygrid, V, -1 .* U, color="red" )
#=
fig = lines(xs, ys_1, color = :blue, linewidth = 4,
    axis = (title = @lift("t = $(round($time, digits = 1))"),))
GLMakie.scatter!(xs, ys_2, color = :red, markersize = 15)
=#
framerate = 1
timestamps = range(0, 5, step=1/framerate)

record(fig, "animation.mp4", timestamps;
        framerate = framerate) do t
    time[] = t
end

using PyPlot 

l_data = read_h5("headon_layer_1.h5")

imshow(l_data[:, :, 15])

img1 = l_data[:, :, 15]

using Images

typeof(Float64.(Gray.(img1)))

typeof(Gray.(img1))

using ImageBinarization
using TestImages
using FileIO

#testimg = testimage("earth_apollo17.jpg")
testimg = load("F:\\Downloads\\cells.jpg")
imgb1 = binarize(Gray.(testimg), Intermodes())

imgb2 = binarize(img1, Polysegment())
imshow(imgb2)
typeof(imgb2)
#t = find_threshold(img1, Otsu(); nbins = 256)

using multi_quickPIV

pivimg1 = l_data[:, :, 1]
pivimg1_b = binarize(pivimg1, Polysegment())
pivimg2 = l_data[:, :, 2]
pivimg2_b = binarize(pivimg2, Polysegment())

# Doing PIV
params = multi_quickPIV.setPIVParameters()
IA_size = multi_quickPIV._isize(params)[1:2]
IA_step = multi_quickPIV._step(params)[1:2]

VF, SN = multi_quickPIV.PIV(pivimg1_b, pivimg2_b, params)

U = VF[1, :, :]
V = VF[2, :, :]
    
xgrid = [(x-1)*IA_step[2] + div(IA_size[2],2) for y in 1:size(U,1), x in 1:size(U,2)] 
ygrid = [(y-1)*IA_step[1] + div(IA_size[1],2) for y in 1:size(U,1), x in 1:size(U,2)]

PyPlot.imshow(imgb2)
#PyPlot.quiver(xgrid, ygrid, V, -1 .* U)
PyPlot.quiver( xgrid, ygrid, U, V, color="red" )


