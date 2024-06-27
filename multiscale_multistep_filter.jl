begin # iterative approaches
    """
        Multiple iterative multiscale_multistep_filters can be very effective 
        to get segmentation on noisy data.

        IM1    = ( (20,20,5), )   
        IM2    = ( (10,10,2), )              
        rads   = ( ((0,0,0),(1,1,1)), )    
        steps  = ( ((1,1,1),(2,2,2)), )
        params = ( (IM1,rads,steps), (IM2,rads,steps), (IM1,rads,steps) ); 
        th     = 1.5

        segm = iterative_multiscalestep_segmentations( vol, params, th )
    """
    function iterative_multiscalestep_segmentation( vol, parameters, ths )
        
        segmentation = ones( Float32, size( vol ) ); 
        
        N_iterations = length( parameters ); 

        for i in 1:N_iterations
            params = parameters[i]
            th = ths[i]
            grid_size, scales, steps = params; 
            out_tmp = multiscale_multistep_filter( vol .* segmentation, grid_sizes=grid_size, scales=scales, steps=steps, typ=Float64 );
            segmentation .*= ( vol .^ th ) .< out_tmp
        end
        
        return segmentation
    end

    function iterative_multiscalestep_segmentation_rev( vol_og, parameters, ths )

        vol = copy( vol_og ); 
        vol_mean = sum( vol_og )/length( vol_og )
        
        segmentation = zeros( Float32, size( vol ) ); 
        
        N_iterations = length( parameters );

        for i in 1:N_iterations
            params = parameters[i]
            th = ths[i]
            grid_size, scales, steps = params; 
            out_tmp = multiscale_multistep_filter( vol, grid_sizes=grid_size, scales=scales, steps=steps, typ=Float64 );
            seg_tmp = ( vol .^ th ) .< out_tmp 
            seg_rev = ( seg_tmp .== 0 )
            vol[ seg_rev ] .^= (1/th)
            segmentation .= max.( segmentation, seg_tmp )
        end
        
        return segmentation, vol
    end

    function iterative_multiscalestep_filter( vol, parameters, th )
        
        segmentation = ones( Float32, size( vol ) ); 
        
        N_iterations = length( parameters ); 

        for i in 1:N_iterations-1
            params = parameters[i]
            grid_size, scales, steps = params; 
            out_tmp = multiscale_multistep_filter( vol .* segmentation, grid_sizes=grid_size, scales=scales, steps=steps, typ=Float64 );
            segmentation .*= ( vol .^ th ) .< out_tmp
        end

        params = parameters[end]
        grid_size, scales, steps = params; 
        out_tmp = multiscale_multistep_filter( vol .* segmentation, grid_sizes=grid_size, scales=scales, steps=steps, typ=Float64 );
        
        return out_tmp
    end
end


"""
    This filter highlights structures in the input data by comparing local patches around each pixel at different scales and steps. 

    For example, for any pixel at position "pos"...
    ... we will compare the intensity patterns within a 7x7 grid of pixels around "pos"... 
    ... to the intensity patterns in zoomed-out 7x7 grid around "pos"...
    ... the zoomed-out grid is sampled from a smoothed version of the input image with step size of 2 (for example). 

    Since the main idea is to compare patches at different scales/steps, this function accepts pairs of scales/steps as arguments.
    For instance, the code below results in the following 2 local patch comparisons around each pixel: 

    1-. ( 7x7 grid at scale 2 step 2 ) VS ( 7x7 grid at scale 4 step 4 )
    2-. ( 5x5 grid at scale 4 step 4 ) VS ( 5x5 grid at scale 8 step 6 )

    output = multiscale_multistep_filter( input... )
                                          grid_sizes = ( (7,7), (5,5) ),
                                          scales     = ( (2,4), (4,8) ),
                                          steps      = ( (2,4), (4,6) )
                                        )

    # NOTE: this filter is capable of dealing with anistropic resolutions by accepting tuples of values for each "grid_size", "scale"
    and "step"... The code below leads to the followign 2 comparisons:

    1-. ( 7x5 grid at scale (2,1) step (2,1) ) VS ( 7x5 grid at scale (4,2) step (4,2) )
    2-. ( 5x8 grid at scale (4,6) step (4,6) ) VS ( 5x8 grid at scale (8,8) step (6,6) )

    output = multiscale_multistep_filter( input... )
                                          grid_sizes = ( ((7,5),(7,5)), ((5,8),(5,8)) ),
                                          scales     = ( ((2,1),(4,2)), ((4,6),(8,8)) ),
                                          steps      = ( ((2,1),(4,2)), ((4,6),(6,6)) )
                                        )
"""
function multiscale_multistep_filter( inp::Array{<:Real,N}; 
                                      grid_sizes = ( (8, 8, 2), (8,8,2) ), # IM
                                      scales     = ( ((2,2,1),(3,3,1)), ((3,3,1),(4,4,2) ) ), # rads
                                      steps      = ( ((2,2,1),(3,3,1)), ((3,3,1),(4,4,2) ) ), # steps
                                      dim_downscale=Tuple(ones(Int,N)),
                                      compensate_borders = true,
                                      typ=Float32 
                                    ) where {N}

    println( "hi ")
    
    @assert length(scales) == length(steps) == length(grid_sizes)
    num_pairs = length(scales)

    # Filter output
    out = zeros( typ, size( inp ) );

    # Number of operations at each pixel, which is used for averaging the computed quantities
    Narr = zeros( typ, size( inp ) );
 
    # In case a scale is repeated... we only need to compute it once. Just in case., we will 
    # compute the scale space over the unique scales.This requires that we asign to each input
    # scale an index its corresponding scale in the scale-space.

    scales_vector = [ scales[j][i] for i in 1:2, j in 1:length(scales) ][:]
    unique_scales = unique( scales_vector ); 
    scale_indices = [ ( findfirst( x->x==scales[i][1], unique_scales), findfirst( x->x==scales[i][2], unique_scales) )  for i in 1:num_pairs ]

    # Creating the scale space
    scale_space = scale_space_( inp, unique_scales..., compensate_borders=compensate_borders, intA_typ=typ )

    # Squaring the scale-space
    scale_space_2 = scale_space .* scale_space

    # adding the results from the multiscale_step_filter at each pair of scale/steps
    for i in 1:num_pairs;

        scalestep_op!( out, Narr, 
                       scale_space, scale_space_2, 
                       grid_sizes[i], scale_indices[i], steps[i] ); 

    end

    return out ./ Narr
end

function multiscale_multistep_filter( inp::Array{<:Real,N}, 
                                      mask; 
                                      grid_sizes = ( (8, 8, 2), (8,8,2) ), # IM
                                      scales     = ( ((2,2,1),(3,3,1)), ((3,3,1),(4,4,2) ) ), # rads
                                      steps      = ( ((2,2,1),(3,3,1)), ((3,3,1),(4,4,2) ) ), # steps
                                      dim_downscale=Tuple(ones(Int,N)),
                                      compensate_borders = true,
                                      typ=Float32 
                                    ) where {N}
    
    @assert length(scales) == length(steps) == length(grid_sizes)
    num_pairs = length(scales)

    # Filter output
    out = zeros( typ, size( inp ) );

    # Number of operations at each pixel, which is used for averaging the computed quantities
    Narr = zeros( typ, size( inp ) );
 
    # In case a scale is repeated... we only need to compute it once. Just in case., we will 
    # compute the scale space over the unique scales.This requires that we asign to each input
    # scale an index its corresponding scale in the scale-space.

    scales_vector = [ scales[j][i] for i in 1:2, j in 1:length(scales) ][:]
    unique_scales = unique( scales_vector ); 
    scale_indices = [ ( findfirst( x->x==scales[i][1], unique_scales), findfirst( x->x==scales[i][2], unique_scales) )  for i in 1:num_pairs ]

    # Creating the scale space
    scale_space = scale_space_( inp, unique_scales..., compensate_borders=compensate_borders, intA_typ=typ )
    scale_space_mask = scale_space_( mask, unique_scales..., compensate_borders=compensate_borders, intA_typ=typ )

    # Squaring the scale-space
    scale_space_2 = scale_space .* scale_space

    # adding the results from the multiscale_step_filter at each pair of scale/steps
    for i in 1:num_pairs;

        scalestep_op!( out, Narr, 
                       scale_space, scale_space_2, scale_space_mask,
                       grid_sizes[i], scale_indices[i], steps[i] ); 

    end
    
    return out ./ Narr
end


##### UTILITY FUNCTION 

begin ##### MULTISCALE MULTISTEP PATCH FILTER

    # TODO: Replace imaginary by real ffts 
    """
        This function samples two patches at different scales/step around each pixel, and adds
        their "multiscale multistep score" to the output array. 

        See TODO for an illustrative explanation of what is actually computed with this filter.
    """
    function scalestep_op!( out::AbstractArray{T,N},             # output array 
                            Narr::AbstractArray{T,N},            # number of operations performed at each pixel, it will be used to compute an average
                            scale_space,                         # scale space of the input data
                            scale_space_2,                       # scale_space squared of the input data
                            grid_size::NTuple{N,Int},            # size of both patches 
                            sidx_pair::NTuple{2,Int},            # scale index for each one of two patches
                            step_pair::NTuple{2,NTuple{N,Int}};  # steps for each one of two patches
                            ) where {T,N}

        # Apart from being sampled at different scales, the two patches are sampled with different
        # "grid steps". The difference between the grid steps determines the kernel that we need
        # to convolve in order to compute the "multiscale multistep score" for each pixel. 

        step_dif    = max.( step_pair[2] .- step_pair[1], 1 );
        grid_kernel = zeros( T, 2 .* grid_size .* step_dif .+ 1 );
        grid_kernel[ Base.StepRange.( 1, step_dif, size(grid_kernel) )... ] .= 1
        
        scale_1_idx = sidx_pair[1]
        scale_2_idx = sidx_pair[2]

        # We want to multiply each pixel in "scale_2" by the values from all surrounding elements 
        # within the "grid_kernel" pattern. This is one half of the "multiscale multistep score" 

        sum_1   = ImageAnalysis.FFTConvolution_crop( scale_space[ axes( out )..., scale_1_idx ], grid_kernel )
        sum_1 .*= view( scale_space, ( axes( out )..., scale_2_idx )... )

        # The second part of the computation involves multiplying each pixel in "scale_2" .^ 2 by the
        # number of surrounding elements within the "grid_kernel" pattern. The pixels in "scale_2" 
        # near the borders require spatial attention... this is what "create_N_grid" is for.

        NN = create_N_grid( grid_size, step_dif, size(out) ) .- 1; 
        Narr .+= NN

        # All in all, the computation is: scale_2 .^ 2 .* N .- scale_2 .* conv( scale_1, grid_kernel )

        out .+= NN .* view( scale_space_2, ( axes( out )..., scale_2_idx )... ) .- sum_1; 

        # out .+= view( scale_space_2, ( axes( out )..., scale_2_idx )... ) .* sum_1 ./ NN

        return nothing
    end

    function scalestep_op_2!( out::AbstractArray{T,N},             # output array 
                            Narr::AbstractArray{T,N},            # number of operations performed at each pixel, it will be used to compute an average
                            scale_space,                         # scale space of the input data
                            scale_space_2,                       # scale_space squared of the input data
                            grid_size::NTuple{N,Int},            # size of both patches 
                            sidx_pair::NTuple{2,Int},            # scale index for each one of two patches
                            step_pair::NTuple{2,NTuple{N,Int}};  # steps for each one of two patches
                            ) where {T,N}

        # Apart from being sampled at different scales, the two patches are sampled with different
        # "grid steps". The difference between the grid steps determines the kernel that we need
        # to convolve in order to compute the "multiscale multistep score" for each pixel. 

        step_dif    = max.( step_pair[2] .- step_pair[1], 1 );
        grid_kernel = zeros( T, 2 .* grid_size .* step_dif .+ 1 );
        grid_kernel[ Base.StepRange.( 1, step_dif, size(grid_kernel) )... ] .= 1
        
        scale_1_idx = sidx_pair[1]
        scale_2_idx = sidx_pair[2]

        # We want to multiply each pixel in "scale_2" by the values from all surrounding elements 
        # within the "grid_kernel" pattern. This is one half of the "multiscale multistep score" 

        sum_1   = ImageAnalysis.FFTConvolution_crop( scale_space[ axes( out )..., scale_2_idx ], grid_kernel )

        # The second part of the computation involves multiplying each pixel in "scale_2" .^ 2 by the
        # number of surrounding elements within the "grid_kernel" pattern. The pixels in "scale_2" 
        # near the borders require spatial attention... this is what "create_N_grid" is for.

        NN = create_N_grid( grid_size, step_dif, size(out) ) .- 1; 
        Narr .+= NN

        # All in all, the computation is: scale_2 .^ 2 .* N .- scale_2 .* conv( scale_1, grid_kernel )

        out .+= view( scale_space, ( axes( out )..., scale_1_idx )... ) .- sum_1 ./ NN

        return nothing
    end

    function scalestep_op!( out::AbstractArray{T,N},             # output array 
                            Narr::AbstractArray{T,N},            # number of operations performed at each pixel, it will be used to compute an average
                            scale_space,                         # scale space of the input data
                            scale_space_2,                       # scale_space squared of the input data
                            scale_space_mask,                     
                            grid_size::NTuple{N,Int},            # size of both patches 
                            sidx_pair::NTuple{2,Int},            # scale index for each one of two patches
                            step_pair::NTuple{2,NTuple{N,Int}};  # steps for each one of two patches
                            ) where {T,N}

        # Apart from being sampled at different scales, the two patches are sampled with different
        # "grid steps". The difference between the grid steps determines the kernel that we need
        # to convolve in order to compute the "multiscale multistep score" for each pixel. 

        step_dif    = max.( step_pair[2] .- step_pair[1], 1 );
        grid_kernel = zeros( T, 2 .* grid_size .* step_dif .+ 1 );
        grid_kernel[ Base.StepRange.( 1, step_dif, size(grid_kernel) )... ] .= 1
        
        scale_1_idx = sidx_pair[1]
        scale_2_idx = sidx_pair[2]

        # We want to multiply each pixel in "scale_2" by the values from all surrounding elements 
        # within the "grid_kernel" pattern. This is one half of the "multiscale multistep score" 

        sum_1   = ImageAnalysis.FFTConvolution_crop( scale_space[ axes( out )..., scale_1_idx ], grid_kernel )
        sum_1 .*= view( scale_space, ( axes( out )..., scale_2_idx )... )

        # The second part of the computation involves multiplying each pixel in "scale_2" .^ 2 by the
        # number of surrounding elements within the "grid_kernel" pattern. The pixels in "scale_2" 
        # near the borders require spatial attention... this is what "create_N_grid" is for.

        NN = ImageAnalysis.FFTConvolution_crop( T.( scale_space_mask[ axes( out )..., scale_1_idx ] .> 0 ), grid_kernel )
        Narr .+= NN

        # All in all, the computation is: scale_2 .^ 2 .* N .- scale_2 .* conv( scale_1, grid_kernel )

        out .+= NN .* view( scale_space_2, ( axes( out )..., scale_2_idx )... ) .- sum_1; 

        return nothing
    end


    begin #### UTILS

        """
            convolution simulates 0-padding for pixels that are close to the border. This leads to 
            the situation where the output of the filter for pixels around the border incorporates
            less data points than for the other pixels. 
            
            This can be compensated by finding the number of kernel elements for each pixel, which
            can be done with a another convolution... but it can be computed more efficiently from 
            the step and gridsize of the convolved grid.
        """
        function create_N_grid( IM::NTuple{2,Int}, step_dif, inp_size )

            output = zeros( Float32, inp_size )

            steps  = step_dif .* IM
            c_max  = steps[2]+step_dif[2]
            r_max  = steps[1]+step_dif[1]
            corner = prod( IM .+ 1 )

            for c in 1:div( inp_size[2]+1, 2 ),
                r in 1:div( inp_size[1]+1, 2 )
                
                c_i = min( 2*IM[2]+1, div( c-1, step_dif[2] ) + ( IM[2] + 1 ) )
                r_i = min( 2*IM[1]+1, div( r-1, step_dif[1] ) + ( IM[1] + 1 ) )
                    
                val = c_i * r_i
                
                output[    r   ,    c    ] = val
                output[    r   , end-c+1 ] = val
                output[ end-r+1,    c    ] = val
                output[ end-r+1 ,end-c+1 ] = val
            end

            return output
        end

        function create_N_grid( IM::NTuple{3,Int}, step_dif, inp_size )

            output = zeros( Float32, inp_size )

            steps  = step_dif .* IM
            z_max  = steps[3]+step_dif[3]
            c_max  = steps[2]+step_dif[2]
            r_max  = steps[1]+step_dif[1]
            corner = prod( IM .+ 1 )

            for z in 1:div( inp_size[3]+1, 2 ), 
                c in 1:div( inp_size[2]+1, 2 ), 
                r in 1:div( inp_size[1]+1, 2 )
                
                z_i = min( 2*IM[3]+1, div( z-1, step_dif[3] ) + ( IM[3] +  1 ) )
                c_i = min( 2*IM[2]+1, div( c-1, step_dif[2] ) + ( IM[2] +  1 ) )
                r_i = min( 2*IM[1]+1, div( r-1, step_dif[1] ) + ( IM[1] +  1 ) )
                
                val = r_i * c_i * z_i 
                
                output[    r   ,    c   ,    z    ] = val
                output[    r   , end-c+1,    z    ] = val
                output[ end-r+1,    c   ,    z    ] = val
                output[ end-r+1 ,end-c+1,    z    ] = val
                output[    r   ,    c   , end-z+1 ] = val
                output[    r   , end-c+1, end-z+1 ] = val
                output[ end-r+1,    c   , end-z+1 ] = val
                output[ end-r+1 ,end-c+1, end-z+1 ] = val
            end

            return output
        end
    end
end

begin ##### SCALE-SPACE RELATED FUNCTIONS

    """
        The code below creates a scale space at 3 different scales: (2,2,1), (3,3,1) & (4,4,4). 

        Each scale contains a smoothing radius for each dimension of the input. Each dimensions may
        contain different values, allowing to adapt the smoothing to the resolution of the data. This
        has been included because it is common to have lower resolutions in fluorescence microscopy 
        along the Z dimension, making it to apply less smoothing along the Z axis.s

        scsp = create_scalespace( volume, (2,2,1), (3,3,1), (4,4,4) )
    """
    function scale_space_( input::AbstractArray{T,N}, 
                           unique_scales::Vararg{NTuple{N,Int},ARGC};
                           intA_typ=Float32, 
                           compensate_borders=false
                         ) where {T,N,ARGC}

        # Integral array of the input data, allowing us to compute average intensities around any pixel at arbitrary scales very efficiently.

        intA = ImageAnalysis.integralArray( input, typ=intA_typ );

        # Computing scale space
        num_scales    = length( unique_scales )
        scale_space   = zeros( T, size( input )..., num_scales ); 

        for i in 1:num_scales
            # Extracting a view to the current scale. Views are references/pointers to a region of an array, and are supposed to avoid unecessary memory copies.
            scale_view = view( scale_space, ( axes( input )..., i )... ) 

            # Computing the smoothed input at the current scale and storing it into the scale_view
            integral_local_avg!( intA, unique_scales[i], scale_view, compensate_borders ); 
        end

        return scale_space
    end

    begin ### INTEGRAL ARRAY UTILITIES

        # with bound checks 

        function integralArea( int_arr::Array{<:AbstractFloat,N}, TL, BR ) where {N}
            @assert all( TL .>= 1 ) && all( BR .< size( int_arr ) ) && all( TL .<= BR )
            return integralArea_unsafe( int_arr, TL, BR )
        end

        # avoiding bound checks... but you will crash julia if you go out-of-bounds

        function integralArea_unsafe( int_arr::Array{<:AbstractFloat,2}, TL, BR )
            TL   = TL .+ 1 .- 1;
            BR   = BR .+ 1;
            @inbounds area = int_arr[BR[1],BR[2]] - int_arr[BR[1],TL[2]] - int_arr[TL[1],BR[2]] + int_arr[TL[1],TL[2]]
            return area
        end

        function integralArea_unsafe( int_arr::Array{<:AbstractFloat,3}, TLF, BRB )
            TLF = TLF .+ 1 .- 1; 
            BRB = BRB .+ 1; 
            @inbounds area  = int_arr[BRB[1],BRB[2],BRB[3]] - int_arr[TLF[1],TLF[2],TLF[3]]
            @inbounds area -= int_arr[TLF[1],BRB[2],BRB[3]] + int_arr[BRB[1],TLF[2],BRB[3]] + int_arr[BRB[1],BRB[2],TLF[3]]
            @inbounds area += int_arr[BRB[1],TLF[2],TLF[3]] + int_arr[TLF[1],BRB[2],TLF[3]] + int_arr[TLF[1],TLF[2],BRB[3]]
            return area
        end

        """
            Computes the local average around the pixel/voxel at position "coord" ...
            ... within the square region of size 2 .+ rads .+ 1 around "coord" ...
            ... by using the integral array "intA" to compute the average intensity ...
            ... and stores the result at position "dest" in the array "out". 
        """
        function integral_local_avg!( intA::Array{T,N},         # integral array of the input
                                      coord::NTuple{N,Int},     # coordinates around which we wish to compute the mean
                                      rads::NTuple{N,Int},      # radii in all dimensions to compute local mean
                                      dest::NTuple{N,Int},      # destination coordinates for each "coords"
                                      out::AbstractArray{Y,N},  # destination array or view to an array
                                      sq_size=prod(2 .*rads.+1) # size of the local squre region, we can pass it as an argument to save computations
                                    ) where {T,Y,N}
            
            # Making sure TLF and BRB are within bounds
            TLF = max.( coord .- rads, 1 );
            BRB = min.( coord .+ rads, size(intA) .- 1 );

            # 
            sq_size = ( sq_size == nothing ) ? prod( BRB .- TLF .+ 1 ) : sq_size;

            # Computing the local sum of intensities with an integral sum
            num = integralArea_unsafe( intA, TLF, BRB );

            # Storing the average intensity in the desired destination
            # TODO: check that "dest" is within bounds of "out"? Or rather leave it unsafe?
            out[ dest... ] = num / sq_size;
            
            return nothing
        end

        """
            Instead of a single "coord" and "dest", this function accepts a list or iterable
            of coordinates and destinations. It iterated over this list and calls the function
            above (integral_local_avg!) for each pair of elements in "coords" and "dests".
            
        """
        function integral_local_avg!( intA, coords, rads, dests, out, compensate_borders::Bool=false; sq_size=prod( 2 .* rads .+  1 ) )
            
            @assert length( coords ) == length( dests ) "coords and dests must contain the same number of elements"    
            
            # reseting the output array
            out .= 0.0

            # computing the local average around each "coord" - "dest" 
            sq_size = ( compensate_borders ) ? nothing : prod( 2 .* rads .+  1 ); 
            for i in 1:length( coords );
                coord = Tuple( coords[i] )
                dest  = Tuple(  dests[i] )
                integral_local_avg!( intA, coord, rads, dest, out, sq_size )
            end
            
            return nothing    
        end

        """
            Simplified version of the previous function for the case when we wish to compute
            the local average of each coordinate in the input array, and we wish to store the
            results in the same order within an array of the same size as the input array.

            In other words, this function assumes that "coord" == "dest".
        """
        function integral_local_avg!( intA, rads, out, compensate_borders::Bool=false; sq_size=prod( 2 .* rads .+  1 ) )
            
            @assert size(intA) == size(out) .+ 1 "dimensions arent compatible"

            # reseting the output array
            out .= 0.0

            # computing the local average around each "coord" - "dest" 
            sq_size = ( compensate_borders ) ? nothing : prod( 2 .* rads .+  1 ); 
            for cartesian_coord in CartesianIndices( axes(out) )
                coord = Tuple( cartesian_coord )
                integral_local_avg!( intA, coord, rads, coord, out, sq_size )
            end

            return nothing    
        end
    end
end 


begin ### previous versions 
    
    # replace convolution to compute the number of grid 
    # accepts pairs of (rad-step) comparisons 
    function multiscale_multistep_filter_( inp::Array{<:Real,N}; 
                                        IM=(8,8,2),
                                        rads=((2,3), (3,4)),
                                        steps=((2,3), (3,4)),
                                        step0=0,
                                        dim_downscale=Tuple(ones(Int,N)),
                                        ) where {N}
        
        @assert length(rads) == length(steps)
        # the main idea of this algorithm is to compare intesity patterns at two different scales/steps at a time.
        # The first comparison will be between scale/step rads[1][1]/steps[1][1] and rads[1][2]/steps[1][2]...
        # The second comparison will be between scale/step rads[2][1]/steps[2][1] and rads[2][2]/steps[2][2]... etc
        num_pairs = length( rads )

        # Filter output
        out = zeros( Float32, size( inp ) ); 

        # Integral array of the input data, allowing us to compute average intensities around any pixel at arbitrary scales very efficiently.
        intA = ImageAnalysis.integralArray( inp, typ=Float32 );

        # Creating a scale-space of the input image at each radius in "rads". A scale-space is an array with multiple blurred versions of
        # the input at different scales. 

        # The same radius can be involved in multiple comparisons... and we want to avoid computing the scale space at the same scale multiple 
        # times. Therefore, we need to get the list of unique radii, and assign an index to this list for each element in "rads".
        unique_radii = unique( [ rads[i][j] for i in 1:num_pairs, j in 1:2 ] )
        rads_index   = [ ( findfirst( x->x==rads[i][1], unique_radii), findfirst( x->x==rads[i][2], unique_radii) )  for i in 1:num_pairs ]
        num_scales   = length(unique_radii); 
        scale_space  = zeros( size(inp)..., num_scales ); 

        for i in 1:num_scales
            # Extracting a view to the current scale. Views are references/pointers to a region of an array... and are supposed to avoid 
            # unecessary memory copies.
            scale_view = view( scale_space, ( axes( inp )..., i )... ) 
            # In 3D light-sheet recordings, the z-resolution is usually lower than xy resolutions. The code above is meant to correct for that, 
            # by allowing the user to downscale the size of "rads" in the z-dimension. It actually allows to downscale each dimension by a 
            # different amount. By default, nothing is downscaled. 
            scale_rad = Tuple( div.( ones(Int,N) .* unique_radii[i], dim_downscale ) );
            # Computing the smoothed input at the current scale and storing it into the scale_view
            integral_avg!( intA, scale_rad, scale_view ); 
        end

        # Squaring the scale-space
        scale_space_2 = scale_space .* scale_space

        # Computing the multiscale_multistep filter by convolving grid patterns at each scale of the scale-space
        ones_ = ones(  Float32, size(inp) ); 
        Ns    = zeros( Float32, size(inp) );

        for i in 1:num_pairs; 

            # As noted above, "dim_downscale" is there to correct for spatial resolution anisotropies
            step_1   = max.( 1, div.( Tuple( ones( Int, N ) ) .* steps[i][1], dim_downscale ) )
            step_2   = max.( 1, div.( Tuple( ones( Int, N ) ) .* steps[i][2], dim_downscale ) )
            step_dif = max.( step_2 .- step_1, 1 )

            grid_kernel = zeros( Float32, 2 .* IM .* step_dif .+ 1 );
            grid_kernel[ Base.StepRange.( ones(Int,N), step_dif, size(grid_kernel) )... ] .= 1
            
            scale_idx_1 = rads_index[i][1]
            scale_idx_2 = rads_index[i][2]

            # scale_2 .* conv( scale_1 )
            sum_1   = ImageAnalysis.FFTConvolution_crop( scale_space[ axes( inp )..., scale_idx_1 ], grid_kernel )
            sum_1 .*= view( scale_space , ( axes( inp )..., scale_idx_2 )... )

            # scale_2 .^ 2 .* N 
            NN = create_N_grid( IM, step_dif, size(inp) ); 
            out .+= NN .* view( scale_space_2, ( axes( inp )..., scale_idx_2 )... ) .- sum_1; 
            Ns  .+= NN
        end
        
        return out ./ Ns
    end


    # still uses convolution to compute the number of grid element around each point
    # doesn't accept pairs of (rad-step) comparisons
    function multiscale_multistep_filter__( inp::Array{<:Real,N}; 
                                        IM=(8,8,2), 
                                        rads=(2,3,4), 
                                        steps=(2,3,4), 
                                        step0=0,
                                        dim_downscale=Tuple(ones(Int,N)),
                                        ) where {N}
        
        @assert length(rads) == length(steps)

        # Filter output
        out  = zeros( Float32, size( inp ) ); 

        # Integral array of the input data, allowing us to compute average intensities around any pixel at arbitrary scales very efficiently.
        intA = ImageAnalysis.integralArray( inp, typ=Float32 );

        # Creating a scale-space of the input image at each radius in "rads". A scale-space is an array with multiple blurred versions of
        # the input at different scales.

        num_scales = length(rads); 
        avg_inps   = zeros( size(inp)..., num_scales ); 

        for i in 1:num_scales

            scale_view = view( avg_inps, ( axes( inp )..., i )... ) 
            scale_rad  = Tuple( div.( ones(Int,N) .* rads[i], dim_downscale ) );
            # In 3D light-sheet recordings, the z-resolution is usually lower than xy resolutions. The code above is meant to correct for that, 
            # by allowing the user to downscale the size of "rads" in the z-dimension. It actually allows to downscale each dimension by a 
            # different amount. By default, nothing is downscaled. 
            integral_avg!( intA, scale_rad, scale_view ); 
        end

        # Squaring the scale-space
        avg_inps2 = avg_inps .* avg_inps

        # Computing the multiscale_multistep filter by convolving grid patterns at each scale of the scale-space

        prev = Float32.( deepcopy( inp ) ); 
        out0 = ImageAnalysis.FFTConvolution_crop( inp, ones( 2 .* IM .+ 1 ) )
        
        for i in 1:num_scales; 
            
            step_dif = steps[i] - step0
            kern_    = zeros( Float32, 2 .* IM .* step_dif .+ 1 );    
            ones_pos = Base.StepRange.( ones(Int,N), step_dif, size(kern_) ); 
            kern_[ ones_pos... ] .= 1
            
            sum_ = ImageAnalysis.FFTConvolution_crop( prev, kern_ )

            scale_view  = view( avg_inps , ( axes( inp )..., i )... ) 
            scale_view2 = view( avg_inps2, ( axes( inp )..., i )... ) 

            out .+= scale_view2 .* sum( kern_ ) .- sum_ .* scale_view; 
            prev .= scale_view; 

            #=
            if N == 3
                out .+= avg_inps2[:,:,:,i] .* sum(kern_) .- sum_ .* avg_inps[:,:,:,i]
            else
                out .+= avg_inps2[:,:,i] .* sum(kern_) .- sum_ .* avg_inps[:,:,i]
            end
            
            if N == 3
                prev .= avg_inps[ :, :, :, i ]
            else
                prev .= avg_inps[ :, :, i ]
            end
            =#
            
            step0 = steps[i]
        end
        
        return out
    end

    # I do not know why this code is here...

    function linear2cartesian( lidx, input_size::Dims{2} )
        h = input_size[1]
        x = ceil( Int, lidx/h )
        y = lidx - (x-1)*h;
        return (y,x)
    end
    
    function linear2cartesian( lidx, input_size::Dims{3} )
        h = input_size[1]; 
        w = input_size[2];
        z = ceil( Int, lidx/(h*w) )
        x = ceil( Int, (lidx - (z-1)*h*w)/h )
        y = lidx - (x-1)*h - (z-1)*h*w;
        return (y,x,z)
    end

    #=
        # TEST 2D

        using PyPlot

        # 2D

        inp      = zeros( 201, 201 ); 
        ones_    = ones(  Float32, size(inp) ); 
        step_dif = ( 3, 3 );
        IM       = ( 7, 7 );
        N        = 2;

        grid_kernel = zeros( Float32, 2 .* IM .* step_dif .+ 1 );
        grid_kernel[ Base.StepRange.( ones(Int,N), step_dif, size(grid_kernel) )... ] .= 1

        @time NN   = ImageAnalysis.FFTConvolution_crop( ones_, grid_kernel )  # N
        @time NN_2 = create_N_grid( IM, step_dif, size(inp) )

        PyPlot.figure( dpi=200 )
        subplot( 1, 3, 1 ); PyPlot.title( "conv_grid" ); PyPlot.imshow( NN );
        subplot( 1, 3, 2 ); PyPlot.title( "const_grid" ); PyPlot.imshow( NN_2 );
        subplot( 1, 3, 3 ); PyPlot.title( "error" ); PyPlot.imshow( ( NN .- NN_2 ) ); colorbar( fraction=0.04 )

        # TEST 3D

        using PyPlot

        # 3D

        inp      = zeros( 200, 200, 100 ); 
        ones_    = ones(  Float32, size(inp) ); 
        step_dif = ( 3, 3, 1 );
        IM       = ( 7, 7, 2 );
        N        = 3;

        grid_kernel = zeros( Float32, 2 .* IM .* step_dif .+ 1 );
        grid_kernel[ Base.StepRange.( ones(Int,N), step_dif, size(grid_kernel) )... ] .= 1

        @time NNN  = ImageAnalysis.FFTConvolution_crop( ones_, grid_kernel ); 
        @time NN_3 = create_N_grid( IM, step_dif, size(inp) )

        println( "max difference = ", maximum( abs.( NNN .- NN_3 ) ) )
    =#

end


begin ### gradient multistep multiscale segmentation

    function multiscale_multistep_filter_dot( inp::Array{<:Real,N}; 
                                        grid_sizes = ( (8, 8, 2), (8,8,2) ), # IM
                                        scales     = ( ((2,2,1),(3,3,1)), ((3,3,1),(4,4,2) ) ), # rads
                                        steps      = ( ((2,2,1),(3,3,1)), ((3,3,1),(4,4,2) ) ), # steps
                                        dim_downscale=Tuple(ones(Int,N)),
                                        compensate_borders = true,
                                        typ=Float32, axis=1
                                        ) where {N}
        
        @assert length(scales) == length(steps) == length(grid_sizes)
        num_pairs = length(scales)

        # Filter output
        out = zeros( typ, size( inp ) );

        # Number of operations at each pixel, which is used for averaging the computed quantities
        Narr = zeros( typ, size( inp ) );
    
        # In case a scale is repeated... we only need to compute it once. Just in case., we will 
        # compute the scale space over the unique scales.This requires that we asign to each input
        # scale an index its corresponding scale in the scale-space.

        scales_vector = [ scales[j][i] for i in 1:2, j in 1:length(scales) ][:]
        unique_scales = unique( scales_vector ); 
        scale_indices = [ ( findfirst( x->x==scales[i][1], unique_scales), findfirst( x->x==scales[i][2], unique_scales) )  for i in 1:num_pairs ]

        # Creating the scale space
        scale_space = scale_space_( inp, unique_scales..., compensate_borders=compensate_borders, intA_typ=typ )

        # adding the results from the multiscale_step_filter at each pair of scale/steps
        for i in 1:num_pairs;

            # scalestep_op!( out, Narr, 
            scalestep_op_dot!( out, Narr, 
                            scale_space,
                            grid_sizes[i], scale_indices[i], steps[i],
                            ax=axis ); 

        end

        # return out ./ Narr
        return out, Narr
    end

    function scalestep_op_dot!( out::AbstractArray{T,N},             # output array 
                                Narr::AbstractArray{T,N},            # number of operations performed at each pixel, it will be used to compute an average
                                scale_space,                         # scale space of the input data
                                grid_size::NTuple{N,Int},            # size of both patches 
                                sidx_pair::NTuple{2,Int},            # scale index for each one of two patches
                                step_pair::NTuple{2,NTuple{N,Int}};  # steps for each one of two patches
                                ax=1
                                ) where {T,N}

        # Apart from being sampled at different scales, the two patches are sampled with different
        # "grid steps". The difference between the grid steps determines the kernel that we need
        # to convolve in order to compute the "multiscale multistep score" for each pixel. 

        if     ax==1
            grid_kernel = grid_ydot( T, step_pair, grid_size )
        elseif ax==2
            grid_kernel = grid_xdot( T, step_pair, grid_size )
        elseif ax==3
            grid_kernel = grid_zdot( T, step_pair, grid_size )
        else
            grid_kernel =  grid_dot( T, step_pair, grid_size )
        end
        step_dif = max.( step_pair[2] .- step_pair[1], 1 );

        scale_1_idx = sidx_pair[1]
        scale_2_idx = sidx_pair[2]

        # We want to multiply each pixel in "scale_2" by the values from all surrounding elements 
        # within the "grid_kernel" pattern. This is one half of the "multiscale multistep score" 

        dot = scale_space[ axes( out )..., scale_2_idx ] .* ImageAnalysis.FFTConvolution_crop( scale_space[ axes( out )..., scale_1_idx ], grid_kernel )
        # dot = ImageAnalysis.FFTConvolution_crop( scale_space[ axes( out )..., scale_1_idx ], grid_kernel )

        # The second part of the computation involves multiplying each pixel in "scale_2" .^ 2 by the
        # number of surrounding elements within the "grid_kernel" pattern. The pixels in "scale_2" 
        # near the borders require spatial attention... this is what "create_N_grid" is for.

        NN = create_N_grid( grid_size, step_dif, size(out) ); 
        Narr .+= NN

        # All in all, the computation is: scale_2 .^ 2 .* N .- scale_2 .* conv( scale_1, grid_kernel )

        out .+= dot

        return nothing
    end

    function grid_dot( typ, step_pair, grid_size  )
        step_dif    = max.( step_pair[2] .- step_pair[1], 1 );
        grid_kernel = zeros( typ, 2 .* grid_size .* step_dif .+ 1 );
        grid_kernel[ Base.StepRange.( 1, step_dif, size(grid_kernel) )... ] .= 1
        return grid_kernel
    end

    function grid_ydot( typ, step_pair, grid_size  )
        step_dif    = max.( step_pair[2] .- step_pair[1], 1 );
        grid_kernel = zeros( typ, 2 .* grid_size .* step_dif .+ 1 );
        grid_kernel[ Base.StepRange.( 1, step_dif, size(grid_kernel) )... ] .= 1
        for y in 1:size(grid_kernel,1), x in 1:size(grid_kernel,2), z in 1:size(grid_kernel,3) 
            mag = sqrt( sum( ( (y,x,z) .- div.( size(grid_kernel), 2 ) ) .^ 2 ) )
            ynorm = ( mag == 0 ) ? 0 : ( y - div( size(grid_kernel,1), 2 ) ) ./ mag
            grid_kernel[y,x,z] *= ynorm
        end
        return grid_kernel
    end

    function grid_xdot( typ, step_pair, grid_size  )
        step_dif    = max.( step_pair[2] .- step_pair[1], 1 );
        grid_kernel = zeros( typ, 2 .* grid_size .* step_dif .+ 1 );
        grid_kernel[ Base.StepRange.( 1, step_dif, size(grid_kernel) )... ] .= 1
        for y in 1:size(grid_kernel,1), x in 1:size(grid_kernel,2), z in 1:size(grid_kernel,3) 
            mag = sqrt( sum( ( (y,x,z) .- div.( size(grid_kernel), 2 ) ) .^ 2 ) )
            xnorm = ( mag == 0 ) ? 0 : ( x - div( size(grid_kernel,2), 2 ) ) ./ mag
            grid_kernel[y,x,z] *= xnorm
        end
        return grid_kernel
    end

    function grid_zdot( typ, step_pair, grid_size  )
        step_dif    = max.( step_pair[2] .- step_pair[1], 1 );
        grid_kernel = zeros( typ, 2 .* grid_size .* step_dif .+ 1 );
        grid_kernel[ Base.StepRange.( 1, step_dif, size(grid_kernel) )... ] .= 1
        for  y in 1:size(grid_kernel,1), x in 1:size(grid_kernel,2), z in 1:size(grid_kernel,3) 
            mag = sqrt( sum( ( (y,x,z) .- div.( size(grid_kernel), 2 ) ) .^ 2 ) )
            znorm = ( mag == 0 ) ? 0 : ( z - div( size(grid_kernel,3), 2 ) ) ./ mag
            grid_kernel[y,x,z] *= znorm
        end
        return grid_kernel
    end

    function high_derivative_y( scale, step, vol, v1=true; intA_typ=Float32 )
        if v1
            h = size( vol, 1 ); 
            der = zeros( Float32, size( vol ) )
            der1 = vol[ 1:h-2*step   , :, : ] .- vol[ 1+step:h-step, :, : ] 
            der2 = vol[ 1+step:h-step, :, : ] .- vol[ 1+2*step:h, :, : ] 
            der[ 1+step:h-step, :, :, ] .= sign.( der2 .- der1 ) .* max.( abs.( der1 ), abs.( der2 ) )
            scsp = scale_space_( der, scale  );
            return scsp[ :, :, :, 1 ] 
        else
            h = size( vol, 1 ); 
            scsp = scale_space_( vol, scale, intA_typ=intA_typ  )[:,:,:,1];
            der  = zeros( intA_typ, size( vol ) )
            der1 = scsp[ 1:h-2*step   , :, : ] .- scsp[ 1+step:h-step, :, : ] 
            der2 = scsp[ 1+step:h-step, :, : ] .- scsp[ 1+2*step:h, :, : ] 
            der[ 1+step:h-step, :, :, ] .= sign.( der2 .- der1 ) .* max.( abs.( der1 ), abs.( der2 ) )
            return der
        end
    end
    
    function high_derivative_x( scale, step, vol, v1=true; intA_typ=Float32 )
        if v1
            w = size( vol, 2 ); 
            der = zeros( Float32, size( vol ) )
            der1 = vol[ :, 1:w-2*step   , : ] .- vol[ :, 1+step:w-step, : ] 
            der2 = vol[ :, 1+step:w-step, : ] .- vol[ :, 1+2*step:w, : ] 
            der[ :, 1+step:w-step, : ] .= sign.( der2 .- der1 ) .* max.( abs.( der1 ), abs.( der2 ) )
            scsp = scale_space_( der, scale  );
            return scsp[ :, :, :, 1 ]     
        else
            w = size( vol, 2 ); 
            scsp = scale_space_( vol, scale, intA_typ=intA_typ  )[:,:,:,1];
            der  = zeros( intA_typ, size( vol ) )
            der1 = scsp[ :, 1:w-2*step   , : ] .- scsp[ :, 1+step:w-step, : ] 
            der2 = scsp[ :, 1+step:w-step, : ] .- scsp[ :, 1+2*step:w, : ] 
            der[ :, 1+step:w-step, : ] .= sign.( der2 .- der1 ) .* max.( abs.( der1 ), abs.( der2 ) )
            return der 
        end
    end
    
    function high_derivative_z( scale, step, vol, v1=true; intA_typ=Float32 )
        if v1
            z = size( vol, 3 ); 
            der = zeros( Float32, size( vol ) )
            der1 = vol[ :, :, 1:z-2*step    ] .- vol[ :, :, 1+step:z-step ] 
            der2 = vol[ :, :, 1+step:z-step ] .- vol[ :, :, 1+2*step:z ] 
            der[ :, :, 1+step:z-step ] .= sign.( der2 .- der1 ) .* max.( abs.( der1 ), abs.( der2 ) )
            scsp = scale_space_( der, scale  );
            return scsp[ :, :, :, 1 ]   
        else
            z = size( vol, 3 ); 
            scsp = scale_space_( vol, scale, intA_typ=intA_typ  )[:,:,:,1];
            der  = zeros( intA_typ, size( vol ) )
            der1 = scsp[ :, :, 1:z-2*step    ] .- scsp[ :, :, 1+step:z-step ] 
            der2 = scsp[ :, :, 1+step:z-step ] .- scsp[ :, :, 1+2*step:z ] 
            der[ :, :, 1+step:z-step ] .= sign.( der2 .- der1 ) .* max.( abs.( der1 ), abs.( der2 ) )
            return der 
        end  
    end
    
    function multiscale_multigrid_gradient_dots( vol_i; 
                                                 grad_smoothing = (4,4,2), 
                                                 grad_distance = (6,6,3),
                                                 grid_size = ( ( 2, 2, 1 ), ), 
                                                 scale = ( ((0,0,0), (1,1,1)), ), 
                                                 steps = ( ((1,1,1), (2,2,2)), ), 
                                                 v1=true, 
                                                 norm_ders=false, typ=Float32
                                               )
        
        der_1 = high_derivative_y( grad_smoothing, grad_distance[1], vol_i, v1, intA_typ=typ )
        der_2 = high_derivative_x( grad_smoothing, grad_distance[2], vol_i, v1, intA_typ=typ )
        der_3 = high_derivative_z( grad_smoothing, grad_distance[3], vol_i, v1, intA_typ=typ )
        
        if norm_ders
            mags = sqrt.( der_1 .^ 2 .+ der_2 .^ 2 .+ der_3 .^ 2 )
            der_1[ mags .> 0 ] ./= mags[ mags .> 0 ]
            der_2[ mags .> 0 ] ./= mags[ mags .> 0 ]
            der_3[ mags .> 0 ] ./= mags[ mags .> 0 ]
        end
    
        sum_u, _ = multiscale_multistep_filter_dot( der_1, grid_sizes=grid_size, scales=scale, steps=steps, axis=1 ); 
        sum_v, _ = multiscale_multistep_filter_dot( der_2, grid_sizes=grid_size, scales=scale, steps=steps, axis=2 ); 
        sum_w, _ = multiscale_multistep_filter_dot( der_3, grid_sizes=grid_size, scales=scale, steps=steps, axis=3 ); 
    
        dots = ( sum_u .+ sum_v .+ sum_w )
        
        return dots, sum_u, sum_v, sum_w, der_1, der_2, der_3
    end
    
    function multiscale_multigrid_gradient_segment( vol_i; 
                                                    grad_smoothing = (4,4,2), 
                                                    grad_distance = (6,6,3),
                                                    grad_grid_size = ( (2,2,1), ), 
                                                    grad_scale = ( ((0,0,0), (1,1,1)), ), 
                                                    grad_steps = ( ((1,1,1), (2,2,2)), ),
                                                    seg_grid_size = ( (3,3,1), ),
                                                    seg_scale = ( ((0,0,0), (1,1,0)), ),
                                                    seg_steps = ( ((1,1,1), (2,2,1)), ), 
                                                    th = 1.1, v1=true, typ=Float64, norm_ders=false
                                                  )
        
        dots, sum_u, sum_v, sum_w, der_1, der_2, der_3 = multiscale_multigrid_gradient_dots( vol_i,                                              
                                                   grad_smoothing = grad_smoothing, 
                                                   grad_distance  = grad_distance,
                                                   grid_size = grad_grid_size, 
                                                   scale = grad_scale, 
                                                   steps = grad_steps, v1=v1, typ=typ, norm_ders=norm_ders ); 
    
        ii  = max.( 0, dots ) .+ 1 
        out = multiscale_multistep_filter( ii, grid_sizes=seg_grid_size, scales=seg_scale, steps=seg_steps, typ=Float64 );
        seg = out .> ii .^ th

        return seg, dots, out, sum_u, sum_v, sum_w, der_1, der_2, der_3
    end


end