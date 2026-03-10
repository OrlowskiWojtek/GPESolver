using GLMakie
using LaTeXStrings
# This files includes sets of plots and animations that require time evolution of BEC
# TODO refactor from common code

include("context.jl")
include("files.jl")
include("units.jl")

set_theme!(theme_latexfonts(), rowgap = 0, colgap = 0)

function animate_iso_bce(data_dir::String, output_file::String)
    files = filter(f -> occursin("wavefunction_", f) && endswith(f, ".gpe.dat"), readdir(data_dir))
    files = sort(files, by = f -> parse(Int, split(split(f, "_")[2], ".")[1]))
    n_frames = length(files)

    fig = Figure(size = (1024*2, 768), backgroundcolor = :gray97)
    ax = Axis3(fig[1, 1],
               xlabel="X [nm]",
               ylabel="Y [nm]",
               zlabel="Z [nm]",
               aspect=:data,
            )

    BCEContext = load_from_text(joinpath(data_dir, files[1]))
    rho = Observable(Array{Float64,3}(abs.(BCEContext.psi)))
    x, y, z = BCEContext.x, BCEContext.y, BCEContext.z

    max_val = maximum(rho[])
    isovals = [0.2 * max_val, 0.5 * max_val, 0.8 * max_val]
    alphas = [0.4, 0.7, 0.9]

    for (i, isovalue) in enumerate(isovals)
        volume!(ax,
                (x[begin], x[end]),
                (y[begin], y[end]),
                (z[begin], z[end]),
                rho,
                algorithm = :iso,
                isovalue = isovalue,
                alpha = alphas[i],
                colormap = :YlGn,
                transparency = true,
                isorange = 0.1 * max_val)
    end

    record(fig, output_file, 1:n_frames; framerate=10) do frame
        println("loading frame", joinpath(data_dir, files[frame]))
        rho[] = abs.(load_from_text(joinpath(data_dir, files[frame])).psi)
    end
end

function plot_local_maxima_evolution(psi_vec::Vector{IsoBECContext})
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "x [nm]", ylabel = "y [nm]")

    n_frames = length(psi_vec)
    slice = interpolate_slice(get_BEC_slice(psi_vec[begin]))
    maxima = find_local_maxima(slice)
    coords = get_coordinates(slice, maxima)

    n_becs = length(coords)
    all_maxima = [LocalMaximaPhysical[] for n in 1:n_becs]

    for (idx, context) in enumerate(psi_vec)
        slice = interpolate_slice(get_BEC_slice(context))
        maxima = find_local_maxima(slice)
        coords = get_coordinates(slice, maxima)
        println("Processing context nr: ", idx)

        for (bec_idx, coord) in enumerate(coords)
            push!(all_maxima[bec_idx], coord)
        end 
    end

    frames = collect(1:length(psi_vec)) .* STEPS_PER_SAVE * time_au_to_ms(TIME_STEP_AU)
    for lmax in all_maxima
        xs = [m.x for m in lmax]
        ys = [m.y for m in lmax]

        scatter!(ax, xs, ys, color = frames, markersize = 5)
    end

    Colorbar(fig[1, 2], label="Time [ms]", limits = (0, frames[end]))

    #xlims!(ax, (psi_vec[begin].x[begin], psi_vec[begin].x[end]))
    #ylims!(ax, (psi_vec[begin].y[begin], psi_vec[begin].y[end]))

    return fig
end

function plot_local_maxima_evolution(slice_vec::Vector{IsoBECSlice})
    fig = Figure()
    ax = Axis(fig[1,1], xlabel = "x [nm]", ylabel = "y [nm]")

    n_frames = length(slice_vec)
    slice = interpolate_slice(slice_vec[begin])
    maxima = find_local_maxima(slice)
    coords = get_coordinates(slice, maxima)

    n_becs = length(coords)
    all_maxima = [LocalMaximaPhysical[] for n in 1:n_becs]

    for (idx, _slice) in enumerate(slice_vec)
        slice = interpolate_slice(_slice)
        maxima = find_local_maxima(slice)
        coords = get_coordinates(slice, maxima)
        println("Processing slice nr: ", idx)

        for (bec_idx, coord) in enumerate(coords)
            push!(all_maxima[bec_idx], coord)
        end 
    end

    frames = collect(1:length(slice_vec)) .* STEPS_PER_SAVE * time_au_to_ms(TIME_STEP_AU)
    for lmax in all_maxima
        xs = [m.x for m in lmax]
        ys = [m.y for m in lmax]

        scatter!(ax, xs, ys, color = frames, markersize = 5)
    end

    Colorbar(fig[1, 2], label="Time [ms]", limits = (0, frames[end]))

    return fig
end

function plot_local_maxima_coordinates(psi_vec::Vector{IsoBECContext})
    fig = Figure()
    ax_x = Axis(fig[1,1], ylabel = "x [nm]", xlabel = "t [ms]")
    ax_y = Axis(fig[2,1], ylabel = "y [nm]", xlabel = "t [ms]")

    n_frames = length(psi_vec)
    slice = interpolate_slice(get_BEC_slice(psi_vec[begin]))
    maxima = find_local_maxima(slice)
    coords = get_coordinates(slice, maxima)

    n_becs = length(coords)
    all_maxima = [LocalMaximaPhysical[] for n in 1:n_becs]

    for (idx, context) in enumerate(psi_vec)
        println("Processing context: ", idx)
        slice = interpolate_slice(get_BEC_slice(context); size = (500, 500))
        maxima = find_local_maxima(slice)
        coords = get_coordinates(slice, maxima)

        sorted_cords = sort(coords, by = c -> (c.x, c.y))

        for (bec_idx, coord) in enumerate(sorted_cords)
            push!(all_maxima[bec_idx], coord)
        end 
    end

    times = collect(1:length(psi_vec)) .* STEPS_PER_SAVE * time_au_to_ms(TIME_STEP_AU)
    for (idx, lmax) in enumerate(all_maxima)
        xs = [m.x for m in lmax]
        ys = [m.y for m in lmax]

        scatter!(ax_x, times, xs, markersize = 5 + (length(all_maxima) - idx) * 2)
        scatter!(ax_y, times, ys, markersize = 5 + (length(all_maxima) - idx) * 2)
    end

    return fig
end

function plot_local_maxima_coordinates(slice_vec::Vector{IsoBECSlice})
    fig = Figure()

    ax_x = Axis(fig[1,1], ylabel = "x [nm]", xlabel = "t [ms]")
    ax_y = Axis(fig[2,1], ylabel = "y [nm]", xlabel = "t [ms]")

    n_frames = length(slice_vec)
    slice = interpolate_slice(slice_vec[begin])
    maxima = find_local_maxima(slice)
    coords = get_coordinates(slice, maxima)

    n_becs = length(coords)
    all_maxima = [LocalMaximaPhysical[] for n in 1:n_becs]

    for (idx, _slice) in enumerate(slice_vec)
        println("Processing slice: ", idx)
        slice = interpolate_slice(_slice; size = (100, 100))
        maxima = find_local_maxima(slice)
        coords = get_coordinates(slice, maxima)

        sorted_cords = sort(coords, by = c -> (c.x, c.y))

        for (bec_idx, coord) in enumerate(sorted_cords)
            push!(all_maxima[bec_idx], coord)
        end 
    end

    times = collect(1:length(slice_vec)) .* STEPS_PER_SAVE * time_au_to_ms(TIME_STEP_AU)
    for (idx, lmax) in enumerate(all_maxima)
        xs = [m.x for m in lmax]
        ys = [m.y for m in lmax]

        scatter!(ax_x, times, xs, markersize = 5 + (length(all_maxima) - idx) * 2)
        scatter!(ax_y, times, ys, markersize = 5 + (length(all_maxima) - idx) * 2)
    end

    return fig
end

function plot_local_maxima_coordinates_one_ax(slice_vec::Vector{IsoBECSlice}; leg_pos = :rc)
    fig = Figure()

    ax = Axis(fig[1,1],
              ylabel = L"\text{r}$_{\text{max}}$ [nm]",
              xlabel = "t [ms]",
              yautolimitmargin = (0.01f0, 0.01f0),
              xtickalign = 1,
              ytickalign = 1,
              )

    n_frames = length(slice_vec)
    slice = interpolate_slice(slice_vec[begin])
    maxima = find_local_maxima(slice)
    coords = get_coordinates(slice, maxima)

    n_becs = length(coords)
    all_maxima = [LocalMaximaPhysical[] for n in 1:n_becs]

    for (idx, _slice) in enumerate(slice_vec)
        println("Processing slice: ", idx)
        #slice = interpolate_slice(_slice)
        maxima = find_local_maxima(_slice)
        coords = get_coordinates(slice, maxima)

        sorted_cords = sort(coords, by = c -> (c.x, c.y))

        for (bec_idx, coord) in enumerate(sorted_cords)
            push!(all_maxima[bec_idx], coord)
        end 
    end

    times = collect(1:length(slice_vec)) .* STEPS_PER_SAVE * time_au_to_ms(TIME_STEP_AU)
    for (idx, lmax) in enumerate(all_maxima)
        xs = [m.x for m in lmax]
        ys = [m.y for m in lmax]

        scatter!(ax, times, xs, markersize = 5, color = :red, label = "x")
        scatter!(ax, times, ys, markersize = 5, color = :blue, label = "y")

    end
    
    axislegend(position = leg_pos, unique = true)

    return fig
end


function plot_evolution(data_dir::String; step = 40, total_size = 6)
    files = filter(f -> occursin("wavefunction_", f) && endswith(f, ".gpe.dat"), readdir(data_dir))
    files = sort(files, by = f -> parse(Int, split(split(f, "_")[2], ".")[1]))

    files = files[begin:step:end]
    if(length(files) > total_size)
        files = files[begin:total_size]
    end

    n_frames = length(files)

    ncols = ceil(Int, 2)
    nrows = ceil(Int, 3)
    fig = Figure(size = (728, 1024), figure_padding  = 40)

    for (file_idx, file) in enumerate(files)
        # Load file 
        BCEContext = load_from_text(joinpath(data_dir, file))
        BCEslice = get_BEC_slice(BCEContext)

        n = [length(BCEContext.x), length(BCEContext.y), length(BCEContext.z)]
        begin_idx = [floor(Int, n[1] / 3), floor(Int, n[2]/ 3), floor(Int, n[3] / 5)]
        end_idx = n .- begin_idx

        slices = (
            begin_idx[1]:end_idx[1], # x-direction
            begin_idx[2]:end_idx[2], # y-direction
            begin_idx[3]:end_idx[3]  # z-direction
        )

        rho = abs.(BCEContext.psi[slices...])
        x = BCEContext.x[begin_idx[1]:end_idx[1]]
        y = BCEContext.y[begin_idx[2]:end_idx[2]]
        z = BCEContext.z[begin_idx[3]:end_idx[3]]

        hm_rho = abs.(BCEslice.psi[slices[1:2]...])
        hm_x = BCEslice.x[begin_idx[1]:end_idx[1]]
        hm_y = BCEslice.y[begin_idx[1]:end_idx[1]]

        # Prepare values for isosurfaces
        max_val = maximum(rho)
        alphas = [0.3, 0.8]
        isovals = [0.4 * max_val, 0.8 * max_val]

        # Find idx of axis in frame
        row = div(file_idx - 1, ncols) + 1
        col = mod(file_idx - 1, ncols) + 1

        time_ms = (file_idx - 1) * step * STEPS_PER_SAVE * time_au_to_ms(TIME_STEP_AU)
        letter = 'a' + (file_idx - 1)

        ax          = Axis3(fig[row,col],
                            aspect=:data,
                            protrusions = 0.,
                            xautolimitmargin = (0., 0.),
                            yautolimitmargin = (0., 0.),
                            zautolimitmargin = (0., 0.),
                            xlabel = "x [nm]",
                            ylabel = "y [nm]",
                            zlabel = "z [nm]",
                            xticklabelpad = 1.,
                            yticklabelpad = 1.,
                            zticklabelpad = 1.,
                            xlabeloffset = 25.,
                            ylabeloffset = 25.,
                            zlabeloffset = 25.,
                            viewmode = :stretch)

        set_lights!(ax, [DirectionalLight(RGBf(1,1,1), Vec3f(-1,1,-1))])

        # Plot values
        for (i, isovalue) in enumerate(isovals)
            volume!(ax,
                    (x[begin], x[end]),
                    (y[begin], y[end]),
                    (z[begin], z[end]),
                    rho,
                    algorithm = :iso,
                    isovalue = isovalue,
                    alpha = alphas[i],
                    colormap = :inferno,
                    transparency = true,
                    isorange = 0.1 * max_val)
        end


        heatmap!(ax, BCEslice.x, BCEslice.y, abs.(BCEslice.psi), transparency = true, colormap = :inferno, transformation=(:xy, z[begin]))

        if(file_idx != 1)
            hidedecorations!(ax)
            hidespines!(ax)
        end
    end

    for col in 1:ncols
        colsize!(fig.layout, col, Auto(1/ncols))
    end

    return fig
end

function plot_states(files::Vector{String}; total_size = 6)
    if(length(files) > total_size)
        files = files[begin:total_size]
    end

    n_frames = length(files)

    ncols = ceil(Int, 2)
    nrows = ceil(Int, 3)
    fig = Figure(size = (800, 1200))

    for (file_idx, file) in enumerate(files)
        BCEContext = load_from_text(file)
        BCEslice = get_BEC_slice(BCEContext)

        n = [length(BCEContext.x), length(BCEContext.y), length(BCEContext.z)]
        begin_idx = [floor(Int, n[1] / 4), floor(Int, n[2]/ 4), floor(Int, n[3] / 5)]
        end_idx = n .- begin_idx

        N = match(r"(\d+)k_atoms", file)

        slices = (
            begin_idx[1]:end_idx[1], # x-direction
            begin_idx[2]:end_idx[2], # y-direction
            begin_idx[3]:end_idx[3]  # z-direction
        )

        rho = abs.(BCEContext.psi[slices...])
        x = BCEContext.x[begin_idx[1]:end_idx[1]]
        y = BCEContext.y[begin_idx[2]:end_idx[2]]
        z = BCEContext.z[begin_idx[3]:end_idx[3]]

        # Prepare values for isosurfaces
        max_val = maximum(rho)
        alphas = [0.3, 0.8]
        isovals = [0.4 * max_val, 0.8 * max_val]

        # Find idx of axis in frame
        row = div(file_idx - 1, ncols) + 1
        col = mod(file_idx - 1, ncols) + 1

        letter = 'a' + (file_idx - 1)

        ax          = Axis3(fig[row,col],
                            aspect=:data,
                            protrusions = 0.,
                            viewmode = :fitzoom)

        set_lights!(ax, [DirectionalLight(RGBf(1,1,1), Vec3f(-1,1,-1))])

        # Plot values
        for (i, isovalue) in enumerate(isovals)
            volume!(ax,
                    (x[begin], x[end]),
                    (y[begin], y[end]),
                    (z[begin], z[end]),
                    rho,
                    algorithm = :iso,
                    isovalue = isovalue,
                    alpha = alphas[i],
                    colormap = :inferno,
                    transparency = true,
                    isorange = 0.1 * max_val)
        end
        heatmap!(ax, BCEslice.x, BCEslice.y, abs.(BCEslice.psi), transparency = true, colormap = :inferno, transformation=(:xy, z[begin]))

        hidedecorations!(ax)
        hidespines!(ax)
    end

    for col in 1:ncols
        colsize!(fig.layout, col, Auto(1/ncols))
    end

    return fig
end

function plot_single_state(file::String; hide_decs = true)
    BCEContext = load_from_text(file)
    BCEslice = get_BEC_slice(BCEContext)

    n = [length(BCEContext.x), length(BCEContext.y), length(BCEContext.z)]
    begin_idx = [floor(Int, n[1] / 4), floor(Int, n[2]/ 4), floor(Int, n[3] / 4.5)]
    end_idx = n .- begin_idx

    slices = (
        begin_idx[1]:end_idx[1], # x-direction
        begin_idx[2]:end_idx[2], # y-direction
        begin_idx[3]:end_idx[3]  # z-direction
    )

    rho = abs.(BCEContext.psi[slices...])
    x = BCEContext.x[begin_idx[1]:end_idx[1]]
    y = BCEContext.y[begin_idx[2]:end_idx[2]]
    z = BCEContext.z[begin_idx[3]:end_idx[3]]
    
    hm_rho = abs.(BCEslice.psi[slices[1:2]...])
    hm_x = BCEslice.x[begin_idx[1]:end_idx[1]]
    hm_y = BCEslice.y[begin_idx[1]:end_idx[1]]

    # Prepare values for isosurfaces
    max_val = maximum(rho)
    alphas = [0.3, 0.8]
    isovals = [0.4 * max_val, 0.8 * max_val]

    fig = Figure(figure_padding = 0)

    ax          = Axis3(fig[1,1],
                        protrusions = 0,
                        aspect=:data,
                        xautolimitmargin = (0., 0.),
                        yautolimitmargin = (0., 0.),
                        zautolimitmargin = (0., 0.),
                        xlabeloffset = 25.,
                        ylabeloffset = 25.,
                        zlabeloffset = 30.,
                        zlabelrotation = 0.,
                        zlabelsize = 20.,
                        zticklabelsize = 20.,
                        xlabel = "x [nm]",
                        ylabel = "y [nm]",
                        zlabel = "z [nm]",
                        viewmode = :fit)

    set_lights!(ax, [DirectionalLight(RGBf(1,1,1), Vec3f(-1,1,-1))])

    colormap = :inferno
    for (i, isovalue) in enumerate(isovals)
        volume!(ax,
                (x[begin], x[end]),
                (y[begin], y[end]),
                (z[begin], z[end]),
                rho,
                algorithm = :iso,
                isovalue = isovalue,
                alpha = alphas[i],
                colormap = colormap,
                transparency = true,
                isorange = 0.1 * max_val)
    end
    heatmap!(ax, hm_x, hm_y, hm_rho, transparency = true, colormap = colormap, transformation=(:xy, z[begin + 2]))

    origin = Point3f(x[begin], y[end], z[begin + 2])
    direction = Vec3f(0, 0, abs(z[end] - z[begin + 2]))

    if(hide_decs)
        hidedecorations!(ax)
        hidespines!(ax)
    else
        arrows2d!(ax, [origin], [direction]; 
                  color=:black
                  )
        hidespines!(ax)
        hidexdecorations!(ax)
        hideydecorations!(ax)
        #hidexspines!(ax)
        #hideyspines!(ax)
    end

    ax.zticks = [-5000, -2500, 2500,5000]

    fig
end
