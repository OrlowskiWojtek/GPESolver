using GLMakie
# This files includes sets of plots and animations that require time evolution of BEC

include("context.jl")
include("files.jl")
include("units.jl")

function animate_iso_bce(data_dir::String, output_file::String)
    files = filter(f -> occursin("wavefunction_", f) && endswith(f, ".gpe.dat"), readdir(data_dir))
    files = sort(files, by = f -> parse(Int, split(split(f, "_")[2], ".")[1]))
    n_frames = length(files)

    fig = Figure(size = (1024, 768), backgroundcolor = :gray97)
    ax = Axis3(fig[1, 1],
               xlabel="X [nm]",
               ylabel="Y [nm]",
               zlabel="Z [nm]",
               title="Isosurface of BCE",
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
        slice = interpolate_slice(get_BEC_slice(context))
        maxima = find_local_maxima(slice)
        coords = get_coordinates(slice, maxima)

        for (bec_idx, coord) in enumerate(coords)
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



