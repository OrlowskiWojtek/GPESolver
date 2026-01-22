using GLMakie

include("context.jl")
include("files.jl")

function animate_iso_bce(output_file::String)
    files = filter(f -> occursin("wavefunction_", f) && endswith(f, ".gpe.dat"), readdir(DATA_DIR))
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

    BCEContext = load_from_text(joinpath(DATA_DIR, files[1]))
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
        println("loading frame", joinpath(DATA_DIR, files[frame]))
        rho[] = abs.(load_from_text(joinpath(DATA_DIR, files[frame])).psi)
    end
end

function plot_local_maxima_evolution()
    files = filter(f -> occursin("wavefunction_", f) && endswith(f, ".gpe.dat"), readdir(DATA_DIR))
    files = sort(files, by = f -> parse(Int, split(split(f, "_")[2], ".")[1]))
    n_frames = length(files)

    fig = Figure(size = (1024, 768))
    ax = Axis(fig[1, 1],
              xlabel="X [nm]",
              ylabel="Y [nm]",
              title="Local Maxima Evolution",
    )

    all_maxima = Vector{Tuple{Float64, Float64, Float64, Int}}()

    for (frame_idx, file) in enumerate(files)
        println("loading frame", joinpath(DATA_DIR, file))
        context = load_from_text(joinpath(DATA_DIR, file))
        slice = get_BEC_slice(context)
        maxima = find_local_maxima(slice)
        coords = get_coordinates(slice, maxima)
        for m in coords
            push!(all_maxima, (m.x, m.y, m.value, frame_idx))
        end
    end

    if !isempty(all_maxima)
        xs = [m[1] for m in all_maxima]
        ys = [m[2] for m in all_maxima]
        frames = [m[4] for m in all_maxima]

        scatter!(ax, xs, ys;
                 color=frames,
                 colormap=:viridis,
                 markersize=8,
                 strokewidth=0)
        
        Colorbar(fig[1, 2], label="Time Step")
    end

    return fig
end
