using GLMakie

include("context.jl")
include("files.jl")


function animate_iso_bce(output_file::String)
    files = filter(f -> occursin("checkpoint_", f) && endswith(f, ".bin"), readdir(DATA_DIR))
    files = sort(files, by = f -> parse(Int, split(split(f, "_")[2], ".")[1]))
    n_frames = length(files)

    fig = Figure()
    ax = Axis3(fig[1, 1], xlabel="X", ylabel="Y", zlabel="Z", title="Isosurface of BCE Structure", aspect=:data)

    BCEContext = load_from_binary(joinpath(DATA_DIR, files[1]))
    rho = Observable(Array{Float64,3}(abs.(BCEContext.psi)))
    x, y, z = BCEContext.x, BCEContext.y, BCEContext.z

    max_val = maximum(rho[])
    isovals = [0.2 * max_val, 0.5 * max_val, 0.8 * max_val]

    for (i, isovalue) in enumerate(isovals)
        volume!(ax,
                (x[begin], x[end]),
                (y[begin], y[end]),
                (z[begin], z[end]),
                rho,
                algorithm = :iso,
                isovalue = isovalue,
                alpha = 0.8,
                colormap = :plasma,
                transparency = true,
                isorange = 0.1 * max_val)
    end

    record(fig, output_file, 1:n_frames; framerate=10) do frame
        println("loading frame", joinpath(DATA_DIR, files[frame]))
        rho[] = abs.(load_from_binary(joinpath(DATA_DIR, files[frame])).psi)
    end
end
