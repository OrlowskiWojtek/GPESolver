using CairoMakie

include("files.jl")
include("plots.jl")
include("animate.jl")
include("units.jl")

JOINING_DATA_DIR = "../../../data/run_zlepanie/"

struct JoiningRunData
    n_atoms::Int32
    slices::Vector{IsoBECSlice}
end

function plot_x_joining(j_data::Vector{JoiningRunData})
    fig = Figure()
    ax_x = Axis(fig[1,1], ylabel = "x [nm]", xlabel = "t [ms]")

    colors = cgrad(:berlin, length(j_data) ; categorical = true)

    for (run_idx, run) in enumerate(j_data)
        println("Plotting run: ", run_idx)

        for (idx, _slice) in enumerate(run.slices)
            time = idx * STEPS_PER_SAVE * time_au_to_ms(TIME_STEP_AU)

            slice = interpolate_slice(_slice)
            maxima = find_local_maxima(slice)
            coords = get_coordinates(slice, maxima)

            for (bec_idx, coord) in enumerate(coords)
                scatter!(ax_x, [time], [coord.x], color = colors[run_idx], label  = "$(run.n_atoms * 100)")
            end 
        end
    end

    #axislegend("Number of atoms", unique = true)
    Colorbar(fig[1,2]; limits = (j_data[begin].n_atoms, j_data[end].n_atoms), colormap = colors, label = "N / 10²")

    return fig
end

##

j_data = Vector{JoiningRunData}(undef, 0)

for natom in 80:2:98 
    slice_vec = load_directory_slice_from_text(joinpath(JOINING_DATA_DIR, "$(natom)k_atoms/2_max"))
    push!(j_data, JoiningRunData(natom, slice_vec))
end

##

#fig = plot_x_joining(j_data[1:10])
#save("zlepanie.pdf", fig)

