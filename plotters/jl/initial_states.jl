using DataFrames
using CSV
using CairoMakie

DATA_FILE = "energies.txt"

#= LOADS DATAFRAME AND SELECTS SEGMENTS FOR LOWEST ENERGY
#
#
=#
function load_energy_dataframe(filename)  
    df = CSV.read(filename, DataFrame; delim=' ')
    df.e_tot_ext = df.e_int .+ df.e_ext .+ df.e_bmf

    grouped = groupby(df, :atom_number)
    lowest_energy = DataFrame()
    for g in grouped
        idx = argmin(g.e_total)
        push!(lowest_energy, g[idx, :])
    end

    println("Loaded df: ", df)

    segments = Vector{DataFrame}()
    start_idx = 1
    for i in 2:nrow(lowest_energy)
        if lowest_energy.max_number[i] != lowest_energy.max_number[i-1]
            push!(segments, lowest_energy[start_idx:i-1, :])
            start_idx = i
        end
    end
    push!(segments, lowest_energy[start_idx:end, :])  # Add the last segment

    return segments;
end

##

function plot_lowest_energy_segments(segments)
    CairoMakie.activate!()

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Number of atoms/10³", ylabel="Energy per atom (nK)", title="Energies for Lowest Energy State")

    conversion = 27211.4 * 11.6 * 10^9
    labels = ["Total", "Kinetic", "Potential", "Interaction"]
    fields = [:e_total, :e_kin, :e_pot, :e_tot_ext]
    colors = [:black, :orange, :red, :blue]
    legend_drawn = fill(false, length(fields))

    for segment in segments
        for (i, (field, label)) in enumerate(zip(fields, labels))
            lines!(ax,
                   segment.atom_number,
                   segment[!, field] .* conversion ./ (segment.atom_number .* 10^3),
                   color = colors[i],
                   label = legend_drawn[i] ? nothing : label)
            legend_drawn[i] = true
        end
    end

    for segment in segments
        change_point = segment.atom_number[end] + 0.5
        vlines!(ax, [change_point], color=:red, linestyle=:dash, linewidth=2)
    end

    Legend(fig[2, 1], ax, orientation=:horizontal, tellwidth=false)
    return fig
end

segments = load_energy_dataframe(DATA_FILE)
fig = plot_lowest_energy_segments(segments)
#save("lowest_energy_all_components.pdf", fig)
