using DataFrames
using CSV
using CairoMakie

CairoMakie.activate!()

DATA_FILE = "../../data/run_find_initial_state/energies.txt"
##
df = CSV.read(DATA_FILE, DataFrame; delim=' ')
df.e_tot_ext = df.e_int .+ df.e_ext .+ df.e_bmf
# Find, for each atom_number, the row with the lowest e_total
grouped = groupby(df, :atom_number)
lowest_energy = DataFrame()
for g in grouped
    idx = argmin(g.e_total)
    push!(lowest_energy, g[idx, :])
end

##
segments = Vector{DataFrame}()
start_idx = 1
for i in 2:nrow(lowest_energy)
    if lowest_energy.max_number[i] != lowest_energy.max_number[i-1]
        push!(segments, lowest_energy[start_idx:i-1, :])
        start_idx = i
    end
end
push!(segments, lowest_energy[start_idx:end, :])  # Add the last segment

# Plot all energies for the lowest energy configuration at each atom_number, segment by segment
fig = Figure();
ax = Axis(fig[1, 1], xlabel="Number of atoms/10³", ylabel="Energy per atom (nK)", title="Energies for Lowest Energy State")

conv = 27211.4 * 11.6*10^9
labels = ["Total", "Kinetic", "Potential", "Interaction"]
fields = [:e_total, :e_kin, :e_pot, :e_tot_ext]
colors = [:black, :orange, :red, :blue]
# For legend: only label the first segment for each energy type
legend_drawn = fill(false, length(fields))

for seg in segments
    for (i, (field, label)) in enumerate(zip(fields, labels))
        lines!(
            ax,
            seg.atom_number,
            seg[!, field] * conv ./ (seg.atom_number * 10^3),
            color = colors[i],
            label = legend_drawn[i] ? nothing : label
        )
        legend_drawn[i] = true
    end
end

# Add vertical lines at change points
change_points = Float64[]
prev = nothing
for row in eachrow(lowest_energy)
    if prev !== nothing && row.max_number != prev
        push!(change_points, row.atom_number .- 0.5)
    end
    prev = row.max_number
end

for x in change_points
    vlines!(ax, [x], color=:red, linestyle=:dash, linewidth=2)
end

max_energy = maximum(df.e_total ./ (df.atom_number * 10^3)) * conv
text!(ax, 20, max_energy / 3., text="III", color=:black, fontsize=24)
text!(ax, 32, max_energy / 3., text="IV", color=:black, fontsize=24)
text!(ax, 43, max_energy / 3., text="V", color=:black, fontsize=24)

Legend(fig[2, 1], ax, orientation=:horizontal, tellwidth=false)

fig
save("lowest_energy_all_components.pdf", fig)
