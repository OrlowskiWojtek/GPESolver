using DataFrames
using CSV
using Glob

include("context.jl")

# another approach; - load all output files, then find number of maximas assiociated to energy
function gather_energy(dist_dir::String; BCE_THRESHOLD = nothing)
    df = DataFrame(atom_number=Int[],
                   max_number=Int[],
                   e_kin=Float64[],
                   e_pot=Float64[],
                   e_int=Float64[],
                   e_ext=Float64[],
                   e_bmf=Float64[],
                   e_total=Float64[])

    for atom_folder in 1:50
        atom_dir = joinpath(dist_dir, "$(atom_folder)k_atoms")
        if !isdir(atom_dir)
            continue
        end

        max_folders = glob("*_max", atom_dir)
        for max_folder in max_folders
            energy_file = joinpath(max_folder, "energy.gpe.dat")

            if !isfile(energy_file)
                energy_file = joinpath(max_folder, "energy.dat")
            end

            if !isfile(energy_file)
                continue
            end

            psi_file = joinpath(max_folder, "initial_state.gpe.dat")

            if !isfile(psi_file)
                psi_file = joinpath(max_folder, "initial_state.dat")
            end

            if !isfile(psi_file)
                @error "No wavefunction file"
            end

            open(energy_file, "r") do f
                lines = readlines(f)
                line = lines[end]
                line = strip(line)
                if isempty(line)
                    return
                end
                numbers = split(line, '\t')
                if length(numbers) < 7
                    return
                end
                atom_number = parse(Int, match(r"(\d+)k_atoms", atom_dir).captures[1])
                max_number = parse(Int, match(r"(\d+)_max", basename(max_folder)).captures[1])
                psi = load_from_text(psi_file)
                l_maxes = number_of_lmax(psi;n_atoms = atom_number, condensation_threshold = 1e-6)

                if(l_maxes != max_number)
                    @info "Number of condensates changed from $max_number to $l_maxes in $atom_number k atoms"
                    max_number = l_maxes
                end

                push!(df, (atom_number,
                    max_number,
                    parse(Float64, numbers[2]),
                    parse(Float64, numbers[3]),
                    parse(Float64, numbers[4]),
                    parse(Float64, numbers[5]),
                    parse(Float64, numbers[6]),
                    parse(Float64, numbers[7])))
            end
        end
    end

    #CSV.write("energies.txt", df, delim=' ')
    return df
end

function segmentize_dataframe(df)
    grouped = groupby(df, :atom_number)
    df.e_tot_ext = df.e_int .+ df.e_ext .+ df.e_bmf

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

function plot_segments(segments)
    fig = Figure();
    ax = Axis(fig[1, 1],
              xlabel="N/10³",
              ylabel="E/N (nK)")

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

    for seg in segments[begin:end-1]
        change_point = seg.atom_number[end] + 0.5
        vlines!(ax, [change_point], color=:red, linestyle=:dash, linewidth=2)
    end

    Legend(fig[2, 1], ax, orientation=:horizontal, tellwidth=false)

    return fig
end

##

df = gather_energy("../../../../data/run_find_initial_states")
segments = segmentize_dataframe(df)

##

fig = plot_segments(segments)
#save("eps_15_stable_phases.pdf", fig)

##


## Plot in order to find threshhold

dist_dir = "../../../../data/run_find_initial_states"
maxs = Float64[]
maxs_atoms = Float64[]
maxs_atoms_per_bec = Float64[]
atoms = Float64[]
for atom_folder in 1:50
    atom_dir = joinpath(dist_dir, "$(atom_folder)k_atoms")
    if !isdir(atom_dir)
        continue
    end

    max_folders = glob("*_max", atom_dir)
    for max_folder in max_folders
        psi_file = joinpath(max_folder, "initial_state.gpe.dat")

        if !isfile(psi_file)
            psi_file = joinpath(max_folder, "initial_state.dat")
        end

        if !isfile(psi_file)
            @error "No wavefunction file"
        end

        psi = load_from_text(psi_file)
        slice = get_BEC_slice(psi)
        n_lmax = number_of_lmax(psi)      

        push!(maxs, maximum(abs.(slice.psi)))
        push!(maxs_atoms, maximum(abs.(slice.psi)) * atom_folder)
        push!(maxs_atoms_per_bec, maximum(abs.(slice.psi)) * atom_folder / n_lmax)
        push!(atoms, atom_folder)
    end
end

##
GLMakie.activate!()

fig = Figure();
ax = Axis(fig[1,1]);
scatter!(ax, atoms, maxs, label = "max", markersize = 10);
scatter!(ax, atoms, maxs_atoms, label = "max atoms", markersize = 8);
scatter!(ax, atoms, maxs_atoms_per_bec, label = "max atoms / bec", markersize = 6);
#vlines!(ax, [6], color = :black)

axislegend()
display(fig)

##
