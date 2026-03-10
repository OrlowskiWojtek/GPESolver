using DataFrames
using CSV
using Glob

include("context.jl")
include("files.jl")
include("units.jl")

function roman_from_idx(idx::Int)
    romans = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X"]
    if 1 ≤ idx ≤ length(romans)
        return romans[idx]
    else
        error("Index out of range (0–9)")
    end
end

# another approach; - load all output files, then find number of maximas assiociated to energy
function gather_energy(dist_dir::String; BCE_THRESHOLD = nothing)
    df = DataFrame(atom_number=Int[],
                   max_number=Int[],
                   e_kin=Float64[],
                   e_pot=Float64[],
                   e_int=Float64[],
                   e_ext=Float64[],
                   e_bmf=Float64[],
                   e_total=Float64[],
                   filename=String[])

    for atom_folder in 1:54
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
                continue
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
                slice = load_slice_from_text(psi_file)
                l_maxes = number_of_lmax(slice ;n_atoms = atom_number, condensation_threshold = 1e-9)

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
                    parse(Float64, numbers[7]),
                    psi_file))
            end
        end
    end

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

    conv = 27211.4 * 11.6 * 10^9
    labels = ["Total", "Kinetic", "Potential", "Interaction"]
    fields = [:e_total, :e_kin, :e_pot, :e_tot_ext]
    colors = [:black, :orange, :red, :blue]

    for (seg_idx, seg) in enumerate(segments)
        for (i, (field, label)) in enumerate(zip(fields, labels))
            ydata = seg[!, field] * conv ./ (seg.atom_number * 10^3)
            xdata = seg.atom_number

            lines!(
                ax,
                xdata,
                ydata,
                color = colors[i],
            )

            if(seg_idx == length(segments))
                text_x_pos = 45
                text!(
                    ax,
                    text_x_pos,
                    ydata[end-3],
                    text = labels[i],
                    color = colors[i],
                    align = (:center, :bottom),
                    fontsize = 18)
            end
        end 
    end

    for (idx, seg) in enumerate(segments[begin:end-1])
        change_point = seg.atom_number[end] + 0.5
        vlines!(ax, [change_point], color=:red, linestyle=:dash, linewidth=2)
    end

    for (idx, seg) in enumerate(segments[begin:end])
        text!(ax,
              seg.atom_number[begin] + (seg.atom_number[end] - seg.atom_number[begin]) / 2.,
              18.,
              text = roman_from_idx(idx),
              align = (:center, :bottom),
              fontsize = 25)
    end

    ax.xticks = 0:5:50

    return fig
end

function get_rhomax(segments, data_dist)
    all_rho = Vector{Vector{Vector{BECMaxRho}}}()

    for seg in segments
        seg_rhos = Vector{Vector{BECMaxRho}}()

        for row in eachrow(seg)
            atom_number = row.atom_number
            max_number = row.max_number

            max_folder = joinpath(data_dist, "$(atom_number)k_atoms", "$(max_number)_max")
            psi_file = row.filename
            if !isfile(psi_file)
                @info "No such file $psi_file"
                psi_file = joinpath(max_folder, "initial_state.dat")
            end
            if !isfile(psi_file)
                @info "No wavefunction file for atom_number=$atom_number, max_number=$max_number"
                continue
            end
            context = load_from_text(psi_file)
            rhos = get_BEC_maxrho(context)
            
            bec_rhos = [BECMaxRho(rho, atom_number) for rho in rhos]
            push!(seg_rhos, bec_rhos)
        end

        push!(all_rho, seg_rhos)
    end

    return all_rho
end

function get_heights(segments, data_dist)
    all_heights = Vector{Vector{Vector{BECHeight}}}()

    for seg in segments
        seg_heights = Vector{Vector{BECHeight}}()

        for row in eachrow(seg)
            atom_number = row.atom_number
            max_number = row.max_number
            max_folder = joinpath(data_dist, "$(atom_number)k_atoms", "$(max_number)_max")
            psi_file = row.filename

            if !isfile(psi_file)
                @info "No wavefunction file for atom_number=$atom_number, max_number=$max_number"
                continue
            end

            context = load_from_text(psi_file)
            heights = get_BEC_heights(context, atom_number)

            push!(seg_heights, heights)
        end

        push!(all_heights, seg_heights)
    end

    return all_heights
end

function plot_heights(segments, seg_heights)
    fig = Figure(size = (600,600))
    ax = Axis(fig[1, 1],
              xlabel="N/10⁴",
              ylabel="H (μm)",
              xlabelsize = 16,
              ylabelsize = 16,
              xticksize = 14,
              yticksize = 14,
              xtickalign = 1,
              ytickalign = 1
              )

    for (seg_idx, (seg, heights)) in enumerate(zip(segments, seg_heights))
        if isempty(heights)
            continue
        end

        xs = seg.atom_number / 10
        n_becs = length(heights[begin])

        plot_heights = [Float64[] for _i in 1:n_becs]
        for bec_idx in 1:n_becs
            for h in heights
                push!(plot_heights[bec_idx], length_au_to_μm(h[bec_idx].height))
            end
        end 
        for (idx, ph) in enumerate(plot_heights)
            lines!(ax, xs, ph, color = :black)
        end
    end

    for (idx, seg) in enumerate(segments[begin:end-1])
        change_point = seg.atom_number[end] + 0.5
        vlines!(ax, [change_point / 10], color=:red, linestyle=:dash, linewidth=2)
    end

    return fig
end

function plot_rhomax(segments, seg_rhos)
    fig = Figure(size = (600,600))
    ax = Axis(fig[1, 1],
              xlabel="N/10⁴",
              ylabel="ρ_{max}",
              xlabelsize = 16,
              ylabelsize = 16,
              xticksize = 14,
              yticksize = 14,
              xtickalign = 1,
              ytickalign = 1
              )

    for (seg_idx, (seg, rhos)) in enumerate(zip(segments, seg_rhos))
        if isempty(rhos)
            continue
        end

        xs = seg.atom_number / 10
        n_becs = length(rhos[begin])

        plot_rhos = [Float64[] for _i in 1:n_becs]
        for bec_idx in 1:n_becs
            for r in rhos
                push!(plot_rhos[bec_idx], r[bec_idx].value * r[bec_idx].n_atoms)
            end
        end 

        for (idx, pr) in enumerate(plot_rhos)
            lines!(ax, xs, pr, color = :black)
        end
    end

    for (idx, seg) in enumerate(segments[begin:end-1])
        change_point = seg.atom_number[end] + 0.5
        vlines!(ax, [change_point / 10], color=:red, linestyle=:dash, linewidth=2)
    end

    return fig
end

function save_segments(segments, filename::String)
    open(filename, "w") do f
        for (seg_idx, seg) in enumerate(segments)
            for row in eachrow(seg)
                write(f, "$(row.atom_number)\t$(row.e_total)\t$(row.e_kin)\t$(row.e_pot)\t$(row.e_int)\t$(row.e_ext)\t$(row.e_bmf)\n")
            end
            # Empty line separates segments in gnuplot format
            if seg_idx < length(segments)
                write(f, "\n")
            end
        end
    end
end

function save_rhomax(seg_rhos::Vector{Vector{Vector{BECMaxRho}}}, filename::String)
    open(filename, "w") do f
        for (seg_idx, seg) in enumerate(seg_rhos)
            n_becs = length(seg[1])  # Get number of BECs from first row
            for bec_idx in 1:n_becs
                for rhos in seg
                    rho = rhos[bec_idx]
                    write(f, "$(rho.n_atoms)\t$(rho.value * rho.n_atoms)\t$bec_idx\n")
                end
                # Empty line separates different BECs within a segment
                if bec_idx < n_becs
                    write(f, "\n")
                end
            end
            # Empty line separates segments in gnuplot format
            if seg_idx < length(seg_rhos)
                write(f, "\n")
            end
        end
    end
end

function save_heights(seg_heights::Vector{Vector{Vector{BECHeight}}}, segments::Vector{DataFrame}, filename::String)
    open(filename, "w") do f
        for (seg_idx, (seg, seg_heights_data)) in enumerate(zip(segments, seg_heights))
            n_becs = length(seg_heights_data[1])  # Get number of BECs from first row
            for bec_idx in 1:n_becs
                for (row_idx, heights) in enumerate(seg_heights_data)
                    atom_number = seg[row_idx, :atom_number]
                    height = heights[bec_idx]
                    write(f, "$(atom_number)\t$(height.height)\t$(height.value)\t$bec_idx\n")
                end
                # Empty line separates different BECs within a segment
                if bec_idx < n_becs
                    write(f, "\n")
                end
            end
            # Empty line separates segments in gnuplot format
            if seg_idx < length(seg_heights)
                write(f, "\n")
            end
        end
    end
end

function plot_df(df)
    fig = Figure();
    ax = Axis(fig[1,1], xlabel = "N / 10^3", ylabel = "E/N [nK]")

    temp = ["4max", "5max"]
    conv = 27211.4 * 11.6 * 10^9
    for (idx, nmdf) in enumerate(groupby(df, :max_number))
        if(idx < 4)
            continue
        end
        lines!(ax, nmdf[!, :atom_number], nmdf[!, :e_total] .* conv ./ (nmdf[!, :atom_number] * 10^3), label = "$(temp[idx - 3])")
    end

    axislegend()
    xlims!(ax, (46, 54))
    ylims!(ax, (35, 38))

    fig
end

##

df = gather_energy("../../../../data/run_find_initial_states_eps145")
segments = segmentize_dataframe(df)

##

fig = plot_df(df)
save("4vs5_ene_diff.png", fig)

##

fig = plot_segments(segments)
#save("eps_15_stable_phases.pdf", fig)
save_segments(segments, "eps_145_segments.dat")

##

seg_heights = get_heights(segments, "../../../../data/run_find_initial_states_eps_145")

##

seg_rhos = get_rhomax(segments, "../../../../data/run_find_initial_states_eps_145")

##

fig = plot_heights(segments, seg_heights)
#save("eps_145_heights.pdf", fig)
save_heights(seg_heights, segments, "eps_145_height.dat")

##

fig = plot_rhomax(segments, seg_rhos)
save_rhomax(seg_rhos, "eps_145_maxrho.dat")

##

## Plot in order to find threshhold

dist_dir = "../../../../data/run_find_initial_states_single_well"
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
            @info "No wavefunction file"
            continue
        end

        slice = load_slice_from_text(psi_file)
        n_lmax = number_of_lmax(slice)      

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
