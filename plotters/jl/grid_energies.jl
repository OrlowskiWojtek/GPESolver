using GLMakie, CairoMakie, LaTeXStrings

include("files.jl")
include("plots.jl")
include("animate.jl")

##

const BASE_DIR = "../../build/run_check_size"

"""
    find_grid_sizes(results::Dict)

Extracts unique and sorted xy and z sizes from results
"""
function find_grid_sizes(results::Dict)
    xy_sizes = unique(first.(keys(results)))
    z_sizes = unique(last.(keys(results)))
    sort!(xy_sizes)
    sort!(z_sizes)
    return xy_sizes, z_sizes
end

"""
    load_all_grid_energies(base_dir::String)

Recursively searches for energy.gpe.dat files in directories like:
BASE_DIR/128xy/32z/energy.gpe.dat
Returns dictionaries mapping (xy_size, z_size) -> NamedTuple of energies
"""
function load_all_grid_energies(base_dir::String)
    results = Dict{Tuple{Int, Int}, NamedTuple{(:e_kin, :e_pot, :e_int, :e_ext, :e_bmf, :e_tot), Tuple{Float64, Float64, Float64, Float64, Float64, Float64}}}()
    
    # Find all xy directories (e.g., "128xy", "256xy")
    xy_pattern = r"^(\d+)xy$"
    z_pattern = r"^(\d+)z$"
    
    for xy_entry in readdir(base_dir, join=true)
        isdir(xy_entry) || continue
        
        xy_match = match(xy_pattern, basename(xy_entry))
        xy_match === nothing && continue
        
        xy_size = parse(Int, xy_match.captures[1])
        
        # Find all z subdirectories
        for z_entry in readdir(xy_entry, join=true)
            isdir(z_entry) || continue
            
            z_match = match(z_pattern, basename(z_entry))
            z_match === nothing && continue
            
            z_size = parse(Int, z_match.captures[1])
            
            energy_file = joinpath(z_entry, "energy.gpe.dat")
            if isfile(energy_file)
                energies = load_energies(energy_file)
                n = length(energies.e_kin)
                results[(xy_size, z_size)] = (
                    e_kin = energies.e_kin[n],
                    e_pot = energies.e_pot[n],
                    e_int = energies.e_int[n],
                    e_ext = energies.e_ext[n],
                    e_bmf = energies.e_bmf[n],
                    e_tot = energies.e_tot[n]
                )
                println("Loaded: $(xy_size)xy, $(z_size)z -> E_tot = $(energies.e_tot[n])")
            end
        end
    end
    
    return results
end

# Energy conversion constant (same as in phase diagram)
const ENERGY_CONV = 27211.4 * 11.6 * 1e9 / 20_000

# Available energy types
const ENERGY_TYPES = (:e_tot, :e_kin, :e_pot, :e_int, :e_ext, :e_bmf)
const ENERGY_LABELS = Dict(
    :e_tot => L"E_{\mathrm{total}}",
    :e_kin => L"E_{\mathrm{kin}}",
    :e_pot => L"E_{\mathrm{pot}}",
    :e_int => L"E_{\mathrm{int}}",
    :e_ext => L"E_{\mathrm{ext}}",
    :e_bmf => L"E_{\mathrm{bmf}}"
)

"""
    build_energy_matrix(results, xy_sizes, z_sizes; energy_type=:e_tot)

Creates a matrix of energies indexed by (xy_size, z_size)
"""
function build_energy_matrix(results, xy_sizes, z_sizes; energy_type::Symbol=:e_tot)
    n_xy = length(xy_sizes)
    n_z = length(z_sizes)
    energy_matrix = fill(NaN, n_xy, n_z)
    
    for (i, xy) in enumerate(xy_sizes)
        for (j, z) in enumerate(z_sizes)
            if haskey(results, (xy, z))
                energy_matrix[i, j] = results[(xy, z)][energy_type]
            end
        end
    end
    
    return energy_matrix
end

"""
    visualize_grid_convergence(base_dir::String; output_path::Union{Nothing, String}=nothing, energy_type=:e_tot)

Main function to visualize energy convergence with grid size
- energy_type: choose from :e_tot, :e_kin, :e_pot, :e_int, :e_ext, :e_bmf
"""
function visualize_grid_convergence(base_dir::String; output_path::Union{Nothing, String}=nothing, energy_type::Symbol=:e_tot)
    energy_type in ENERGY_TYPES || error("Unknown energy type: $energy_type. Choose from $ENERGY_TYPES")
    
    CairoMakie.activate!()
    results = load_all_grid_energies(base_dir)
    
    if isempty(results)
        error("No energy data found in $base_dir")
    end
    
    xy_sizes, z_sizes = find_grid_sizes(results)
    energy_matrix = build_energy_matrix(results, xy_sizes, z_sizes; energy_type=energy_type)
    
    # Create figure
    fig = Figure()

    ax = Axis(fig[1, 1],
              xlabel="n",
              ylabel=" E/N [nK]",
              xticklabelsize=18,
              yticklabelsize=18,
              xlabelsize=20,
              ylabelsize=20,
              xticks=xy_sizes,
              xtickalign=0,
              ytickalign=0,
              )
    
    for (j, z) in enumerate(z_sizes)
        energies_for_z = Float64[]
        sizes_for_z = Int[]
        for (i, xy) in enumerate(xy_sizes)
            if !isnan(energy_matrix[i, j])
                push!(sizes_for_z, xy)
                push!(energies_for_z, energy_matrix[i, j])
            end
        end
        if !isempty(energies_for_z)
            lines!(ax, sizes_for_z, energies_for_z * ENERGY_CONV,
                label="$(z)",
                linestyle = j % 2 == 1 ? :dash : :solid,
                linewidth = 3
            )
        end
    end
    
    axislegend(ax, L"\mathrm{n}_z", position=(:right, :center), titlesize=18, labelsize=16)
    
    if output_path !== nothing
        save(output_path, fig, pt_per_unit=1)
    end
    
    return fig
end

function visualize_all_energies(base_dir::String; output_path::Union{Nothing, String}=nothing)
    CairoMakie.activate!()
    results = load_all_grid_energies(base_dir)
    
    if isempty(results)
        error("No energy data found in $base_dir")
    end
    
    xy_sizes, z_sizes = find_grid_sizes(results)
    
    # Create figure with 2x3 subplots
    fig = Figure(size=(1200, 800))
    
    for (idx, energy_type) in enumerate(ENERGY_TYPES)
        energy_matrix = build_energy_matrix(results, xy_sizes, z_sizes; energy_type=energy_type)
        
        # Calculate row and column for subplot
        row = (idx - 1) ÷ 3 + 1
        col = (idx - 1) % 3 + 1
        
        ax = Axis(fig[row, col],
                  xlabel=L"\mathrm{n}_{xy}",
                  ylabel="E/N [nK]",
                  xticklabelsize=12,
                  yticklabelsize=12,
                  xlabelsize=14,
                  ylabelsize=14,
                  xticks=xy_sizes,
                  title=String(energy_type),
                  titlesize=16
                  )
        
        for (j, z) in enumerate(z_sizes)
            energies_for_z = Float64[]
            sizes_for_z = Int[]
            for (i, xy) in enumerate(xy_sizes)
                if !isnan(energy_matrix[i, j])
                    push!(sizes_for_z, xy)
                    push!(energies_for_z, energy_matrix[i, j])
                end
            end
            if !isempty(energies_for_z)
                scatter!(ax, sizes_for_z, energies_for_z * ENERGY_CONV,
                    label="$(z)",
                )
            end
        end
        
        axislegend(ax, L"\mathrm{n}_z", position=(:right, :center), titlesize=12, labelsize=10)
    end
    
    if output_path !== nothing
        save(output_path, fig, pt_per_unit=1)
    end
    
    return fig
end

# Convenience functions for each energy type
visualize_grid_total(args...; kwargs...) = visualize_grid_convergence(args...; energy_type=:e_tot, kwargs...)
visualize_grid_kinetic(args...; kwargs...) = visualize_grid_convergence(args...; energy_type=:e_kin, kwargs...)
visualize_grid_potential(args...; kwargs...) = visualize_grid_convergence(args...; energy_type=:e_pot, kwargs...)
visualize_grid_interaction(args...; kwargs...) = visualize_grid_convergence(args...; energy_type=:e_int, kwargs...)
visualize_grid_external(args...; kwargs...) = visualize_grid_convergence(args...; energy_type=:e_ext, kwargs...)
visualize_grid_bmf(args...; kwargs...) = visualize_grid_convergence(args...; energy_type=:e_bmf, kwargs...)

## Run visualization

fig = visualize_grid_total(BASE_DIR);
save("grid.pdf", fig)

##
