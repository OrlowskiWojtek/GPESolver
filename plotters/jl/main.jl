using GLMakie

include("files.jl")
include("plots.jl")
include("animate.jl")

##

GLMakie.activate!()
psi = load_from_text(joinpath(TEMP_DATA_DIR, "initial_state.gpe.dat"))
#psi = load_from_text(joinpath(DATA_DIR, "initial_state.gpe.dat"))
plot_iso_bce(psi)

##

enes = load_energies(joinpath(DATA_DIR, "energy.dat"))
plot_energies(enes)

##
using GLMakie
GLMakie.activate!()
##

using CairoMakie
CairoMakie.activate!()
##
#animate_iso_bce(TEMP_DATA_DIR, "bce_evolution.gif")

##

#psi_vec = load_directory_from_text(TEMP_DATA_DIR)

##
fig = plot_local_maxima_evolution(psi_vec)
save("eps_1_5_atoms_30k_1500nm.pdf", fig)

##

fig = plot_local_maxima_coordinates(psi_vec)
save("coordinates_eps_1_5_atoms_30k_1500nm.pdf", fig)

##

plot_slice(interpolate_slice(get_BEC_slice(psi_vec[begin])))

##

psi = load_from_text(joinpath("../../../data/run_find_initial_states", "6k_atoms/2_max/initial_state.dat"))
#maxs = find_local_maxima(get_BEC_slice(psi))
#plot_iso_bce(psi)
plot_slice(interpolate_slice(get_BEC_slice(psi)))
#maxs = number_of_lmax(psi)
#
##
psi_1 = get_BEC_slice(load_from_text(joinpath("../../../data/run_find_initial_states", "7k_atoms/1_max/initial_state.gpe.dat")))
psi_2 = get_BEC_slice(load_from_text(joinpath("../../../data/run_find_initial_states", "7k_atoms/2_max/initial_state.gpe.dat")))

##
plot_slice(psi_1)
##
dist_dir = "../../../data/run_find_initial_states"
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
