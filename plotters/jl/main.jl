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

psi_vec = load_directory_from_text(TEMP_DATA_DIR)

##
fig = plot_local_maxima_evolution(psi_vec)
save("eps_1_5_atoms_30k_1500nm.pdf", fig)

##

fig = plot_local_maxima_coordinates(psi_vec)
save("coordinates_eps_1_5_atoms_30k_1500nm.pdf", fig)

##

plot_slice(interpolate_slice(get_BEC_slice(psi_vec[begin])))

##

psi = load_from_text(joinpath("../../../../data/run_find_initial_states",
                              "7k_atoms/2_max/initial_state.dat"))
#maxs = find_local_maxima(get_BEC_slice(psi))
plot_iso_bce(psi)
#plot_slice(interpolate_slice(get_BEC_slice(psi)))
#maxs = number_of_lmax(psi)

##
psi_1 = get_BEC_slice(load_from_text(joinpath("../../../../data/run_find_initial_states", "10k_atoms/1_max/initial_state.gpe.dat")))
psi_2 = get_BEC_slice(load_from_text(joinpath("../../../../data/run_find_initial_states", "10k_atoms/2_max/initial_state.gpe.dat")))

plot_slice(psi_2)
