using GLMakie

include("files.jl")
include("plots.jl")
include("animate.jl")

##

GLMakie.activate!()
#slice = load_slice_from_text(joinpath("../../../data/run_find_initial_states_eps145/50k_atoms/4_max", "initial_state.gpe.dat"))
psi = load_from_text(joinpath("../../../data/run_find_initial_states_eps145/40k_atoms/4_max", "initial_state.gpe.dat"))
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

slice_vec = load_directory_slice_from_text(TEMP_DATA_DIR)

##

fig = plot_local_maxima_evolution(slice_vec)
save("eps_1_5_atoms_30k_1500nm.pdf", fig)

##

fig = plot_local_maxima_coordinates(slice_vec)
save("coordinates_eps_1_5_atoms_30k_1500nm.pdf", fig)

##

psi = load_from_text(joinpath("../../../data/run_find_initial_states_eps145/40k_atoms/4_max", "initial_state.gpe.dat"))
##
get_BEC_heights(psi)
