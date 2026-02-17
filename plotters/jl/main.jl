using GLMakie

include("files.jl")
include("plots.jl")
include("animate.jl")

##

GLMakie.activate!()
#slice = load_slice_from_text(joinpath("../../../data/run_find_initial_states_eps145/50k_atoms/4_max", "initial_state.gpe.dat"))
psi = load_from_text(joinpath(TEMP_DATA_DIR, "wavefunction_250.gpe.dat"))
plot_iso_bce(psi)

##

animate_iso_bce("../../../data/run_30_atoms", "animation_30k_eps_150.gif")

##

enes = load_energies(joinpath(TEMP_DATA_DIR, "energy.gpe.dat"))
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
save("eps_145_atoms_27k5_1500nm.pdf", fig)

##

fig = plot_local_maxima_coordinates(slice_vec)
save("coordinates_eps_145_atoms_27k5_1500nm.pdf", fig)

##

psi = load_from_text(joinpath("../../../data/run_find_initial_states_eps145/40k_atoms/4_max", "initial_state.gpe.dat"))

##

get_BEC_heights(psi)

## Making plot composition using GLMakie
# 3D plots are required as they review correct approach to time evulotion

#fig = plot_evolution(TEMP_DATA_DIR; step = 40)
save("time_test.png", fig)
