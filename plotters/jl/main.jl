using GLMakie

include("files.jl")
include("plots.jl")
include("animate.jl")
##

animate_iso_bce("../../../data/run_30_atoms", "animation_30k_eps_150.gif")

##

enes = load_energies(joinpath("../../../data/run_275_atoms/", "energy.gpe.dat"))
fig = plot_energies(enes)

##

using GLMakie
GLMakie.activate!()

##

using CairoMakie
CairoMakie.activate!()

##

animate_iso_bce("../../../data/run_275_atoms/", "bce_evolution.gif")

##

slice_vec = load_directory_slice_from_text("../../../data/run_40_atoms")

##

fig = plot_local_maxima_evolution(slice_vec);
#save("eps_145_atoms_40k_1500nm.pdf", fig);

##

fig = plot_local_maxima_coordinates_one_ax(slice_vec)
#save("coordinates_eps_150_atoms_40k_1500nm.pdf", fig)

##

psi = load_from_text(joinpath("../../../data/run_find_initial_states/20k_atoms/4_max", "initial_state.dat"))

##

get_BEC_heights(psi)

##

get_BEC_maxrho(psi) .* 20000


## Making plot composition using GLMakie
# 3D plots are required as they review correct approach to time evulotion

fig = plot_evolution("../../../"; step = 10, total_size = 6);
save("time_3d_edd_15_dd_1500_atoms_30000.png", fig)

##

fig = plot_evolution(""; step = 10, total_size = 6);
save("time_3d_edd_15_dd_1500_atoms_30000.png", fig)


