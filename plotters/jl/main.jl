using GLMakie

include("files.jl")
include("plots.jl")
include("animate.jl")
##

test = load_xy_cut(joinpath(DATA_DIR, "cut_xy_1000.dat"))
fig = plot_heatmap_cut(test)

##

#psi = load_from_binary(joinpath(DATA_DIR, "checkpoint_1.bin"))
#psi = load_from_binary(joinpath(DATA_DIR, "last_state.bin"))
psi = load_from_text(joinpath("../../data/run_find_initial_states/50k_atoms/2_max", "initial_state.dat"))
#psi = load_from_fort(joinpath(DATA_DIR, "ff.dat"))
plot_iso_bce(psi)

##

enes = load_energies(joinpath(DATA_DIR, "energy.dat"))
plot_energies(enes)

##
using GLMakie
GLMakie.activate!()
##
animate_iso_bce("bce_evolution.gif")
