using GLMakie

include("files.jl")
include("plots.jl")
include("animate.jl")

##

#psi = load_from_binary(joinpath(DATA_DIR, "checkpoint_1.bin"))
GLMakie.activate!()
psi = load_from_text(joinpath(TEMP_DATA_DIR, "30k_atoms/3_max/initial_state.dat"))
plot_iso_bce(psi)

##

enes = load_energies(joinpath(DATA_DIR, "energy.dat"))
plot_energies(enes)

##
using GLMakie
GLMakie.activate!()
##
animate_iso_bce("bce_evolution.gif")

##


