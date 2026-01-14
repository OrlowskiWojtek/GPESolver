using GLMakie

include("files.jl")
include("plots.jl")
include("animate.jl")

##

#psi = load_from_binary(joinpath(DATA_DIR, "checkpoint_1.bin"))
psi = load_from_text(joinpath(DATA_DIR, "initial_state.gpe.dat"))
plot_iso_bce(psi)

##

enes = load_energies(joinpath(DATA_DIR, "energy.dat"))
plot_energies(enes)

##
using GLMakie
GLMakie.activate!()
##
animate_iso_bce("bce_evolution.gif")
