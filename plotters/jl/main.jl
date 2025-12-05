using GLMakie

include("files.jl")
include("plots.jl")
include("animate.jl")
##
test = load_xy_cut(joinpath(DATA_DIR, "cut_xy_1000.dat"))
fig = plot_heatmap_cut(test)
##
using GLMakie
GLMakie.activate!()

#psi = load_from_binary(joinpath(DATA_DIR, "checkpoint_50000.bin"))
psi = load_from_text(joinpath(DATA_DIR, "data_2_cons/initial_state.dat"))
plot_iso_bce(psi)

##
animate_iso_bce("bce_evolution.gif")
