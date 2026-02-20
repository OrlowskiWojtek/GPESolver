using GLMakie

include("files.jl")
include("plots.jl")
include("animate.jl")

##

fig = plot_evolution("../../../data/run_30_atoms"; step = 10, total_size = 6)
save("time_3d_edd_15_dd_1500_atoms_30000.png", fig)

##

fig = plot_evolution("../../../data/run_40_atoms"; step = 10, total_size = 6);
save("time_3d_edd_15_dd_1500_atoms_40000.png", fig)

##

fig = plot_evolution("../../../data/run_40_atoms_eps_145"; step = 10, total_size = 6);
save("time_3d_edd_145_dd_1500_atoms_40000.png", fig)

##

fig = plot_evolution("../../../data/run_275_atoms_eps_145"; step = 10, total_size = 6);
save("time_3d_edd_145_dd_1500_atoms_27500.png", fig)

##

minimas_edd_15 = ["1k_atoms/1_max/initial_state.gpe.dat",
                   "6k_atoms/1_max/initial_state.gpe.dat",
                   "8k_atoms/2_max/initial_state.gpe.dat",
                   "27k_atoms/2_max/initial_state.dat",
                   "30k_atoms/3_max/initial_state.dat",
                   "37k_atoms/4_max/initial_state.dat",
                    ];

paths_edd_15 = joinpath.("../../../data/run_find_initial_states/", minimas_edd_15)
fig = plot_states(paths_edd_15);
save("edd_15_ground_states.png", fig);

##

minimas_edd_145 = ["1k_atoms/2_max/initial_state.gpe.dat",
                   "2k_atoms/2_max/initial_state.gpe.dat",
                   "5k_atoms/1_max/initial_state.gpe.dat",
                   "27k_atoms/2_max/initial_state.gpe.dat",
                   "30k_atoms/3_max/initial_state.gpe.dat",
                   "49k_atoms/5_max/initial_state.gpe.dat",
                    ];

paths_edd_145 = joinpath.("../../../data/run_find_initial_states_eps145/", minimas_edd_145)
fig = plot_states(paths_edd_145);
save("edd_145_ground_states.png", fig);
