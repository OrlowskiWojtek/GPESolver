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
                   "40k_atoms/4_max/initial_state.dat",
                    ];

#paths_edd_15 = joinpath.("../../../data/run_find_initial_states/", minimas_edd_15)
#fig = plot_states(paths_edd_15);
#save("edd_15_ground_states.png", fig);

##

for (idx, minima_edd_15) in enumerate(minimas_edd_15)
    file = joinpath("../../../data/run_find_initial_states", minima_edd_15)
    fig = plot_single_state(file)
    N = match(r"(\d+)k_atoms", file)
    maxs = match(r"(\d+)_max", file)

    dir = "plots/eps15/"
    filename = joinpath(dir, "wavefunction_" * N[1] * "k_atoms_" * maxs[1] * "_condensates" * ".png")

    save(filename, fig)
end

##

minimas_edd_145 = ["10k_atoms/2_max/initial_state.gpe.dat",
                   "20k_atoms/2_max/initial_state.gpe.dat",
                   "27k_atoms/3_max/initial_state.gpe.dat",
                   "40k_atoms/4_max/initial_state.gpe.dat",
                   "50k_atoms/5_max/initial_state.gpe.dat",
                    ];

#paths_edd_145 = joinpath.("../../../data/run_find_initial_states_eps145/", minimas_edd_145)
#fig = plot_states(paths_edd_145);
#save("edd_145_ground_states.png", fig);

##

for (idx, minima_edd_145) in enumerate(minimas_edd_145)
    file = joinpath("../../../data/run_find_initial_states_eps145", minima_edd_145)
    fig = plot_single_state(file)
    N = match(r"(\d+)k_atoms", file)
    maxs = match(r"(\d+)_max", file)

    dir = "plots/eps145/"
    filename = joinpath(dir, "wavefunction_" * N[1] * "k_atoms_" * maxs[1] * "_condensates" * ".png")

    save(filename, fig)
end

##

file = joinpath("../../../data/run_find_initial_states", minimas_edd_15[end])
fig = plot_single_state(file)

##
#
    %\includegraphics[height=0.18\textwidth, trim={1.1cm 3cm 0 0},clip]{1500/1k45/1k45_1e4_1.png} & 2
    %\includegraphics[height=0.18\textwidth, trim={1.1cm 3cm 0 0},clip]{1500/1k45/x15001k45_2e4_1.png} & 2
    %\includegraphics[height=0.18\textwidth, trim={1.1cm 3cm 0 0},clip]{1500/1k45/f1500_1k45_2k75_1.png} & 3
    %\includegraphics[height=0.18\textwidth, trim={1.1cm 3cm 0 0},clip]{1500/1k45/fr1k5_1k45_4e4_1.png} & 4
    %\includegraphics[height=0.18\textwidth, trim={1.1cm 3cm 0 0},clip]{1500/1k45/f1500_1k45_5e4__1.png} 5
