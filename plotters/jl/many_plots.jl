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
                   "38k_atoms/3_max/initial_state.dat",
                   "40k_atoms/4_max/initial_state.dat",
                    ];

paths_edd_15 = joinpath.("../../../data/run_find_initial_states/", minimas_edd_15)
fig = plot_single_state(paths_edd_15[2], hide_decs = false, n_atoms = 6000);
save("edd_15_test.png", fig, px_per_unit = 2);

##

for (idx, minima_edd_15) in enumerate(minimas_edd_15)
    file = joinpath("../../../data/run_find_initial_states", minima_edd_15)
    N = match(r"(\d+)k_atoms", file)
    maxs = match(r"(\d+)_max", file)

    fig = plot_single_state(file; hide_decs = (idx != 1), n_atoms = parse(Float64, N[1]) * 1000)
    dir = "plots/eps15/"
    filename = joinpath(dir, "wavefunction_" * N[1] * "k_atoms_" * maxs[1] * "_condensates" * ".png")

    save(filename, fig)
end

##

file = joinpath("../../../data/run_find_initial_states", minimas_edd_15[begin])
fig = plot_single_state(file, hide_decs = false)

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
    N = match(r"(\d+)k_atoms", file)
    maxs = match(r"(\d+)_max", file)
    fig = plot_single_state(file, hide_decs = (idx != 1), n_atoms = parse(Float64, N[1]) * 1000)

    dir = "plots/eps145/"
    filename = joinpath(dir, "wavefunction_" * N[1] * "k_atoms_" * maxs[1] * "_condensates" * ".png")

    save(filename, fig)
end

##

minimas_single_well = ["5k_atoms/1_max/initial_state.gpe.dat",
    "6k_atoms/1_max/initial_state.gpe.dat",
    "7k5_atoms/1_max/initial_state.gpe.dat",
    "20k_atoms/2_max/initial_state.gpe.dat",
    "27k5_atoms/3_max/initial_state.gpe.dat",
    "40k_atoms/4_max/initial_state.gpe.dat",
];

minimas_single_well_atom_counts = [
    5000,   # 5k_atoms
    6000,   # 6k_atoms
    7500,   # 7k5_atoms
    20000,  # 20k_atoms
    27500,  # 27k5_atoms
    40000,  # 40k_atoms
];

for (idx, minima_single_well) in enumerate(minimas_single_well)
    file = joinpath("../../../data/run_find_initial_states_single_well", minima_single_well)
    N = minimas_single_well_atom_counts[idx]
    maxs = match(r"(\d+)_max", file)
    fig = plot_single_state(file; hide_decs = (idx != 1), n_atoms = N[1])

    dir = "plots/single_well/"
    filename = joinpath(dir, "wavefunction_$(N)_atoms_" * maxs[1] * "_condensates" * ".png")

    save(filename, fig)
end

##

minimas_eps_140 = [
    "20k_atoms/2_max/initial_state.gpe.dat",
    "40k_atoms/4_max/initial_state.gpe.dat"
];

minimas_eps_140_atom_counts = [
    20000,   # 5k_atoms
    40000
];

for (idx, minima_eps_140) in enumerate(minimas_eps_140)
    file = joinpath("../../../data/run_find_initial_states_eps140", minima_eps_140)
    N = minimas_eps_140_atom_counts[idx]
    maxs = match(r"(\d+)_max", file)

    fig = plot_single_state(file; hide_decs = (idx != 1), n_atoms = N)

    dir = "plots/eps140/"
    filename = joinpath(dir, "wavefunction_$(N)_atoms_" * maxs[1] * "_condensates" * ".png")
    display(fig)

    #save(filename, fig)
end
