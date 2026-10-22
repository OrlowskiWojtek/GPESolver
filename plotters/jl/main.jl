include("files.jl")
include("plots.jl")
include("animate.jl")

##

animate_slice_bce("../../build/long_from_1_35_to_1_45", "long_with_random_noise.gif")

##

enes = load_energies(joinpath("../../build/ground_states/eps_1_45/gauss", "energy.gpe.dat"))
fig = plot_energies(enes)

##

##

animate_iso_bce("../../build/long_from_1_35_to_1_45", "long_with_random_noise_3d.gif")

##

slice_vec = load_directory_slice_from_binary("../../build/run_full_time_n30000_edd15")

##

fig = plot_local_maxima_evolution(slice_vec);
#save("eps_145_atoms_40k_1500nm.pdf", fig);

##

#fig = plot_local_maxima_coordinates_one_ax(slice_vec)
##save("coordinates_eps_150_atoms_40k_1500nm.pdf", fig)

save_local_maxima_to_file(slice_vec, "n30000_time.dat", 30000)


##

psi = load_from_text(joinpath("../../build/", "initial_state.gpe.dat"))
#psi = load_from_binary(joinpath("../../build/test_binary", "wavefunction_1.gpe.bin"))

#slice = get_BEC_slice(psi)
#f = plot_slice(slice, 40000.)
plot_iso_bce(psi)

#save("tmp.png", f)
##

x0_idx = 129
y0_idx = 129
rho = abs.(psi.psi[x0_idx, y0_idx, :])

FHWM = maximum(rho) / 2.
##
##

##
##

animate_iso_bce_interactive_prefetch("../../build/test_binary")
#animate_slice_bce("../../build/speed_test", "tmp.gif")

##

pote = load_pote_from_text("../../build/initial_potential.gpe.dat")
pote2 = load_pote_from_text("../../build/noised_pote.gpe.dat")
conv = 27211.4 * 11.6 * 10^9
fig = plot_potential_comparison(pote, pote2; conversion = conv)
