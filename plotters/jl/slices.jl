## Here I am making 1D slices along x and z directions in (0,0,z) and (x, 0, 0)

using CairoMakie

include("files.jl")
include("plots.jl")
include("animate.jl")

CairoMakie.activate!()
using DelimitedFiles

##

filenames = [
    "../../build/harmonic_1x1.gpe.dat",
    "../../build/harmonic_2x1.gpe.dat",
    "../../build/regular_1x1.gpe.dat",
    "../../build/regular_2x1.gpe.dat",
    "../../build/mexican_2x1.gpe.dat"
]

##


function load_and_save(filename)
    data = load_from_text(filename)

    x_0_index = round(Int32, data.nx / 2 + 1)
    y_0_index = round(Int32, data.ny / 2 + 1)
    z_0_index = round(Int32, data.nz / 2 + 1)

    psi_x = data.psi[:, y_0_index, z_0_index]
    psi_z = data.psi[x_0_index, y_0_index, :]
    rho_x = abs.(psi_x) .^ 2 ./ length_au3_to_μm3(1.) * 10000
    rho_z = abs.(psi_z) .^ 2 ./ length_au3_to_μm3(1.) * 10000

    fig = Figure();
    ax = Axis(fig[1,1],
              xlabel = "x [nm]",
              ylabel = "ρ [μm⁻³]");

    lines!(ax, data.z, rho_z, color = :red,  label = "z")
    lines!(ax, data.x, rho_x, color = :blue, label = "x")

    axislegend()

    output = filename[begin+12:end-7]*"pdf"

    save(output, fig)

    output_x = filename[begin+12:end-8]*"_x.txt"
    output_z = filename[begin+12:end-8]*"_z.txt"

    open(output_x, "w") do io
        writedlm(io, [data.x rho_x], '\t')
    end;
    open(output_z, "w") do io
        writedlm(io, [data.z rho_z], '\t')
    end;
end


for file in filenames
    load_and_save(file)
end
