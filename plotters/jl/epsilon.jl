using CairoMakie
using DataFrames

function gather_energy(dist_dir::String)
    df = DataFrame(epsilon=Float64[],
                   e_kin=Float64[],
                   e_pot=Float64[],
                   e_int=Float64[],
                   e_ext=Float64[],
                   e_bmf=Float64[],
                   e_total=Float64[],
                   filename=String[])

    for eps_folder in 140:160
        eps_path = joinpath(dist_dir, "$(eps_folder)eps")
        if !isdir(eps_path)
            continue
        end

        energy_file = joinpath(eps_path, "energy.gpe.dat")

        if !isfile(energy_file)
            continue
        end

        open(energy_file, "r") do f
            lines = readlines(f)
            line = lines[end]
            line = strip(line)
            if isempty(line)
                return
            end
            numbers = split(line, '\t')
            if length(numbers) < 7
                return
            end

            push!(df, (eps_folder,
                parse(Float64, numbers[2]),
                parse(Float64, numbers[3]),
                parse(Float64, numbers[4]),
                parse(Float64, numbers[5]),
                parse(Float64, numbers[6]),
                parse(Float64, numbers[7]),
                eps_path))
        end
    end

    return df
end

function plot_energies(df)
    fig = Figure();
    ax  = Axis(fig[1,1],
               title = "N = 4e4",
               xlabel = "εdd",
               ylabel = "E/N [nK]");

    df.e_tot_ext = df.e_int .+ df.e_ext .+ df.e_bmf

    labels = ["Total", "Kinetic", "Potential", "Interaction"]
    fields = [:e_total, :e_kin, :e_pot, :e_tot_ext]
    colors = [:black, :orange, :red, :blue]
    conv = 27211.4 * 11.6 * 10^9 / 40000

    for idx in eachindex(labels)
        lines!(ax, df[!, :epsilon] ./ 100,
               df[!, fields[idx]] .* conv,
               color = colors[idx],
               label = "$(fields[idx])"
               )
    end

    Legend(fig[2,1], ax, orientation = :horizontal)
    fig
end

##

df = gather_energy("../../../data/run_find_epsilon")

##

fig = plot_energies(df)
save("ene_vs_epsilon_n4e4.pdf", fig)
