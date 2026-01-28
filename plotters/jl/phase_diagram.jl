using DataFrames
using CSV
using Glob

include("context.jl")

# another approach; - load all output files, then find number of maximas assiociated to energy
function gather_energy(dist_dir::String; BCE_THRESHOLD = nothing)
    df = DataFrame(atom_number=Int[],
                   max_number=Int[],
                   e_kin=Float64[],
                   e_pot=Float64[],
                   e_int=Float64[],
                   e_ext=Float64[],
                   e_bmf=Float64[],
                   e_total=Float64[])

    for atom_folder in 1:50
        atom_dir = joinpath(dist_dir, "$(atom_folder)k_atoms")
        if !isdir(atom_dir)
            continue
        end

        max_folders = glob("*_max", atom_dir)
        for max_folder in max_folders
            energy_file = joinpath(max_folder, "energy.gpe.dat")

            if !isfile(energy_file)
                energy_file = joinpath(max_folder, "energy.dat")
            end

            if !isfile(energy_file)
                continue
            end

            psi_file = joinpath(max_folder, "initial_state.gpe.dat")

            if !isfile(psi_file)
                psi_file = joinpath(max_folder, "initial_state.dat")
            end

            if !isfile(psi_file)
                @error "No wavefunction file"
            end

            open(energy_file, "r") do f
                line = readline(f)
                line = strip(line)
                if isempty(line)
                    return
                end
                numbers = split(line, '\t')
                if length(numbers) < 7
                    return
                end
                atom_number = parse(Int, match(r"(\d+)k_atoms", atom_dir).captures[1])
                max_number = parse(Int, match(r"(\d+)_max", basename(max_folder)).captures[1])
                psi = load_from_text(psi_file)
                l_maxes = number_of_lmax(psi;n_atoms = atom_number, condensation_threshold = 1e-6)

                if(l_maxes != max_number)
                    @info "Number of condensates changed from $max_number to $l_maxes in $atom_number k atoms"
                    max_number = l_maxes
                end

                push!(df, (atom_number,
                    max_number,
                    parse(Float64, numbers[2]),
                    parse(Float64, numbers[3]),
                    parse(Float64, numbers[4]),
                    parse(Float64, numbers[5]),
                    parse(Float64, numbers[6]),
                    parse(Float64, numbers[7])))
            end
        end
    end

    #CSV.write("energies.txt", df, delim=' ')
    return df
end

##

gather_energy("../../../data/run_find_initial_states")
