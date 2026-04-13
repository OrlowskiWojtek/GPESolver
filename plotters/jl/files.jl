include("context.jl")

const DATA_DIR = "../../build/"
TEMP_DATA_DIR = "../../../data/run_30_atoms"

function load_from_binary(file_path::String)
    file    = open(file_path, "r")
    nx      = read(file, Int32)
    ny      = read(file, Int32)
    nz      = read(file, Int32)

    dx      = read(file, Float64)
    dy      = read(file, Float64)
    dz      = read(file, Float64)

    array_3d = zeros(ComplexF64, nx, ny, nz)

    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                real_part = read(file, Float64)
                imag_part = read(file, Float64)
                array_3d[i, j, k] = ComplexF64(real_part, imag_part)
            end
        end
    end

    x = ((1:nx) .- (div(nx, 2) + 1)) .* dx
    y = ((1:ny) .- (div(ny, 2) + 1)) .* dy
    z = ((1:nz) .- (div(nz, 2) + 1)) .* dz

    close(file)
    return IsoBECContext(array_3d, x, y, z, nx, ny, nz, dx, dy, dz)
end

function load_from_text(file_path::String)
    file    = open(file_path, "r")
    nx      = parse(Int32, readline(file))
    ny      = parse(Int32, readline(file))
    nz      = parse(Int32, readline(file))

    dx      = parse(Float64, readline(file))
    dy      = parse(Float64, readline(file))
    dz      = parse(Float64, readline(file))

    array_3d = zeros(ComplexF64, nx, ny, nz)

    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                line = readline(file)
                splitted = split(line, "\t")
                real_part = parse(Float64, splitted[1])
                imag_part = parse(Float64, splitted[2])
                array_3d[i, j, k] = ComplexF64(real_part, imag_part)
            end
        end
    end

    x = ((1:nx) .- (div(nx, 2) + 1)) .* dx
    y = ((1:ny) .- (div(ny, 2) + 1)) .* dy
    z = ((1:nz) .- (div(nz, 2) + 1)) .* dz

    close(file)
    return IsoBECContext(array_3d, x, y, z, nx, ny, nz, dx, dy, dz)
end

function load_pote_from_text(file_path::String)
    file    = open(file_path, "r")

    nx      = parse(Int32, split(readline(file), '=')[2])
    ny      = parse(Int32, split(readline(file), '=')[2])
    nz      = parse(Int32, split(readline(file), '=')[2])

    dx      = parse(Float64, split(readline(file), '=')[2])
    dy      = parse(Float64, split(readline(file), '=')[2])
    dz      = parse(Float64, split(readline(file), '=')[2])

    array_3d = zeros(Float64, nx, ny, nz)

    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                line = readline(file)
                array_3d[i, j, k] =  parse(Float64, line)
            end
        end
    end

    x = ((1:nx) .- (div(nx, 2) + 1)) .* dx
    y = ((1:ny) .- (div(ny, 2) + 1)) .* dy
    z = ((1:nz) .- (div(nz, 2) + 1)) .* dz

    close(file)
    return PoteContext(array_3d, x, y, z, nx, ny, nz, dx, dy, dz)
end

function save_to_text(psi::IsoBECContext, file_path::String)
    file = open(file_path, "w")

    # Write header information
    println(file, psi.nx)
    println(file, psi.ny)
    println(file, psi.nz)
    println(file, psi.dx)
    println(file, psi.dy)
    println(file, psi.dz)

    # Write 3D array data
    for i in 1:psi.nx
        for j in 1:psi.ny
            for k in 1:psi.nz
                val = psi.psi[i, j, k]
                println(file, "$(real(val))\t$(imag(val))")
            end
        end
    end

    close(file)
end

function load_energies(filename::String)
    iter = Int32[]
    e_kin = Float64[]
    e_pot = Float64[]
    e_int = Float64[]
    e_ext = Float64[]
    e_bmf = Float64[]
    e_tot = Float64[]

    open(filename, "r") do file
        for line in eachline(file)
            parts = split(line)
            push!(iter, parse(Int32, parts[1]))
            push!(e_kin, parse(Float64, parts[2]))
            push!(e_pot, parse(Float64, parts[3]))
            push!(e_int, parse(Float64, parts[4]))
            push!(e_ext, parse(Float64, parts[5]))
            push!(e_bmf, parse(Float64, parts[6]))
            push!(e_tot, parse(Float64, parts[7]))
        end
    end

    return EnergiesContext(iter, e_kin, e_pot, e_int, e_ext, e_bmf, e_tot)
end

# Function to load data from quick directory.
function load_directory_from_text(data_dir::String)
    files = filter(f -> occursin("wavefunction_", f) && endswith(f, ".gpe.dat"), readdir(data_dir))
    files = sort(files, by = f -> parse(Int, split(split(f, "_")[2], ".")[1]))
    n_frames = length(files)

    bec_data_vec = Vector{IsoBECContext}(undef, n_frames)

    for (frame_idx, file) in enumerate(files)
        println("loading frame", joinpath(data_dir, file))
        context = load_from_text(joinpath(data_dir, file))
    
        bec_data_vec[frame_idx] = context 
    end

    return bec_data_vec
end

function load_slice_from_text(file_path::String)
    file = open(file_path, "r")
    nx = parse(Int32, readline(file))
    ny = parse(Int32, readline(file))
    nz = parse(Int32, readline(file))

    dx = parse(Float64, readline(file))
    dy = parse(Float64, readline(file))
    dz = parse(Float64, readline(file))

    # z=0 is at k0:
    k0 = div(nz, 2) + 1

    slice_2d = zeros(ComplexF64, nx, ny)

    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                line = readline(file)
                if k == k0
                    splitted = split(line, "\t")
                    real_part = parse(Float64, splitted[1])
                    imag_part = parse(Float64, splitted[2])
                    slice_2d[i, j] = ComplexF64(real_part, imag_part)
                end
            end
        end
    end

    x = ((1:nx) .- (div(nx, 2) + 1)) .* dx
    y = ((1:ny) .- (div(ny, 2) + 1)) .* dy

    close(file)
    return IsoBECSlice(slice_2d, x, y, nx, ny, dx, dy)
end

function load_directory_slice_from_text(data_dir::String)
    files = filter(f -> occursin("wavefunction_", f) && endswith(f, ".gpe.dat"), readdir(data_dir))
    files = sort(files, by = f -> parse(Int, split(split(f, "_")[2], ".")[1]))
    n_frames = length(files)

    slice_data_vec = Vector{IsoBECSlice}(undef, n_frames)

    for (frame_idx, file) in enumerate(files)
        println("loading frame", joinpath(data_dir, file))
        context = load_slice_from_text(joinpath(data_dir, file))
    
        slice_data_vec[frame_idx] = context 
    end

    return slice_data_vec
end
