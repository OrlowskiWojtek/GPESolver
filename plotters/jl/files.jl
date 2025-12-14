include("context.jl")

const DATA_DIR = "../../build/"

function load_from_binary(file_path::String)
    file    = open(file_path, "r")
    nx      = read(file, Int32)
    ny      = read(file, Int32)
    nz      = read(file, Int32)

    dx      = 200
    dy      = 200
    dz      = 500

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

    x = collect(-div(nx, 2):div(nx, 2)) * dx
    y = collect(-div(ny, 2):div(ny, 2)) * dy
    z = collect(-div(nz, 2):div(nz, 2)) * dz

    close(file)
    return IsoBECContext(array_3d, x, y, z, nx, ny, nz, dx, dy, dz)
end

function load_from_text(file_path::String)
    file    = open(file_path, "r")
    nx      = parse(Int32, readline(file))
    ny      = parse(Int32, readline(file))
    nz      = parse(Int32, readline(file))

    dx      = 200
    dy      = 200
    dz      = 500

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

    x = collect(-div(nx, 2):div(nx, 2)) * dx
    y = collect(-div(ny, 2):div(ny, 2)) * dy
    z = collect(-div(nz, 2):div(nz, 2)) * dz

    close(file)
    return IsoBECContext(array_3d, x, y, z, nx, ny, nz, dx, dy, dz)
end

function load_xy_cut(file_path::String)
    file = open(file_path, "r")

    nx::Int32 = 0
    ny::Int32 = 0
    while(true)
        line = readline(file)

        if contains(line, "# X	Y	|Psi|^2")
            break
        end

        if(contains(line, "nx"))
            splitted = split(line, ':')
            nx = parse(Int32, splitted[2])
        end
        if(contains(line, "ny"))
            splitted = split(line, ':')
            ny = parse(Int32, splitted[2])
        end
        if(contains(line, "nz"))
            splitted = split(line, ':')
            nz = parse(Int32, splitted[2])
        end
    end

    rho = zeros(Float64, nx, ny)
    println("nx: $nx, ny: $ny")

    for i in 1:nx
        for j in 1:ny
            line = readline(file)
            splitted = split(line)
            rho[i, j] = parse(Float64, splitted[3])
        end
    end

    close(file)
    return rho
end

function load_from_fort(file_path::String)
    file    = open(file_path, "r")

    nline = readline(file)
    splitted = split(nline)
    nx = 2 * parse(Int32, splitted[1]) + 1
    ny = 2 * parse(Int32, splitted[2]) + 1
    nz = 2 * parse(Int32, splitted[3]) + 1

    # default - to change

    x = zeros(Float64, nx)
    y = zeros(Float64, ny)
    z = zeros(Float64, nz)
    array_3d = zeros(ComplexF64, nx, ny, nz)
    for i in 1:nx
        for j in 1:ny
            for k in 1:nz
                line = readline(file)
                splitted = split(line)
                x[i] = parse(Float64, splitted[1])
                y[j] = parse(Float64, splitted[2])
                z[k] = parse(Float64, splitted[3])
                real_part = parse(Float64, splitted[4])
                imag_part = parse(Float64, splitted[5])
                array_3d[i, j, k] = ComplexF64(real_part, imag_part)
            end
        end
    end

    dx = x[2] - x[1]
    dy = y[2] - y[1]
    dz = z[2] - z[1]
    
    return IsoBECContext(array_3d, x, y, z, nx, ny, nz, dx, dy, dz)
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
