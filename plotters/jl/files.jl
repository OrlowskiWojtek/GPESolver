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
    return IsoBCEContext(array_3d, x, y, z, nx, ny, nz, dx, dy, dz)
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
