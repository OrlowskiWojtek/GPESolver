include("units.jl")
using Interpolations

X_INTERPOLATION_SIZE = 2000
Y_INTERPOLATION_SIZE = 2000
Z_INTERPOLATION_SIZE = 500

struct IsoBECContext
    psi::Array{ComplexF64,3}  # 3D array for the wavefunction data
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    nx::Int32
    ny::Int32
    nz::Int32
    dx::Float64
    dy::Float64
    dz::Float64
end

struct PoteContext
    pote::Array{Float64,3}  # 3D array for the potential data
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    nx::Int32
    ny::Int32
    nz::Int32
    dx::Float64
    dy::Float64
    dz::Float64
end

struct EnergiesContext
    iter::Vector{Int32}
    e_kin::Vector{Float64}
    e_pot::Vector{Float64}    
    e_int::Vector{Float64}    
    e_ext::Vector{Float64}    
    e_bmf::Vector{Float64}    
    e_tot::Vector{Float64}    
end

struct LocalMaximaGrid
    i::Int32
    j::Int32
    value::Float64
end

struct LocalMaximaPhysical
    x::Float64
    y::Float64
    value::Float64
end

struct BECHeight
    i::Int32
    j::Int32
    x::Float64
    y::Float64
    height::Float64
    value::Float64
end

struct BECMaxRho
    value::Float64
    n_atoms::Float64
end

struct IsoBECSlice
    psi::Array{ComplexF64,2}
    x::Vector{Float64}
    y::Vector{Float64}
    nx::Int32
    ny::Int32
    dx::Float64
    dy::Float64
end

function get_BEC_slice(context::IsoBECContext; z_0_index_offset = 1)
    #@assert context.nz % 2 == 1 "Wrong grid in 'z' dimension"
    z_0_index = Int64(floor(context.nz / 2) + z_0_index_offset)
    
    psi_slice = context.psi[:, :, z_0_index]  # 2D slice at z=0
    return IsoBECSlice(
        psi_slice,
        context.x,
        context.y,
        context.nx,
        context.ny,
        context.dx,
        context.dy
    )
end

function find_local_maxima(slice::IsoBECSlice)
    rho = abs.(slice.psi)
    max_val = maximum(rho)
    threshold_val = 0.6 * max_val
    
    local_maxima = Vector{LocalMaximaGrid}()
    
    for j in 2:(slice.ny-1)
        for i in 2:(slice.nx-1)
            current_val = rho[i, j]
            
            if current_val < threshold_val
                continue
            end
            
            is_maximum = true
            for dj in -1:1
                for di in -1:1
                    if di == 0 && dj == 0
                        continue
                    end
                    if rho[i+di, j+dj] > current_val
                        is_maximum = false
                        break
                    end
                end
                if !is_maximum
                    break
                end
            end
            
            if is_maximum
                push!(local_maxima, LocalMaximaGrid(i, j, current_val))
            end
        end
    end
    
    return local_maxima
end

function get_coordinates(slice::IsoBECSlice, maxima::Vector{LocalMaximaGrid})
    coords = [LocalMaximaPhysical(slice.x[m.i], slice.y[m.j], m.value) for m in maxima]
    return coords
end

function number_of_lmax(context::IsoBECContext; n_atoms::Int64 = 1, condensation_threshold = nothing)
    threshold = 0.6
    n_max = 0;

    slice   = get_BEC_slice(context)
    rho     = abs.(slice.psi)

    max_val = threshold * maximum(rho)
    total_maximas = find_local_maxima(slice)

    for m in total_maximas
        if(condensation_threshold != nothing)
            if(m.value * n_atoms < condensation_threshold)
                @info "Low number of atoms in maximum: $(m.value)"
                continue 
            end
        end
        if(m.value > max_val)
            n_max += 1
        end
    end

    return n_max
end

function number_of_lmax(slice::IsoBECSlice; n_atoms::Int64 = 1, condensation_threshold = nothing)
    threshold = 0.6
    n_max = 0;

    rho     = abs.(slice.psi)

    max_val = threshold * maximum(rho)
    total_maximas = find_local_maxima(slice)

    for m in total_maximas
        if(condensation_threshold != nothing)
            if(m.value * n_atoms < condensation_threshold)
                @info "Low number of atoms in maximum: $(m.value)"
                continue 
            end
        end
        if(m.value > max_val)
            n_max += 1
        end
    end

    return n_max
end

# WARNING : changes phase of slice, do not use in phase calculations
function interpolate_slice(slice::IsoBECSlice; size = (2001, 2001))
    nx = size[1]
    ny = size[2]
    interpolated_psi = Array{Float64, 2}(undef, nx, ny)

    xs_old = LinRange(slice.x[begin], slice.x[end], slice.nx)
    ys_old = LinRange(slice.y[begin], slice.y[end], slice.ny)

    xs = LinRange(slice.x[begin], slice.x[end], nx)
    ys = LinRange(slice.y[begin], slice.y[end], ny)

    dx = (xs[end] - xs[begin]) / nx
    dy = (ys[end] - ys[begin]) / ny

    itp = cubic_spline_interpolation((xs_old, ys_old), abs.(slice.psi))

    for ix in eachindex(xs)
        for iy in eachindex(ys)
            interpolated_psi[ix, iy] = itp(xs[ix], ys[iy])
        end
    end

    return IsoBECSlice(
        interpolated_psi,
        xs,
        ys,
        nx,
        ny,
        dx,
        dy
    )
end

function get_itp_context(context::IsoBECContext)

    xs_old = LinRange(context.x[begin], context.x[end], context.nx)
    ys_old = LinRange(context.y[begin], context.y[end], context.ny)
    zs_old = LinRange(context.z[begin], context.z[end], context.nz)

    itp = cubic_spline_interpolation((xs_old, ys_old, zs_old), abs.(context.psi))

    return itp
end

function interpolate_context(context::IsoBECContext; size = (X_INTERPOLATION_SIZE, Y_INTERPOLATION_SIZE, Z_INTERPOLATION_SIZE))
    itp = get_itp_context(context)

    nx = size[1]
    ny = size[2]
    nz = size[3]

    interpolated_psi = Array{Float64, 3}(undef, nx, ny, nz)

    xs = LinRange(context.x[begin], context.x[end], nx)
    ys = LinRange(context.y[begin], context.y[end], ny)
    zs = LinRange(context.z[begin], context.z[end], nz)

    dx = (xs[end] - xs[begin]) / nx
    dy = (ys[end] - ys[begin]) / ny
    dz = (zs[end] - zs[begin]) / nz

    for ix in eachindex(xs)
        for iy in eachindex(ys)
            for iz in eachindex(zs)
                interpolated_psi[ix, iy, iz] = itp(xs[ix], ys[iy], zs[iz])
            end
        end
    end

    return IsoBECContext(
        interpolated_psi,
        xs,
        ys,
        zs,
        nx,
        ny,
        nz,
        dx,
        dy,
        dz
    )
end

function get_BEC_heights_rms(psi::IsoBECContext, max)
    # Calculate standard deviation of the BEC in z-direction
    z_rho = abs.(psi.psi[max.i, max.j, :])
    z_norm = z_rho ./ sum(z_rho)
    z_positions = LinRange(psi.z[begin], psi.z[end], psi.nz)
    mean_z = sum(z_norm .* z_positions)
    rms_width = sqrt(sum(z_norm .* (z_positions .- mean_z).^2))
    return 2 * rms_width  # ~2σ covers ~95% of distribution
end

function get_BEC_heights(psi::IsoBECContext, n_atoms; FWXM = 0.1)
    slice  = get_BEC_slice(psi)
    maxima = find_local_maxima(slice)
    coords = get_coordinates(slice, maxima)
    n_becs = length(maxima)

    heights = Vector{BECHeight}(undef, n_becs)
    for idx in eachindex(maxima)
        max     = maxima[idx]
        coord   = coords[idx]
        
        z_rho   = abs.(psi.psi[max.i, max.j, :])
        z_iter  = LinRange(psi.z[begin], psi.z[end], psi.nz)
        z_new   = LinRange(psi.z[begin], psi.z[end], 2000)

        itp = cubic_spline_interpolation(z_iter, z_rho)
        rho_itp = [itp(z) for z in z_new]
        max_z   = maximum(rho_itp)

        #cz_rho  = rho_itp .> (max_z * FWXM)
        #
        atoms_density = abs.(rho_itp).^2 * n_atoms * 1000 / length_au3_to_μm3(1.)
        cz_rho  = atoms_density .> (1)

        begin_idx = findfirst(cz_rho)
        end_idx = findlast(cz_rho)

        height  = z_new[end_idx] - z_new[begin_idx]

        heights[idx] = BECHeight(max.i, max.j, coord.x, coord.y, height, coord.value)
    end

    return heights
end

function get_BEC_maxrho(context::IsoBECContext)
    z_0_index_offset = findmax(abs.(context.psi))[2][3] - div(context.nz, 2)

    _slice = get_BEC_slice(context; z_0_index_offset = z_0_index_offset)

    slice = interpolate_slice(_slice)
    maxima = find_local_maxima(slice)

    coords = get_coordinates(slice, maxima)
    n_becs = length(maxima)

    maxrhos = Vector{Float64}(undef, n_becs)

    for (idx, max) in enumerate(maxima) 
        max_value = max.value^2 
        maxrhos[idx] = max_value / length_au3_to_μm3(1.)
    end

    return maxrhos
end

Base.show(io::IO, context::IsoBECContext) = print(io,
                                                  "BEC data with size ($(context.nx), $(context.ny), $(context.nz)) and spacings ($(context.dx) x $(context.dy) x $(context.dz))")
Base.show(io::IO, context::PoteContext) = print(io,
                                                  "Potential data with size ($(context.nx), $(context.ny), $(context.nz)) and spacings ($(context.dx) x $(context.dy) x $(context.dz))")
Base.show(io::IO, context::IsoBECSlice) = print(io, "BEC slice with size ($(context.nx), $(context.ny))")
