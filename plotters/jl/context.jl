using Interpolations

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

struct IsoBECSlice
    psi::Array{ComplexF64,2}
    x::Vector{Float64}
    y::Vector{Float64}
    nx::Int32
    ny::Int32
    dx::Float64
    dy::Float64
end

function get_BEC_slice(context::IsoBECContext)
    @assert context.nz % 2 == 1 "Wrong grid in 'z' dimension"
    z_0_index = Int64(floor(context.nz / 2) + 1)
    
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
    threshold_val = 0.3 * max_val
    
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
    threshold = 0.5
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

# WARNING : changes phase of slice, do not use in calculations
function interpolate_slice(slice::IsoBECSlice; size = (5000, 5000))
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

Base.show(io::IO, context::IsoBECContext) = print(io, "BEC data with size ($(context.nx), $(context.ny), $(context.nz))")
Base.show(io::IO, context::IsoBECSlice) = print(io, "BEC slice with size ($(context.nx), $(context.ny))")
