struct CutContext
    data::Array{Float64,2}  # 2D array for the cut data
    x::Vector{Float64}
    y::Vector{Float64}
    nx::Int
    ny::Int
    dx::Float64
    dy::Float64
end

struct IsoBECContext
    psi::Array{ComplexF64,3}  # 3D array for the wavefunction data
    x::Vector{Float64}
    y::Vector{Float64}
    z::Vector{Float64}
    nx::Int
    ny::Int
    nz::Int
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

Base.show(io::IO, context::IsoBECContext) = print(io, "BEC data with size ($(context.nx), $(context.ny), $(context.nz))")
