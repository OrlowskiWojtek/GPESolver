struct CutContext
    data::Array{Float64,2}  # 2D array for the cut data
    x::Vector{Float64}
    y::Vector{Float64}
    nx::Int
    ny::Int
    dx::Float64
    dy::Float64
end

struct IsoBCEContext
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
